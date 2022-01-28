import warnings
from typing import Dict

import neurokit2 as nk
import numpy as np
import pandas as pd
from context_aware.interfaces.empatica import EMPATICA_CONFIG
from context_aware.modules.tsflex.chunking import chunk_data
# -------------------------------- TIME SERIES --------------------------------
from context_aware.modules.tsflex.processing import (
    dataframe_func,
    SeriesProcessor,
    SeriesPipeline,
)
from context_aware.modules.tsflex.processing.utils import process_chunks_multithreaded
from context_aware.time_series_tools.processing.gsr.clean_gsr import clean_eda


# ------------------------------ GSR PROCESSING ---------------------------------
def clean_empatica_eda(
        gsr_series: pd.Series,
        valid_window: int = 3,
        max_interpolate: int = 30,
        min_valid_len: int = 5 * 60,
        eda_sqi_smoothen_win_size_s: float = 4,
) -> pd.DataFrame:
    fs = EMPATICA_CONFIG.SENSOR_FREQUENCIES["gsr"]

    out = clean_eda(
        gsr_series.to_frame(),
        "EDA",
        fs=fs,
        valid_window=valid_window,
        max_interpolate=max_interpolate,
        min_valid_len=min_valid_len,
    )

    out["EDA_SQI_smoothened"] = (
            out["EDA_SQI"]
            .rolling(
                int(fs * eda_sqi_smoothen_win_size_s + 1),  # center=True
            )
            .sum()
            / (fs * eda_sqi_smoothen_win_size_s + 1)
    )
    # if SQI is 0 -> retain zero
    out["EDA_SQI_smoothened"] *= out["EDA_SQI"].map(int)
    return out


def decompose_eda_empatica(df_cleaned: pd.Series, method="cvxEDA") -> pd.DataFrame:
    df_cleaned_ = df_cleaned.dropna()
    if len(df_cleaned_):
        sanitized_eda: pd.Series = nk.signal_sanitize(df_cleaned_)
        return nk.eda_phasic(
            sanitized_eda,
            sampling_rate=EMPATICA_CONFIG.SENSOR_FREQUENCIES["gsr"],
            method=method,
        ).set_index(df_cleaned_.index)
    else:
        print(
            "\tcleaned shape: ",
            df_cleaned.shape,
            "\tcleaned dropna shape: ",
            df_cleaned_.shape,
            flush=True,
        )
        warnings.warn(
            f"no not-Nan data `decompose_eda_empatica` will return an "
            + "empty dataframe"
        )
        return pd.DataFrame(columns=["EDA_Phasic", "EDA_Tonic"])
    # TODO -> add slope features


def find_peaks_empatica(df_phasic: pd.Series, method="neurokit") -> pd.DataFrame:
    # note -> GAMBOA2008 it's peaks seem to be exactly the same as neurokit!
    df_phasic = df_phasic.dropna()
    if not len(df_phasic):
        warnings.warn("find_peaks_empatica -> `df_phasic` has no not-nan values!")
        return pd.DataFrame(
            columns=[
                "SCR_RiseTime",
                f"SCR_Peaks_{method}",
                "SCR_RecoveryTime",
                "SCR_Amplitude",
            ]
        )

    peak_signal, _ = nk.eda_peaks(
        df_phasic,
        sampling_rate=EMPATICA_CONFIG.SENSOR_FREQUENCIES["gsr"],
        amplitude_min=0,  # TODO -> play more with this parameter
        method=method,
    )
    return peak_signal.set_index(df_phasic.index).rename(
        columns={"SCR_Peaks": f"SCR_Peaks_{method}"}
    )


@dataframe_func
def remove_false_positives_acc(
        df,
        min_rise_time_s=1,
        min_recovery_time_s=1,
        max_rise_time_s=15,
        max_recovery_time_s=7,
        min_scr_amplitude=0.02,
        scr_amplitude_ratio=0.03,
        eda_sqi_smoothened_threshold=0.6,
        max_acc_std=45,
) -> pd.DataFrame:
    assert min_rise_time_s < max_rise_time_s
    assert min_recovery_time_s < max_recovery_time_s
    acc_cols = ["ACC_x", "ACC_y", "ACC_z"]
    scr_cols = [
        "SCR_RiseTime",
        "SCR_Peaks_neurokit",  # not that loosely coupled
        "SCR_RecoveryTime",
        "SCR_Amplitude",
        # additional EDA metadata
        "EDA_Tonic",
        "EDA_SQI_smoothened",
    ]

    df_scr = df[scr_cols].dropna(how="all").copy()
    df_scr = df_scr.fillna(0)
    if not len(df_scr):
        warnings.warn(
            "remove_false_positives_acc -> `df_scr` has no not-nan values!"
        )
        return pd.DataFrame(columns=scr_cols).add_suffix("_reduced_acc")

    df_acc = df[acc_cols].dropna()
    # calculate the std
    win_duration_s = 1
    win_size = win_duration_s * EMPATICA_CONFIG.SENSOR_FREQUENCIES["acc"]
    win_size += 1 if win_size % 2 == 0 else 0  # make sure if win_size is odd
    fs_gsr = 4
    df_acc_std = (
        pd.concat(
            [df_acc[col].rolling(win_size, center=True).std() for col in acc_cols],
            axis=1,
            keys=acc_cols,
        )
            .max(axis=1)
            .rename("ACC_std")
            .resample(f"{int(1000 // fs_gsr)}ms")
            .mean()
            .rolling(2 * fs_gsr + 1, center=True)
            .max()
            .dropna()
    )

    df_scr = pd.merge_asof(
        df_scr, df_acc_std, left_index=True, right_index=True, direction="nearest"
    )

    df_scr.loc[df_scr["SCR_RiseTime"] < min_rise_time_s, scr_cols] = 0
    df_scr.loc[df_scr["SCR_RiseTime"] > max_rise_time_s, scr_cols] = 0

    # df_scr.loc[df_scr["SCR_RecoveryTime"] < min_recovery_time_s, scr_cols] = 0
    df_scr.loc[df_scr["SCR_RecoveryTime"] > max_recovery_time_s, scr_cols] = 0

    # additional acc & amplitude related processing (compared to vic's processing)
    # 1. relative SCR_Amplitude processing
    df_scr.loc[
        df_scr["SCR_Amplitude"]
        < np.clip(
            scr_amplitude_ratio * df_scr["EDA_Tonic"],
            a_min=min_scr_amplitude,
            a_max=None,
        ),
        scr_cols,
    ] = 0

    # 2. Remove peaks detected near low EDA sqi data
    df_scr.loc[
        df_scr["EDA_SQI_smoothened"] < eda_sqi_smoothened_threshold, scr_cols
    ] = 0

    df_scr.loc[df_scr["ACC_std"] > max_acc_std, scr_cols] = 0
    # todo -> do not add _reduce_acc suffix to the 'add_std_col'
    return pd.concat(
        [
            df_scr.drop(
                columns=["EDA_Tonic", "ACC_std", "EDA_SQI_smoothened"]
            ).add_suffix("_reduced_acc"),
            df_scr["ACC_std"],
        ],
        axis=1,
        ignore_index=False,
    )


# -------------------------- The processing pipelines
gsr_processing_pipeline = SeriesPipeline(
    processors=[
        SeriesProcessor(series_names=["EDA"], function=clean_empatica_eda),
        SeriesProcessor(
            series_names=["EDA"],
            function=decompose_eda_empatica,
        ),
        SeriesProcessor(
            series_names=["EDA_Phasic"],
            function=find_peaks_empatica,
            # additional kwargs
            method="neurokit",
        ),
        SeriesProcessor(
            series_names=tuple(
                [
                    "SCR_RiseTime",
                    "SCR_Peaks_neurokit",
                    "SCR_RecoveryTime",
                    "SCR_Amplitude",
                    "EDA_Tonic",
                    "EDA_SQI_smoothened",
                    "ACC_x",
                    "ACC_y",
                    "ACC_z",
                ]
            ),
            function=remove_false_positives_acc,
        ),
    ]
)


# ------------- PIPELINES WRAPPERS
def process_gsr_pipeline(
        data_dict: Dict[str, pd.DataFrame], multiprocessing=True, show_progress=True
) -> pd.DataFrame:
    if multiprocessing:
        df_out_gsr = pd.concat(
            process_chunks_multithreaded(
                same_range_chunks_list=chunk_data(
                    # only process a sub-chunk of the data dict
                    data=[data_dict[k] for k in ["gsr", "acc"]],
                    fs_dict=EMPATICA_CONFIG.SIGNAL_FREQUENCIES,
                    chunk_range_margin_s=10,
                    min_chunk_dur_s=60 * 5,
                    max_chunk_dur_s=60 * 60,
                    sub_chunk_overlap_s=60 * 5,
                    verbose=False,
                    copy=False,
                ),
                series_pipeline=gsr_processing_pipeline,
                n_jobs=16,
                show_progress=show_progress,
                drop_keys=["ACC_x", "ACC_y", "ACC_z"],
                return_all_series=True,
                return_df=True,
            )
        )
    else:
        df_out_gsr = gsr_processing_pipeline.process(
            data=list(data_dict.values()),
            return_all_series=False,
            return_df=True,
            drop_keys=["ACC_x", "ACC_y", "ACC_z"],
        )

    df_out_gsr = (
        df_out_gsr.rename(columns={"EDA": "EDA_Processed"})
            .dropna(axis=0, subset=["EDA_SQI_smoothened"],
                    how="all")  # must be all, right?
            .reset_index()
            .drop_duplicates(subset="timestamp", keep="first")
            .set_index("timestamp", drop=True)
            .sort_index()
            # TODO -> fix `Na` investigation in `clean_eda`
            .fillna(method="backfill")
    )

    # note -> with real-not continuous data we will need to work on a gap-like basis ...
    # assert pd.infer_freq(df_out_gsr.index) is not None
    return df_out_gsr
