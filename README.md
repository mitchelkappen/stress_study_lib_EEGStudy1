# Acoustic speech features in social comparison: how stress impacts the way you sound

This repository withholds the code base for the publication mentioned above.

## Project structure

```txt
├── loc_data          <── data generated / used by the scripts
├── notebooks         <── Python notebooks for data processing + exploratory analysis
│   └── speech_study  <── helper functions for the notebooks
├── reports           <── images / additional information related to the project.
└── scripts           <── R scripts for data analysis
```

Data processing is done in Python and corresponding notebooks can be found in subdirectory ```notebooks/```.
Subsequent analyses and visualizations are done in R and corresponding scripts are found in subdirectory ```scripts/```.
All notebooks and scripts are numbered in chronological order and run out of the box using relative paths.
Therefore, we advise to download this entire repository to reproduce our results exactly as we performed them.

## Python scripts

[Poetry](https://python-poetry.org/) is used as Python package management system.
To reconstruct the same environment in which the packages were installed:

1. Make sure you have poetry installed (follow
[these instructions](https://python-poetry.org/docs/#installation) if not )

2. Create the environment with the following command
```sh
poetry install
```
3. Activate the environment
```sh
poetry shell
```

## R scripts
R installation and package version numbers are available in the SupplementalMaterial.pdf file in the subdirectory ```Supplemental Material/```.

## Data availability
Raw audio files are only available upon request due to privacy implications.
All processed data is freely available in subdirectory ```loc_data/```.

## Notes

The following terms:
* `control` & `neutral` feedback
* `stress` & `negative` feedback

Are used interchangeably throughout the code.

<br>

---

<p align="center">
👤 <i> Mitchel Kappen, Jonas Van Der Donckt</i>
</p>
