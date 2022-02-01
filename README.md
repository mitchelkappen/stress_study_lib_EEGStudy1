# Speech features in social comparison: how stress impacts the way you sound

This repository withholds the code base for the publication mentioned above.

## Project structure 

```txt
├── loc_data          <── data generated / used by the scripts
├── notebooks         <── Python notebooks for data processing + exploratory analysis
│   └── speech_study  <── helper functions for the notebooks
├── reports           <── images / additional information related to the project.
└── scripts           <── R scripts for data analysis
```

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