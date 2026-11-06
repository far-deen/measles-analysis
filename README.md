# Measles Analysis, Bangladesh

Exploratory analysis of measles surveillance data from the Directorate General of Health Services (DGHS), Bangladesh, performed in R.

## Data

- **Source:** Directorate General of Health Services (DGHS), Bangladesh. https://dghs.gov.bd/views/latest-news
- **Accessed:** 22 May 2026
- **Period covered:** 15 April 2026 to 15 May 2026
- **Method:** Data were manually extracted from the DGHS website into Excel (.xlsx) files. No further cleaning was required.
- **Limitations:** Manual extraction may contain transcription errors. Please compare against the original source before relying on any figure.

All data remain the property of DGHS. This repository only redistributes the extracted tables for reproducibility.

## Repository contents

- `data/` extracted DGHS data
- `code/` R scripts used for the analysis
- `images/` figures produced by the analysis

## Requirements

- R 4.6.0
- Packages:
  - Data handling: `tidyverse`, `readxl`, `rio`, `here`, `lubridate`, `aweek`
  - Epidemiology: `incidence2`, `i2extras`
  - Statistics: `MASS`, `fitdistrplus`
  - Spatial: `sf`, `spdep`, `geodata`
  - Plotting: `patchwork`, `cowplot`, `scales`, `ggrepel`, `viridis`, `RColorBrewer`
  - Helper: `pacman`

Install all packages with:

    install.packages(c("tidyverse", "readxl", "rio", "here", "lubridate",
                       "aweek", "incidence2", "i2extras", "MASS",
                       "fitdistrplus", "sf", "spdep", "geodata",
                       "patchwork", "cowplot", "scales", "ggrepel",
                       "viridis", "RColorBrewer", "pacman"))

Note: `geodata` may download map boundary files when the scripts run, so an internet connection could be needed.

## How to run

1. Download this repository (green **Code** button, then **Download ZIP**) and unzip it.
2. Open R or RStudio and set the working directory to the unzipped folder.
3. Run the scripts in `code/` in order.

## Session information

Analysis was run on R 4.6.0 (2026-04-24 ucrt), Windows 11 x64.

## License

- Code: MIT License (see `LICENSE`)
- Figures in `images/`: CC BY 4.0
- Data: property of DGHS, see source above

## Author

Nawaz Hossain Fardin ([nawazfardin3-art](https://github.com/nawazfardin3-art))
