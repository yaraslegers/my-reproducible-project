# Dynamic Linear Models for E. coli Monitoring in Broilers

## Description

In this project, I tried to fit a dynamic linear model to the incidence and prevalence of E. coli in necropsy findings from broilers sent in to Royal GD. The goal is to provide an early warning if the incidence increases in the Netherlands. (spoiler alert: it doesn't work).

## Usage

In the 'src' folder, there are three scripts. Ecoli.R is the main script, which sources from the other two. Side note: I only made the Ecoli.R myself.
This repository also contains a dummy dataset Necropsy.xlsx, so that you can run the script and produce graphs.

## Installation Requirements

This project is created in R version 4.0.3 (2020-10-10).

The following packages need to be installed:

-   TempPackage_1.0
-   docstring_1.0.0
-   lubridate_1.8.0
-   readxl_1.3.1
-   tidyr_1.2.0
-   dplyr_1.0.8

## Project organization

-   PG = project-generated
-   HW = human-writable
-   RO = read only

```{=html}
<!-- -->
```
    .
    ├── .gitignore
    ├── CITATION.md
    ├── LICENSE.md
    ├── README.md
    ├── requirements.txt
    ├── bin                <- Compiled and external code, ignored by git (PG)
    │   └── external       <- Any external source code, ignored by git (RO)
    ├── config             <- Configuration files (HW)
    ├── data               <- All project data, ignored by git
    │   ├── processed      <- The final, canonical data sets for modeling. (PG)
    │   ├── raw            <- The original, immutable data dump. (RO)
    │   └── temp           <- Intermediate data that has been transformed. (PG)
    ├── docs               <- Documentation notebook for users (HW)
    │   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
    │   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
    ├── results
    │   ├── figures        <- Figures for the manuscript or reports (PG)
    │   └── output         <- Other output for the manuscript or reports (PG)
    └── src                <- Source code for this project (HW)

## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)
