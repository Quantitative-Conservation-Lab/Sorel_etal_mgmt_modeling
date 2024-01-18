## *Management modeling of salmon habitat restoration and hatchery supplementation*

#### Mark H. Sorel, Richard W. Zabel, Andrew R. Murdoch, and Sarah J. Converse

##### Please contact Mark Sorel for questions about the code or data: [mark.sorel\@dfw.wa.gov](mailto:mark.sorel@dfw.wa.gov){.email}

##### Secondary contact: Sarah Converse ([sconver\@usgs.gov](mailto:sconver@usgs.gov){.email})

------------------------------------------------------------------------

## Abstract

A primary challenge in conservation is optimally allocating resources to achieve desired outcomes. There are a variety of management interventions available for the conservation of imperiled salmonids in North America, and tools are needed to inform the allocation of management effort to activities such as habitat restoration and hatchery supplementation. We convened a workshop of decision makers and experts on Chinook salmon (*Oncorhynchus tshawytscha*) conservation in the Wenatchee River Basin of Washington State. Participants helped define candidate strategies involving habitat restoration and hatchery management. We evaluated the effects of these strategies on a variety of metrics -- including persistence and maintenance of wild genetic structure -- using simulations built on an integrated population model. Simulations indicated that restoration of natal streams would result in greater natural production than restoration of downstream habitats. However, downstream restoration would benefit fish from multiple natal streams as well as fish that spawn in downstream habitats. Reducing hatchery broodstock sizes from current targets would result in relatively small changes in natural productivity and reduce the risk of the hatchery program affecting adaptations of the wild population. Using principals of decision analysis and management modeling enabled evaluation of alternative management strategies that are most relevant to conservation decisions.

### Table of Contents

### [Scripts](./scripts)

Contains scripts to run all analyses.

**IPM_4\_3_hatch.R** Constructs the integrated population model from an Rdata object `data/dat_IPM_hatch.Rdata` with data and a TMB model `src/IPM_non_centered_hatchery_scenarios.cpp`. Then posterior samples are drawn using tmbstan.

**PVA-preliminary-results-scenarios.Rmd** Conducts population projection simulations assuming different habitat and hatchery management strategies. The outputs from the simulations are saved in the results folder. The script also summarizes the results of the simulations in text and figures. This Rmarkdown script can be knitted to generate the figures in the manuscript.

### [Src](./src)

**IPM_non_centered_hatchery_scenarios.cpp** The integrated population model written in TMB.

### [Data](./data)

Contains processed data and model inputs. For raw data and to see the steps for generating model inputs, see [this repository](https://github.com/Quantitative-Conservation-Lab/Sorel_etal_2023_CJFAS).

**dat_IPM_hatch.Rdata** Data inputs to the integrated population model.

**par_IPM_hatch.Rdata** Initial parameter inputs to the model.

**map_IPM_hatch.Rdata** Parameters to be treated as random effects.

**rand_par_IPM_hatch.Rdata** Parameters to be fixed at initial values during fitting

**broodstock_remova.xlsx** Historical numbers of natural origin fish removed for hatchery broodstock.

**proj_arrays_hatch_8\_11_2022.Rdata** Simulated future trajectories of environmental variables included as covariates in the population .

### [Results](./results)

Contains raw and processed results

**ipm_fit_non_centered_4\_13_hatch_scen.Rdata** The TMB integrated population model object. This is used to conducted the population projections using posterior samples of parameters, simulated environmental variables, and random draws of random effects of year.

**list_of_draws.rda** This file contains posterior samples from the integrated population model.

**list_of_sims_8\_11_22.Rdata** Raw results of population projection simulations under alternative management strategies. This file is too big for GitHub but is available for download on [Zenodo](https://zenodo.org/records/10526151). The 'PVA-preliminary-results-scenarios.Rmd' file has code (line 477) to download the file from Zenodo into the correct place in the file directory.

**summarized_sims_01_17_24.Rdata** Processed and summarized results of population projection simulations.

**2_panel_map_elev_9242021.png** Map used in the manuscript.

### Required Packages and Versions Used

base_4.2.2\
dplyr_1.1.0\
forcats_1.0.0\
ggplot2_3.4.1\
graphics_4.2.2\
grDevices_4.2.2\
here_1.0.1\
lubridate_1.9.1\
methods_4.2.2\
purrr_1.0.1\
readxl_1.4.3\
stats_4.2.2\
stringr_1.5.0\
tibble_3.1.8\
tidyr_1.3.0\
tidyverse_2.0.0\
TMB_1.9.3\
TMBhelper_1.4.0\
tmbstan_1.0.7\
utils_4.2.2\
viridisLite_0.4.1

### Details of Article

Sorel MH, RW Zabel, AR Murdoch, and SJ Converse. In review. Management modeling of salmon habitat restoration and hatchery supplementation.

### How to Use this Repository

The Rmarkdown file located in the 'scripts' folder is the main scripts for the analysis and can be knit to generate the results in the paper. Several of the steps in that script of circumvented to save time, and results are simply read in from an Rdata file. However, those code chunks that are set to eval=FALSE, could be turned back on to re-generate the results. The IPM_4\_3_hatch.R file in th 'scripts' folder draws posterior samples using tmbstan, but is very computationally and time intensive, so samples have been saved in the results and archived on [Zenodo](https://zenodo.org/records/10526151).
