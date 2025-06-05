## *Management strategy evaluation for salmon habitat restoration and hatchery supplementation *

#### Mark H. Sorel, Richard W. Zabel, Andrew R. Murdoch, and Sarah J. Converse

##### Please contact Mark Sorel for questions about the code or data: [mark.sorel\@dfw.wa.gov](mailto:mark.sorel@dfw.wa.gov){.email}

##### Secondary contact: Sarah Converse ([sconver\@usgs.gov](mailto:sconver@usgs.gov){.email})

------------------------------------------------------------------------

## Abstract

A primary challenge in conservation is the evaluation of management actions. For species with complex life histories, this evaluation is particularly challenging, and quantitative models are needed to understand how effects on particular life stages translate to objectives such as an acceptably low extinction risk or large abundance. We modeled the effects of increasing juvenile survival (as would result from habitat restoration) and hatchery supplementation on the achievement of recovery objectives for endangered Chinook salmon (*Oncorhynchus tshawytscha*) in the Wenatchee River Basin of Washington State. We evaluated these actions using simulations with a model that accounted for density-dependence operating in juvenile rearing habitat. Simulations indicated that there were compounding benefits of increasing survival in both natal streams, where juveniles are born, and in downstream areas where a portion of juveniles rear. We also found that the marginal benefit of producing more hatchery fish declined as juvenile survival improved. Simulations suggested that achieving recovery targets for abundance of spring-run Chinook salmon in the Wenatchee River would require more than a 150% increase in juvenile survival. Therefore, comprehensive management focusing on multiple factors affecting survival may be required to meet recovery criteria. recovery criteria.

### Table of Contents

### [Scripts](./scripts)

Contains scripts to run all analyses.

**IPM_4\_3_hatch.R** Constructs the integrated population model from an Rdata object `data/dat_IPM_hatch.Rdata` with data and a TMB model `src/IPM_non_centered_hatchery_scenarios.cpp`. Then posterior samples are drawn using tmbstan.

**new_plots.R** Conducts population projection simulations assuming different habitat and hatchery management strategies. The outputs from the simulations are saved in the results folder. The script also summarizes the results of the simulations in figures.

### [Src](./src)

**IPM_non_centered_hatchery_scenarios.cpp** The integrated population model written in TMB.

### [Data](./data)

Contains processed data and model inputs. For raw data and to see the steps for generating model inputs, see [this repository](https://github.com/Quantitative-Conservation-Lab/Sorel_etal_2023_CJFAS).

**dat_IPM_hatch.Rdata** Data inputs to the integrated population model.

**par_IPM_hatch.Rdata** Initial parameter inputs to the model.

**map_IPM_hatch.Rdata** Parameters to be treated as random effects.

**rand_par_IPM_hatch.Rdata** Parameters to be fixed at initial values during fitting

**broodstock_remova.xlsx** Historical numbers of natural origin fish removed for hatchery broodstock.

**proj_arrays_hatch_8_11_2022.Rdata** Simulated future trajectories of environmental variables included as covariates in the population .

### [Results](./results)

Contains raw and processed results

**ipm_fit_non_centered_5_02_hatch_scen.Rdata** The TMB integrated population model object. This is used to conducted the population projections using posterior samples of parameters, simulated environmental variables, and random draws of random effects of year.

**list_of_draws.rda** This file contains posterior samples from the integrated population model.


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

Sorel MH, RW Zabel, AR Murdoch, and SJ Converse. In prep. Management modeling of salmon habitat restoration and hatchery supplementation.

### How to Use this Repository

The `new_plots.R` file located in the 'scripts' folder is the main scripts for the analysis and can be used to generate the results in the paper. The IPM_4_3_hatch.R file in th 'scripts' folder draws posterior samples using tmbstan, but is very computationally and time intensive, so samples have been saved in the results and archived on [Zenodo](https://zenodo.org/records/10526151).
