# Hutson_etal_2025_LeSMicrophys

author: Abby Hutson, Cooperative Institute for Great Lakes Research, University of Michigan
contact: hutsona at umich dot edu

This repository contains scripts and namelists used to initialize and post-process the model simulations found in Hutson et al. 2025.  

Hutson, A., C. Pettersen, A. Fujisaki-Manome, G. E. Mann, S. W. Nesbitt, and M. S. Kulie, 2025: Evaluating the Performance of Microphysics Schemes against Observations during High-Impact Lake-Effect Snow Events. Wea. Forecasting, 40, 2273–2292, https://doi.org/10.1175/WAF-D-24-0249.1.

ABSTRACT:
In the Great Lakes region (GLR), lake-effect snow (LeS) events are a common occurrence, in which narrow, intense bands of convection cause snowfall downwind of the lakes. The shallow convection associated with LeS is dynamically different from deeper, synoptically driven snow, and the particle size distributions (PSDs) of the precipitation have different shapes, as well. This work considers whether or not the Thompson–Eidhammer microphysics scheme, which includes single-moment prediction of snow, is accurate in estimating the PSDs of LeS convection. The High-Resolution Rapid Refresh (HRRR) configuration of the Weather Research and Forecasting (WRF) Model is used to simulate three different LeS events in the GLR using two different microphysics schemes: the Thompson–Eidhammer “aerosol-aware” scheme and the Morrison double-moment scheme. Model-estimated PSDs are calculated and compared to observed PSDs at three locations in the region: Marquette, Michigan; Gaylord, Michigan; and Buffalo, New York. Model-predicted liquid water equivalent snowfall and snow density are also compared to observed products. It is found that parameterization performance varies depending on location, with Thompson struggling to create the correct PSD shape for Marquette. Both microphysics schemes do not perform well in predicting particles greater than 6 mm in diameter except in Buffalo, where both simulated and observed PSDs contain snow particles greater than 10 mm in diameter.

FILE CONTENTS:
 - namelist.input.HRRR: the namelist used to initialize the WRF model used to simulate three lake-effect snow cases. The WRF model is configured to replicate the High Resolution Rapid Refresh (HRRR) model, so we use WRF version 3.9.

 - les_functions.py: this file contains python code for functions used repeatedly throughout the study. That includes the functions used to recreate the particle size distribution shapes as predicted by the Thompson and Morrison schemes.
 -  observed_psd.ipynb: this iPython notebook was used to create the model-observation PSD comparisons, along with other aggregations and tests. While this notebook is to serve as an example for the applications of the functions defined in les_functions.py, reader beware! It is not well commented, and used solely for research purposes. Please reach out to the email listed above if you have questions or concerns.

 -  scaling_morrison_reflectivity.ipynb: this notebook contains the methodology used to find the fit needed to map the model output of radar reflectivity from the Morrison microphysics scheme onto the appropriate range of levels that are typically output using the Thompson scheme. See Hutson et al. 2025 for further details on why this was needed. 
