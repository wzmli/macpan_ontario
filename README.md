# MacPan

## Current model

Accounts for vaccination, ``invading voc'', and piecewise break dates for reopening and public health measures


__WARNING__: Before running, make sure you have the latest version of McMasterPandemic package. https://github.com/mac-theobio/McMasterPandemic.git


### Calibration

- Need to gather latest information on all (6) PTs and implement to the model. The specific location is in the `params` subdirectory: `calibration_params_by_province.csv` and all the macpan parameters for each province are in their own separate csv files.

- run `scripts/calibate_xx.R`. This will automatically download the data off `wzmli/COVID-19_CANADA/clean.Rout.csv` the up to date reports/hospitalization data. 

- There are 3 options (need to tweak and refine a bit more) for calibrating the betas (transmission parameters):
  1. last/recent `rel_beta` (current)
  2. last/recent `rel_beta` with `beta0` (this will have CIs for the past as well)
  3. all betas

- currently skipping calibrating to hospitalization parameters and doing a posthoc fitting manually.

### Forecasting

- first run `scripts/project_vaccine.R` to make the vaccine projections
- this will automatically pull in the latest up to date vaccine by province data and project it and save the projections into the `params` folder
- run `scripts/forecast_xx_hosp.R` (**Warning**: _change the dates_ in each forecast script)
- run `scripts/combo_hosp_forecast.R` (remember to change the dates)

__WARNINGS__

- make sure your computer has enough cores to run in parallel
- check if the projections make sense (may need to tweak)






