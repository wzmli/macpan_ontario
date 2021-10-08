time1 <- Sys.time()


## LOAD LIBRARIES
library(dplyr)
library(readr)
library(stringr)
library(gsheet)
library(tidyr)
library(zoo)
library(lubridate)
library(tidyr)
library(parallel)
library(McMasterPandemic)
source("get_vacdat.R")

options(macpan_pfun_method = "grep")

## LOAD DATA

## calibrate inputs
calibration_params_all <- read_csv("calibration_params_by_province.csv")

## load data to which we're fitting (various types of reports, e.g. infection reports, hospitalizations, etc.)

source("clean.R")

reports_all <- (all_sub
	%>% filter(var %in% c("report"))
)
								
								
## get province list

## GETTING PARAMETERS FROM THE PAST AS STARTING VALUES
province_list <- calibration_params_all$province
province_list <- c("ON")
calibration_params_all <- calibration_params_all %>% filter(province == province_list)

## setting initial values using past
# calibration_params_all$initial_relative_beta0 <- paste0(as.character(coef(past$fit,"fitted")[["time_params"]]),collapse = ";")

voc_start <- as.Date("2021-03-15")
est_date <- as.Date("2021-08-30")

voc_p0 <- 0.001
r_diff <- 0.07

# datevec <- seq.Date(as.Date("2021-03-15"),Sys.Date(),by=1)
# dd2 <- data.frame(date= datevec
# 	, prop = plogis(qlogis(voc_p0) + r_diff*0:(length(dd)-1))
# )



## FUNCTION TO CALIBRATE FOR A GIVEN PROVINCE
calibrate_province <- function(prov
	, calibration_params_all
	, reports_all
	, all_beta = TRUE
	, voc_start = voc_start
	, voc_p0 = voc_p0
	, r_diff = r_diff
	, est_date = est_date
	){
  
	# prov <- "ON"
  print(paste0("starting ", prov, " calibration..."))
  calib_time1 <- Sys.time()
  
  ## GET CALIBRATION PARAMS AND DATA TO FIT FOR THIS PROVINCE
  calibration_params <- calibration_params_all %>% filter(province == prov)
  report_vars <- unlist(strsplit(calibration_params[["report_vars"]],"/"))
  print(calibration_params)
  reports <- (reports_all 
    %>% filter(province == prov) 
    %>% select(-province)
    %>% filter(var %in% report_vars)
    %>% mutate(value = ifelse(value == 0, NA, value))
  )
  
  # SET UP TIME WINDOW
  start_date <- calibration_params[["trim_date"]]
  end_date <- as_date(max(reports$date))
  start_date_offset <- 60
  date_vec <- as_date(start_date:end_date)
  
  ## GET MODEL PARAMS
  model_params <- fix_pars(
    read_params(paste0("model_params_", prov,".csv")),
    target = c(R0 = calibration_params[["R0"]],
               Gbar = calibration_params[["Gbar"]])
  )
  model_params[["beta0"]] <- calibration_params$beta0 
  ## update beta0 based on different initial value for beta0 
  ## (used as initial condition for  calibration)

  ## add vaccination structure
  model_params <- expand_params_vax(
    model_params,
    model_type = "twodose",
    vax_doses_per_day = 1,
    vax_prop_first_dose = 1,
    vax_efficacy_dose1 = 0.6,
    vax_efficacy_dose2 = 0.9,
    vax_avg_response_time = 14,
    vax_avg_response_time_R = 7,
    vax_alpha_dose1 = 0.6,
    vax_alpha_dose2 = 0.9,
    vax_mu_dose1 = 1,
    vax_mu_dose2 = 1
  )
  ## ADD VARIANT STRUCTURE AFTER VAX STRUCTURE; USING VERY LOW VARIANT PROP, CAN'T SET IT TO ZERO/ONE LIKE VAXS
  model_params <- expand_params_variant(
  	model_params,
  	variant_prop = 1e-7,
  	variant_advantage = 1.5,
  	variant_vax_efficacy_dose1 = 0.3,
  	variant_vax_efficacy_dose2 = 0.8
  )
  
  
  print(summary(model_params))
  
  ## SET UP OPT_PARS AND TIME_ARGS
  
  ## get break dates
  break_dates <- as.Date(unlist(strsplit(
    calibration_params[["break_dates"]],
    split = ";")))
  bd <- break_dates
  if(!all_beta){
  	bd <- break_dates[break_dates>est_date]
  }
  
  
  n_break_dates = length(bd)
  
  ## get calibration starting points for beta0 relative values 
  ## (trim if truncating break dates)
  initial_relative_beta0 <- as.numeric(unlist(strsplit(
    calibration_params[["initial_relative_beta0"]],
    split = ";"
  )))
  
  
  # if(n_break_dates != length(initial_relative_beta0)) stop("number of break dates does not match number of initial conditions for calibration of beta0 relative values; please check specifications in params/calibration_params_by_province.csv")
  
  ## get vaccination data
  vacdat <- get_vacdat(calibration_params[["vacdat_url"]], end_date)
  vacdat <- (vacdat
  	%>% transmute(Date,Symbol,Value=Relative_value
  				, Type = "rel_orig"
  	)
  )
  
  ## SET UP CALIBRATION PARAMETERS
  
  ## set up all time-varying parameters 
  ## TODO, DO IT PROGRAMMATICALLY TO CHANGE ALL BETAS AFTER XX TO NA
  betadat <-data.frame(Date = break_dates
  	, Symbol = "beta0"
  	, Value = initial_relative_beta0
  	, Type = "rel_orig"
  )
  
  if(!all_beta){
  betadat <- (betadat
  	%>% mutate(Value = ifelse(Date >= est_date,NA,Value))						
  )
  }
  
  ## SET UP VOC TIME VARYING PARAMETERS 
  ## TODO: PROBABLY DO IT IN A R SCRIPT AND THEN SOURCE IT INSTEAD OF HARD CODING
  voc_date <- seq.Date(voc_start,end_date,by=1)
  vocdat <- data.frame(Date = voc_date
  	, Symbol = "variant_prop"
  	, Value = 1e7*plogis(qlogis(voc_p0) + r_diff*0:(length(voc_date)-1))
  	, Type = "rel_orig"
  )
  
  ## BIND ALL THE TIME VARYING DATA
  
  time_pars <- bind_rows(betadat,vocdat,vacdat)
  
  ## check whether first vaccination date is after simulation start date
  ## if so, shut off vaccination on first date of the sim
  first_vax_date <- min(vacdat$Date)
  if(start_date < first_vax_date){
    time_pars <- rbind(
      time_pars,
      data.frame(
        Date = start_date-start_date_offset, ## shut off vax initially
        Symbol = "vax_doses_per_day",
        Value = 0,
        Type = "rel_orig"
))
    
  }
  ## reorder by symbol, then date
  ## (i think calibrate/run_sim_break does this to ensure breaks are implemented
  ## in order, so this may not be necessary...)
  time_pars <- time_pars %>% arrange(Symbol, Date)
  
  ## set up calibration initial conditions
  ## IF ALL BETA, THEN CHANGE ALL RELATIVE BETA TO NA AND ADD BETA0
  
  
  opt_pars <- list(params=c(log_beta0 = log(model_params[["beta0"]]))
    , time_params = tail(initial_relative_beta0,length(bd))
  	)
  
  opt_pars <- list(time_params = tail(initial_relative_beta0,length(bd))
  )
  
  if(all_beta){
  	opt_pars <- list(params=c(log_beta0 = log(model_params[["beta0"]]))
  		, time_params = initial_relative_beta0
  		)
  }
  
  # opt_pars <- list(params=c(nonhosp_mort  = model_params[["rho"]]))
  # time_pars <- bind_rows(vocdat,vacdat)
  
  
  ## set up time-varying parameters for calibration
  ## (keep fixed ones and set ones that will be calibrated to NA)
  params_timevar <- time_pars
  if(all_beta){
  params_timevar$Value[which(params_timevar$Symbol == "beta0")] <- NA
  }
  ## GET CALIBRATION INPUT DATA
  
  ## format input data for calibration
  fitdat <- (data.frame(
    date = rep(date_vec,length(unique(reports[["var"]]))),
    var  = rep(report_vars ,each=length(date_vec))
  )
    %>% left_join(reports, by = c("date", "var"))
  )
  
  ## CALIBRATE
  res <- calibrate(
    base_params  = model_params
    , data       = fitdat
    , opt_pars   = opt_pars
    , sim_args   = list(ndt = 1,
                        step_args = list(do_hazard = TRUE))
    , time_args = list(params_timevar = params_timevar)
    , mle2_control = list(maxit = 1e3,
                          reltol = calibration_params[["reltol"]],
                          ## here beta and gamma are nelder-mead parameters,
                          ## not macpan model parameters!!!
                          beta = 0.5,
                          gamma = 2
                          )
    , start_date_offset = start_date_offset
    , use_DEoptim = FALSE
    , debug_plot = FALSE
    
  )
  
  calib_time2 <- Sys.time()
  
  ## SAVE CALIBRATION
  res_list <- list(fit = res,
                   calibration_params = calibration_params,
                   fitdat = fitdat,
                   vacdat = vacdat,
                   calib_time = print(calib_time2-calib_time1))
  
  filename <- paste0("calibration_",end_date,"_", prov, ".RDS")
  if(all_beta){
  	filename <- paste0("calibration_",end_date,"_", prov, "_full.RDS")
  }
  
  saveRDS(object=res_list, file=filename)
}

calibrate_province(province_list
	, calibration_params_all
	, reports_all
	, all_beta = FALSE
	, voc_start = voc_start 
	, voc_p0 = 0.001
	, r_diff = 0.07
	, est_date = est_date)

## calibrate each province in parallel
# mclapply(province_list,
#          calibrate_province,
#          calibration_params_all = calibration_params_all,
#          reports_all = reports_all,
#          mc.cores = 4)

time2 <- Sys.time()
print(time2-time1)
