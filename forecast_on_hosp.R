time1 <- Sys.time()

library(tidyverse)
library(McMasterPandemic)
options(macpan_pfun_method = "grep")
source("forecast_newhelper.R")

end_date <- "2021-10-04"

n_sim <- 50
n_cores <- 6

flist <- "calibration_2021-10-04_ON.RDS"

#reopen_date is to tweak last parameter
#reopen_date2 is the actual the scenario projection date

lift_frame <- data.frame(province = c("ON")
												 , reopen_date = c("2021-09-02")
												 , reopen_date2 = c("2021-09-25")
												 , voc_start = c("2021-03-15")
												 , voc_p0 = c(0.001)
												 , r_diff = c(0.07)
												 , keep_last = c(1)
												 , scale_factor = c(1)
)

forecast_with_VoC <- function(x
															, forecast_params_all
															, n_days_forecast = 90
															, reopen_factor = 1.5
															, n_sim = 200
															, n_cores = 4
															, total_doses_factor = 1
															, reports_only = FALSE
															, use_new = FALSE
)
{
	## load calibration
	# x <- flist
	mod <- readRDS(file.path("results",x))
	
	reopen_date <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"reopen_date"] %>% as.Date()
	reopen_date2 <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"reopen_date2"] %>% as.Date()
	
	voc_start <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"voc_start"] %>% as.Date()
	voc_p0 <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"voc_p0"]
	r_diff <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"r_diff"]
	scale_factor <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"scale_factor"]
	keep_last <- (lift_frame %>% filter(province == mod$calibration_params$province))[,"keep_last"]
	
	
	## convert fit end date to date for calculations
	fit_end_date <- max(mod$fitdat$date)
	
	## get fitted calibration params and process close/reopen factors
	cc <- coef(mod$fit,"fitted")
	
	## get beta0 break dates
	beta0dat <- mod$fit$forecast_args$time_args$params_timevar %>% filter(Symbol == "beta0")
	beta0dat <- (beta0dat)

	mudat <-data.frame(Date = as.Date(c("2020-10-01","2020-11-01","2021-01-25","2021-03-01","2021-05-01","2021-07-30","2021-09-01"))
										 , Symbol = "mu"
										 , Value = c(0.99,0.99,0.995,0.985,0.99,0.965,0.99)
										 , Type = "rel_orig"
	)

	start_date <- mod$fit$forecast_args$start_date
	forecast_end_date <- fit_end_date + n_days_forecast
	
	## get time-varying vaccination params	
	vax_timevar <- mod$fit$forecast_args$time_args$params_timevar %>% filter(str_detect(Symbol, "^vax_"))
	
	vacproject <- read_csv("params/vaccine_projection_only.csv")
	vacproject_clean <- (vacproject
											 %>% mutate(province = toupper(province))
											 %>% filter(province == mod$calibration_params$province)
											 %>% filter(date > max(vax_timevar$Date))
											 %>% filter(date <= forecast_end_date)
	)
	
	vax_projection_jabs <- (vacproject_clean
													%>% transmute(Date = date
																				, Symbol = "vax_doses_per_day"
																				, Value = jabs
																				, Type = "rel_orig"
													)
													%>% filter(!is.na(Value))
	)
	
	vax_projection_prop <- (vacproject_clean
													%>% transmute(Date = date
																				, Symbol = "vax_prop_first_dose"
																				, Value = firstJabs/jabs
																				, Type = "rel_orig"
													)
													%>% filter(!is.na(Value))
	)
	vax_forecast_timevar <- rbind(vax_projection_prop,vax_projection_jabs)
	
	

	# ## add on vax projection for forecast period
	vax_timevar <- bind_rows(
		vax_timevar,
		vax_forecast_timevar
	)
	
	if(use_new){
		beta0_timevar <- (mod$fit$forecast_args$time_args$params_timevar
											%>% filter(Symbol == "beta0")
											# %>% mutate(Value = ifelse(Date == as.Date("2020-10-03"),0.88,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2020-11-23"),0.85,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2020-12-30"), 0.8928553 ,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2021-01-14"),0.65,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2021-06-11"),0.6,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2021-06-26"),0.9398743,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2021-07-16"),1.45,Value))
											# %>% mutate(Value = ifelse(Date == as.Date("2021-08-18"), 1.5 ,Value))

		)
		beta0_timevar <- rbind(beta0_timevar
													 , data.frame(Date = reopen_date
													 						 , Symbol = "beta0"
													 						 , Value = 1.05
													 						 , Type = "rel_prev"
													 )
													 , data.frame(Date = reopen_date2
													 						, Symbol = "beta0"
													 						, Value = reopen_factor
													 						, Type = "rel_prev"
													 )
		)
	}
	
	## get full parameters list
	# pp <- coef(mod$fit,"all")
	pp <- mod$fit$forecast_args$base_params
	
	
	voc_date <- seq.Date(voc_start,forecast_end_date,by=1)
	vocdat <- data.frame(Date = voc_date
											 , Symbol = "variant_prop"
											 , Value = 1e7*plogis(qlogis(voc_p0) + r_diff*0:(length(voc_date)-1))
											 , Type = "rel_orig"
	)
	
	
	## update forecast args
	fa <- mod$fit$forecast_args
	fa$start_date <- "2020-07-17"
	fa$end_date <- forecast_end_date
	fa$opt_pars <- mod$fit$forecast_args$opt_pars
	fa$base_params <- mod$fit$forecast_args$base_params
	fa$base_params["mu"] <- 0.98
	# fa$base_params["gamma_s"] <- 1/4
	fa$base_params["c_prop"] <- 0.33
	fa$base_params["beta0"] <- 0.502
	
	
	fa$time_args$params_timevar <- bind_rows(beta0_timevar,mudat,vocdat,vax_timevar)
	# aa <- bind_rows(beta0_timevar,vocdat,vax_timevar)
	
	## keep all state variables (split by vax cat) to inspect
	fa$sim_args$condense_args <- list(keep_all = TRUE)
	
	## check metadata in calibrate object before feeding into forecast_ensemble
	rr <- forecast_ensemble(mod$fit
													, nsim=n_sim 
													, forecast_args = fa
													, qvec = c(0.025,0.5,0.975)
													, scale_Sigma = scale_factor
													, seed = 1
													, parallel = TRUE
													, n_cores = n_cores
	)
	
	aa <- rr %>% filter(var %in% c("H"))
	rr2 <- (rr 
					%>% mutate(province = mod$calibration_params$province
										 , new_strain_fraction = voc_p0
										 , voc_start = voc_start
										 , reopen_factor = reopen_factor
										 , reopen_date = reopen_date
					)
					%>% left_join(., beta0_timevar, by=c("date"="Date"))
					%>% left_join(., pivot_wider(vax_timevar,
																			 id_cols = Date,
																			 names_from = "Symbol",
																			 values_from = "Value"), by = c("date" = "Date"))
	)
	
	if(reports_only){
		rr2 <- rr2 %>% filter(var %in% c("H"))
	}
	
	return(rr2)
}

sim1 <- lapply(flist,function(y){forecast_with_VoC(x=y,reopen_factor = 1, n_sim = n_sim, n_cores = n_cores,use_new=TRUE)})
sim2 <- lapply(flist,function(y){forecast_with_VoC(x=y,reopen_factor = 1.15, n_sim = n_sim, n_cores = n_cores,use_new=TRUE)})
sim3 <- lapply(flist,function(y){forecast_with_VoC(x=y,reopen_factor = 0.85, n_sim = n_sim, n_cores = n_cores,use_new=TRUE)})
betaforecast_dat <- bind_rows(sim1,sim2,sim3)


# betaforecast_dat <- bind_rows(sim1)

reportdat <- (betaforecast_dat %>% filter(var == "report"))
hospdat <- (betaforecast_dat
	%>% filter(grepl("H_",var) | grepl("H2_",var))
	%>% group_by(date,province,reopen_date,reopen_factor)						
	%>% summarise(lwr = sum(lwr)
				, value = sum(value)
				, upr = sum(upr)
				, var = "H"
	)
)

Xdat <- (betaforecast_dat
				 %>% filter(grepl("X_",var))
				 %>% group_by(date,province,reopen_date,reopen_factor)						
				 %>% summarise(lwr = sum(lwr)
				 							, value = sum(value)
				 							, upr = sum(upr)
				 							, var = "X"
				 )
				 %>% arrange(reopen_factor,date)
				 %>% group_by(reopen_factor)
				 %>% mutate(lwr = diff(c(0,lwr))
				 					 , value = diff(c(0,value))
				 					 , upr = diff(c(0,upr))
				 					 , var = "Hosp Admission"
				 )
)

combodat <- bind_rows(reportdat,hospdat,Xdat)
# betaforecast_dat <- bind_rows(sim1)
# betaforecast_dat <- bind_rows(sim3)


use_local_data_repo <- FALSE
source("scripts/clean.R")

betaforecast_dat2 <- (all_sub
											%>% transmute(date, var, province,obs = value)
											%>% left_join(combodat,.)
											%>% filter(date >= as.Date("2020-09-15"))
											# %>% filter(var %in% c("report","H"))
											%>% mutate(scale_factor = reopen_factor)
											# %>% filter(var %in% c("report","death","ICU","H"))
											
											# , obs = ifelse(date > as.Date("2021-02-21"),NA,obs)
)
vac <- (betaforecast_dat2
	%>% select(province,date,var,lwr,med=value,upr,obs,scale_factor)
)

write.csv(vac,"2021-10-04_ON_forecast.csv")


gg <- (ggplot(vac,aes(x=date))
	+ geom_line(aes(y=med))
	+ geom_line(aes(y=lwr),lty="dashed")
	+ geom_line(aes(y=upr),lty="dashed")
	+ geom_point(aes(y=obs))
	+ facet_grid(scale_factor~var,scale="free")
	# + facet_wrap(~var,scale="free")
	+ xlim(as.Date(c("2021-07-01","2021-10-30")))
	+ ylim(c(0,1000))
	+ theme_bw()
)

print(gg)

