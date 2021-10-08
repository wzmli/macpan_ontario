library(gsheet)

get_vacdat <- function(url, end_date){
  vacdat <- data.frame(gsheet2tbl(url))
  
  vacdat <- (vacdat
    ## get cumulative counts
    %>% transmute(date = as.Date(data...date)
                 , administered = data...total_vaccinations
                 , second = data...total_vaccinated
                 , first = administered - second
   )
   ## drop empty obs
    %>% filter(!is.na(administered),
               administered != 0)
   ## assume NAs for cummulative first and second doses at the beginning for the time series are there because all administered were first doses; replace NAs accordingly
   %>% mutate(
     second = case_when(
       is.na(first) ~ 0,
       T ~ second
     ),
     first = case_when(
     is.na(first) ~ administered,
     T ~ first
   ))
   ## compute daily rates
    %>% transmute(date
                 , vax_doses_per_day = diff(c(0,administered))
                 , vax_prop_first_dose = case_when(
                   vax_doses_per_day != 0 ~ diff(c(0,first))/vax_doses_per_day,
                   T ~ 1
   ))
    %>% filter(date <= end_date)
    %>% pivot_longer(-date)
    %>% rename(Date = date,
               Symbol = name,
               Relative_value = value)
  )
  
  return(vacdat)
}