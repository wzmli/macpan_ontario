suppressMessages({
    library(McMasterPandemic)
    library(readr)
    library(dplyr)
    library(tidyr)
})

use_local_data_repo <- FALSE

if(!exists("use_public")){
    use_public <- FALSE
}
if(!use_public){
    
    
    url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/git_push/clean.Rout.csv"
    
    if(use_local_data_repo) 
	url <- "../COVID19-Canada/git_push/clean.Rout.csv"

    dd <- read_csv(url)
    
    all <- (dd
            # %>% filter(Province=="ON")
            %>% select(Province,Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations,bestTotal)
            %>% group_by(Province)
            %>% mutate(newDeaths=c(NA,diff(deceased))
                       ## ON hosp includes ICU, our model compartment is just acute care
                       # , Hospitalization=ifelse(Province == "ON",Hospitalization,Hospitalization-ICU)
                       , newtotal = bestTotal
                       , newTests = diff(c(NA,newtotal))
                       , newConfirmations = ifelse((newTests == 0) & (newConfirmations == 0), NA, newConfirmations)
            )
            # %>% mutate(newConfirmations = ifelse((Province == "BC") & (weekdays(Date)%in%c("Monday","Sunday")), NA, newConfirmations))
            %>% select(-c(deceased,newTests,newtotal))
            %>% pivot_longer(names_to="var",-c(Date,Province))
            %>% setNames(tolower(names(.)))
            %>% ungroup()
    )
    
    ## translate variable names to internally used values
    ## drop unused variables
    keep_vars <- c("H","ICU","death","report","bestTotal")
    
    all_sub <- (all
                %>% mutate_at("var",trans_state_vars)
                %>% filter(var %in% keep_vars)
    )
    
    #print(all_sub,n=Inf)
}

if(use_public){
    province_dat <- data.frame(province_name=c("British Columbia", "Alberta"
                                               ,"Ontario","Quebec"
                                               , "Saskatchewan","Manitoba","New Brunswick", "Newfoundland and Labrador"
                                               , "Nova Scotia", "Prince Edward Island", "Northwest Territories", "Nunavut","Yukon"
    )
    , province = c("BC","AB","ON","QC","SK","MB","NB","NL","NS","PEI","NT","NU","YU")
    )
    
    public_url <- "https://health-infobase.canada.ca/src/data/covidLive/covid19-download.csv"
    public_dat <- read_csv(public_url)
    
    public_clean <- (public_dat
                     %>% select(province_name = prname, date
                                , cumreport=numtotal, cumdeath=numdeaths,newTests=numtestedtoday)	
                     %>% arrange(province_name,date)
                     %>% mutate(report = diff(c(0,cumreport))
                                , death = diff(c(0,cumdeath))
                     )
                     %>% select(-c(cumreport,cumdeath))
                     %>% gather(key="var",value="value",-date,-province_name)
    )
    
    all_sub <- (province_dat
                %>% left_join(.,public_clean)
    )
}


