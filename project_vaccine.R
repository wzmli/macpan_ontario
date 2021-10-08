## get vacdat

library(jsonlite)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)
library(lubridate)

end_date <- as.Date("2021-09-30")

get_data <- function(province){
	url <- paste0('https://api.covid19tracker.ca/reports/province/'
								, province
								, '?fill_dates=true'
	)
	upd <- fromJSON(url)
	dat <- (upd$data
					%>% mutate(province = province)
					%>% select(province, date,total_vaccinations,total_vaccinated)
	)
	return(dat)
}

provinces <- c("bc","ab","sk","mb","on","qc","nb","ns","pe","nl","yt","nt","nu")

listdat <- lapply(provinces,get_data)

print(listdat)


## clean dat


start_date <- as.Date("2020-12-01")
on_pop <- 14.57e6
eli_pop <- 12.93e6

## NAs in cumulative columns should all be zeroes
dd <- (bind_rows(listdat)
			 %>% mutate(NULL
			 					 , total_vaccinated = ifelse(is.na(total_vaccinated), 0, total_vaccinated)
			 					 , total_vaccinations = ifelse(is.na(total_vaccinations), 0, total_vaccinations)
			 )
)

## Sensible names and simple calcs
dd <- (dd
			 %>% group_by(province)
			 %>% transmute(date = as.Date(date)
			 							, secondVaxPop = total_vaccinated
			 							, vaxPop = total_vaccinations - secondVaxPop
			 							, dailyJabs = diff(c(0,total_vaccinations))
			 							, dailySecond = diff(c(0,total_vaccinated))
			 							, dailySecondProp = dailySecond/dailyJabs
			 )
			 %>% filter(date > start_date)
			 %>% filter(date < end_date)
			 %>% ungroup()
)

print(dd)

## vaccine plot 

popdat <- read_csv("params/pop.csv")

## Start with just ON
start_date <- as.Date("2020-12-01")

## Calculating daily dose prop and first dose

## Looking at the second dose prop for the last couple of days and the average vaccine administered

print(ggplot(dd, aes(date, dailySecondProp))
			+ geom_point()
			+ facet_wrap(~province,scale="free")
)

cum_dat <- (dd
						%>% select(date, vaxPop, secondVaxPop,province)
						%>% group_by(date,province)
						%>% pivot_longer(names_to = "population", values_to = "count",-c(date,province))
						%>% left_join(.,popdat)
)

canada_cum <- (cum_dat
							 %>% group_by(date, population)
							 %>% summarise(count = sum(count,na.rm=TRUE)     
							 							, pop = sum(pop,na.rm=TRUE)
							 							, eli_pop = sum(eli_pop,na.rm=TRUE)
							 							, province = "Canada")
)

print(cum_dat)
print(canada_cum)

cum_dat <- rbind(cum_dat, canada_cum)
cum_dat <- (
	filter(cum_dat, province %in% c("bc","ab","sk","mb","on","qc"))
	# filter(cum_dat, province %in% c("Canada","bc","ab","sk","mb","on","qc","nb","ns"))
	%>% mutate(province2 = ifelse(province == "Canada","Canada",province2))
)


## Cumulative first and second dose with population and eligible population

print(cum_dat)

gg_pop <- (ggplot(cum_dat,aes(date,y=count,color=population))
					 + geom_line()
					 + facet_wrap(~province2, scale="free")
					 + scale_color_manual(values=c("blue","black"))
					 + ylab("Cumulative count")
					 + geom_hline(aes(yintercept = pop))
					 + geom_hline(aes(yintercept = eli_pop),color="red")
					 + theme(legend.position="bottom")
					 + ggtitle("Black = pop; red = eligible pop")
)

print(gg_pop)


gg_pop2 <- (ggplot(cum_dat,aes(date,y=count/pop,color=population))
						+ geom_line()
						+ facet_wrap(~province2, scale="free")
						+ scale_color_manual(values=c("blue","black"))
						+ ylab("Cumulative count")
						+ geom_hline(aes(yintercept = pop/pop))
						+ geom_hline(aes(yintercept = eli_pop/pop),color="red")
						+ theme(legend.position="bottom")
						+ ggtitle("Black = pop; red = eligible pop")
)

print(gg_pop2)

## project

cutDays <- 7
futureSteps <- 100
satScale <- 0.5 ## Look into this; should it be different for doses or provinces?
hesitancy = c(0.1, 0.15)
hesitancy[[2]] <- 1 - (1-hesitancy[[2]])/(1-hesitancy[[1]])
print(hesitancy)

## Simple saturating model

startdat <- (dd
						 %>% group_by(province)
						 %>% filter(date >= (max(date)-7))
						 %>% mutate(across(.fns=mean)
						 					 , dailySecondProp = ifelse(is.na(dailySecondProp),dailySecond/dailyJabs,dailySecondProp)
						 					 #		, date = max(dd$date)
						 )
						 %>% distinct()
						 %>% left_join(popdat)
)



print(startdat)



canstart <- (startdat
						 %>% group_by(date)
						 %>% summarise(secondVaxPop = sum(secondVaxPop,na.rm=TRUE)
						 							, vaxPop = sum(vaxPop, na.rm=TRUE)
						 							, dailyJabs = sum(dailyJabs,na.rm=TRUE)
						 							, dailySecond = sum(dailySecond,na.rm=TRUE)
						 							, dailySecondProp = mean(dailySecondProp,na.rm=TRUE)
						 							, pop = sum(pop)
						 							, eli_pop = sum(eli_pop)
						 							, province = "Canada"
						 							, province2 = "Canada"
						 )
)

startdat <- rbind(startdat,canstart)


vfun <- function(vpop, steps, start, tpop, scale=1){
	if(length(tpop)==1){
		tpop = rep(tpop, steps)
	}
	stopifnot(length(tpop)==steps)
	stopifnot(tpop[[1]] > vpop + start*scale)
	maxvacc = (tpop[[1]]-vpop)*start/(tpop[[1]]-vpop-start*scale) 
	v <- numeric(steps)
	v[[1]] <- vpop
	for(i in 2:steps){
		pool <- tpop[[i]] - vpop
		vacc <- pool*maxvacc/(pool+maxvacc*scale)
		v[[i]] <- vpop <- vpop + vacc
	}
	return(v)
}

## The most naive saturating approach

provinces <- c("Canada","bc","ab","sk","mb","on","qc","nb","ns","nl","nt","yt","pe","nu")
provinces <- provinces[2:7]
# provinces <- "mb"



vacproject <- function(pp){
	
	start <- startdat %>% filter(province == pp) 
	
	vacc_project <- tibble(NULL
												 , date = seq(start$date, length.out=futureSteps, by=1)
												 , province = pp
												 , vaxPop=vfun(
												 	start$vaxPop, futureSteps, start$dailyJabs - start$dailySecond
												 	, ((1-hesitancy[[1]])*start$eli_pop)
												 	, satScale
												 )
												 , secondVaxPop=vfun(
												 	start$secondVaxPop, futureSteps, start$dailySecond
												 	, ((1-hesitancy[[2]])*vaxPop)
												 	, satScale
												 )
												 , firstJabs = diff(c(NA, vaxPop))
												 , secondJabs = diff(c(NA, secondVaxPop))
												 , jabs = firstJabs + secondJabs
												 , eli_pop = start$eli_pop
												 , pop = start$pop
												 
	)
	
	return(vacc_project)
	
}

projectList <- lapply(provinces,vacproject)
# projectList <- lapply(provinces[6],vacproject)

projectdat <- (bind_rows(projectList))


print(projectdat)

## project plot and save

projectdat <- (projectdat
							 %>% filter(province %in% c("bc","ab","sk","mb","on","qc","nb","ns","nl","Canada"))
							 %>% mutate( province2 = province
							 						, province2 = ifelse(province2 == "bc", "BC", province2)
							 						, province2 = ifelse(province2 == "ab", "AB", province2)
							 						, province2 = ifelse(province2 == "sk", "SK", province2)
							 						, province2 = ifelse(province2 == "mb", "MB", province2)
							 						, province2 = ifelse(province2 == "on", "ON", province2)
							 						, province2 = ifelse(province2 == "qc", "QC", province2)
							 						, province2 = ifelse(province2 == "nb", "NB", province2)
							 						, province2 = ifelse(province2 == "nl", "NL", province2)
							 						, province2 = ifelse(province2 == "ns", "NS", province2)
							 						, province2 = ifelse(province2 == "Canada", "Canada", province2)
							 						
							 						, province2 = factor(province2,level=c("BC","AB","SK","MB","ON","QC","NB","NS","NL","Canada"))
							 )
							 
)

gg_proj <- (gg_pop
						+ geom_line(data=projectdat, aes(x=date,y=vaxPop), color="black", lty="dashed")
						+ geom_line(data=projectdat, aes(x=date,y=secondVaxPop), color="blue", lty="dashed")
						+ geom_vline(aes(xintercept=as.Date("2021-10-04")))
						+ scale_x_date(date_labels="%b", date_breaks="1 month")
						+ facet_wrap(~province2, scale="free",ncol=2)
						+ ggtitle("")
)

print(gg_proj)
ggsave("vac_proj.png",width=10,height=8)


print(gg_proj2 <- (gg_pop2
									 + geom_line(data=projectdat, aes(x=date,y=vaxPop/pop), color="black", lty="dashed")
									 + geom_line(data=projectdat, aes(x=date,y=secondVaxPop/pop), color="blue", lty="dashed")
									 + geom_vline(aes(xintercept=as.Date("2021-10-04")))
									 + scale_x_date(date_labels="%b", date_breaks="1 month")
									 + facet_wrap(~province2, scale="free",ncol=2)
									 + ggtitle("")
									 + ylab("Proportion")
))

ggsave("vac_proj_prop.png",width=10,height=8)

spreadDat <- (cum_dat
							%>% spread(key=population,value=count)
)

print(spreadDat)

print(projectdat)

outputdat <- (projectdat
							#	%>% select(-c(firstJabs, secondJabs, jabs))
							%>% bind_rows(spreadDat)
							%>% arrange(province,date)
							%>% rename(at_least_1_dose = vaxPop
												 , double_dose = secondVaxPop
												 , population = pop
												 , eligible_population = eli_pop 
							)
							%>% mutate(type = ifelse(date <= as.Date("2021-10-04"), "data", "projection")
												 , unvaccinated_population = population - at_least_1_dose
												 , partially_dose = at_least_1_dose - double_dose
							)
)
print(outputdat)

write.csv(projectdat,"params/vaccine_projection_only.csv")

write.csv(outputdat,"vaccine_projection.csv")

