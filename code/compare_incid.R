
t_deb  = round((per_g[1]-1)/dt/12) + 2 # number of time steps period 1
t_fin  = round(per_g[2]/dt/12) + 1  # number of time steps period 1
days = seq(td,td+t_fin-t_deb,1)
col = c(2:(N+1))
colc = c(2:(N+2))
scens = paste0(sctry,scenv)
thresh_g<-rep(NA,nscen)
# ************************************************************************************* #
#Plot outputs
# [[1]] Infections
# [[2]] Symptomatics cases
# [[3]] Model compartments
# [[4]] transition rate
# [[5]] births
# [[6]] death all causes
# [[7]] outbreak
# [[8]] Death covid
# [[9]] Hospitalisation COVID

# ************************************************************************************* #

#### RUN MODELS AND WRITE CSV FILES

nb = 1

# all compartments
ctry_par = read.csv(paste0(scenF,sctry[1],".csv"), stringsAsFactors = F)

ctryF = paste0(ctry_par[1,3],"/")
obsF = read.csv(paste0(inputF,ctryF,ctry_par[1,2],".csv"), stringsAsFactors = F)
obsF$date = as.Date(obsF$date)

ctryNPI = read.csv(paste0(inputF,ctry_par[10,3],"/",ctry_par[10,2],".csv"), stringsAsFactors = F)
thresh_g[1] = ctryNPI$value[3]*incid

state <-read.csv(paste0(resF,ctryF,scens[1],file = "_state.csv"))
modeltime<-startyear - dt +state[,1]

state <- state[t_deb:t_fin,]
susc_scen <- data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(state[,1+c(varind[1,],varind[7,])])))
pop_tot <- rowSums(state) / incid
susc_scen$value = susc_scen$value /pop_tot

pop_age<-NULL
for (n in 1:NC){
  pop_age[n]<-sum(state[,varind[,n]])
}
pop_age <- pop_age / nrow(state) /incid

# infections
inf <- read.csv(file = paste0(resF,ctryF,scens[1],"_infections.csv"))
inf <-data.frame(cbind(modeltime,inf))
colnames(inf) = c("time",catname)

inf <- inf[t_deb:t_fin,]
inf_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(inf[,col])/pop_tot))
inf_age_scen <-data.frame(scenario = scens_lbl[1],  outcome = "infection", age = catname, value = colSums(inf[,colc])/pop_age,stringsAsFactors = F)


# cases
cases <- read.csv(file = paste0(resF,ctryF,scens[1],"_cases.csv"))
cases <-data.frame(cbind(modeltime,cases))
colnames(cases) = c("time",catname)
cases <- cases[t_deb:t_fin,]
cases_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(cases[,col])/pop_tot))
cases_age_scen <-data.frame(scenario = scens_lbl[1],  outcome = "cases", age = catname, value = colSums(cases[,colc])/pop_age,stringsAsFactors = F)


# hospit
hosp <- read.csv(file = paste0(resF,ctryF,scens[1],"_hospitalisations.csv"))
hosp <-data.frame(cbind(modeltime,hosp))
colnames(hosp) = c("time",catname)
hosp <- hosp[t_deb:t_fin,]
hosp_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(hosp[,col])/pop_tot))
hosp_age_scen <-data.frame(scenario = scens_lbl[1], outcome = "hospitalisation", age = catname, value = colSums(hosp[,colc])/pop_age,stringsAsFactors = F)

# covid deaths
deathC <- read.csv(file = paste0(resF,ctryF,scens[1],"_death COVID.csv"))
deathC <-data.frame(cbind(modeltime,deathC))
colnames(deathC) = c("time",catname)
deathC <- deathC[t_deb:t_fin,]
deathC_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(deathC[,col])/pop_tot))
deathC_age_scen <-data.frame(scenario = scens_lbl[1],outcome = "death", age = catname, value = colSums(deathC[,colc])/pop_age,stringsAsFactors = F)

# number vaccinated
nvac <- read.csv(file = paste0(resF,ctryF,scens[1],"_vaccinated.csv"))
nvac <-data.frame(cbind(modeltime,nvac))
colnames(nvac) = c("time",catname)
nvac <- nvac[t_deb:t_fin,]
nvac_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(nvac[,col])/pop_tot))

# distancing
distage <- read.csv(file = paste0(resF,ctryF,scens[1],"_distancing.csv"))
distage <-data.frame(cbind(modeltime,distage))
colnames(distage) = c("time",catname)
distage <- distage[t_deb:t_fin,]
# wdist = apply(distage[,col], 1, function(x) weighted.mean(x, pop_age[1:N]))
dist_scen <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=distage[,2]))

nb = 2
while (nb <= length(scens))
{
  ctry_par = read.csv(paste0(scenF,sctry[nb],".csv"), stringsAsFactors = F)
  ctryF = paste0(ctry_par[1,3],"/")
  ctryNPI = read.csv(paste0(inputF,ctry_par[10,3],"/",ctry_par[10,2],".csv"), stringsAsFactors = F)
  thresh_g[nb] = ctryNPI$value[3]*incid
  
  state <-read.csv(paste0(resF,ctryF,scens[nb],file = "_state.csv"))
  modeltime<-startyear - dt +state[,1]
  
  state <- state[t_deb:t_fin,]
  susc_tot <- data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(state[,1+c(varind[1,],varind[7,])])))
  susc_tot$value = susc_tot$value /pop_tot
  susc_scen = rbind(susc_scen,susc_tot)
  
  
  # infections
  inf <- read.csv(file = paste0(resF,ctryF,scens[nb],"_infections.csv"))
  inf <-data.frame(cbind(modeltime,inf))
  colnames(inf) = c("time",catname)
  
  inf <- inf[t_deb:t_fin,]
  inf_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(inf[,col])/pop_tot))
  inf_scen = rbind(inf_scen,inf_tot)
  inf_age <-data.frame(scenario = scens_lbl[nb],  outcome = "infection", age = catname, value = colSums(inf[,colc])/pop_age,stringsAsFactors = F)
  inf_age_scen = rbind(inf_age_scen,inf_age)
  
  # cases
  cases <- read.csv(file = paste0(resF,ctryF,scens[nb],"_cases.csv"))
  cases <-data.frame(cbind(modeltime,cases))
  colnames(cases) = c("time",catname)
  cases <- cases[t_deb:t_fin,]
  cases_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(cases[,col])/pop_tot))
  cases_scen = rbind(cases_scen,cases_tot)
  cases_age <-data.frame(scenario = scens_lbl[nb],  outcome = "cases", age = catname, value = colSums(cases[,colc])/pop_age,stringsAsFactors = F)
  cases_age_scen = rbind(cases_age_scen,cases_age)
  
  
  # hospit
  hosp <- read.csv(file = paste0(resF,ctryF,scens[nb],"_hospitalisations.csv"))
  hosp <-data.frame(cbind(modeltime,hosp))
  colnames(hosp) = c("time",catname)
  hosp <- hosp[t_deb:t_fin,]
  hosp_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(hosp[,col])/pop_tot))
  hosp_scen = rbind(hosp_scen,hosp_tot)
  hosp_age <-data.frame(scenario = scens_lbl[nb], outcome = "hospitalisation", age = catname, value = colSums(hosp[,colc])/pop_age,stringsAsFactors = F)
  hosp_age_scen = rbind(hosp_age_scen,hosp_age)
  
  # covid deaths
  deathC <- read.csv(file = paste0(resF,ctryF,scens[nb],"_death COVID.csv"))
  deathC <-data.frame(cbind(modeltime,deathC))
  colnames(deathC) = c("time",catname)
  deathC <- deathC[t_deb:t_fin,]
  deathC_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(deathC[,col])/pop_tot))
  deathC_scen = rbind(deathC_scen,deathC_tot)
  deathC_age <-data.frame(scenario = scens_lbl[nb],outcome = "death", age = catname, value = colSums(deathC[,colc])/pop_age,stringsAsFactors = F)
  deathC_age_scen = rbind(deathC_age_scen,deathC_age)
  
  # number vaccinated
  nvac <- read.csv(file = paste0(resF,ctryF,scens[nb],"_vaccinated.csv"))
  nvac <-data.frame(cbind(modeltime,nvac))
  colnames(nvac) = c("time",catname)
  nvac <- nvac[t_deb:t_fin,]
  nvac_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(nvac[,col])/pop_tot))
  nvac_scen = rbind(nvac_scen,nvac_tot)
  
  # distancing
  distage <- read.csv(file = paste0(resF,ctryF,scens[nb],"_distancing.csv"))
  distage <-data.frame(cbind(modeltime,distage))
  colnames(distage) = c("time",catname)
  distage <- distage[t_deb:t_fin,]
  # wdist = apply(distage[,col], 1, function(x) weighted.mean(x, pop_age[1:N]))
  dist_tot <-data.frame(scenario = scens_lbl[nb],cbind(modeltime = modeltime[t_deb:t_fin],value=distage[,2]))
  dist_scen = rbind(dist_scen,dist_tot)
  
  nb = nb +1 
}


##########################################################################
# plot
##########################################################################


  # hospitalisations - all population age groups
  
  tmp <- hosp_scen %>% mutate(modeltime = as.Date(td + (modeltime-startyear) * (tf-td)))

  g_tmp <- ggplot() +
    theme_classic() +
    geom_line(aes(x = tmp$modeltime , y = tmp$value, col = tmp$scenario),size = 1.5) +
    # geom_point(aes(obsF=hosp2$date,y=hosp2$mhosp),shape = 23 , fill = "grey50" ) +
    scale_color_manual(values=cb1) +
    labs(x = "Time",  y = "hospitalisations per 1M", title = "Hospitalisations", col = element_blank()) +
    scale_x_date(labels = date_format("%m/%y"),date_breaks = datebreak) +
    theme(legend.position= "bottom") +
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  if (thresh_line){
    tmp <- tmp  %>% mutate(thresh = NA)
    for (i in 1:nscen){
      tmp$thresh[tmp$scenario== scens_lbl[i]] <- thresh_g[i]
    }
    g_tmp <- g_tmp +
      geom_line(aes(x = tmp$modeltime , y = tmp$thresh, group  = tmp$scenario),size = 0.5, linetype = "dashed", color = "black") 
  }
  
  ggsave(file=paste0(figF,ctryF,scname,file = "_hosp_inc.png"), plot=g_tmp, width=10.5, height=7.5, device ='png', scale = 0.75)

  #  cases
  
  tmp <- cases_scen %>% mutate(modeltime = as.Date(td + (modeltime-startyear) * (tf-td)))
  
  g_tmp <- ggplot() +
    theme_classic() +
    geom_line(aes(x = tmp$modeltime , y = tmp$value, col = tmp$scenario),size = 1.5) +
    geom_point(aes(x=obsF$date,y=obsF$cases_inc),shape = 23 , fill = "grey50" ) +
    labs(x = "Time",  y = "cases per 1M", title = "Symptomatic cases", col = element_blank()) +
    scale_x_date(labels = date_format("%m/%y"),date_breaks = datebreak) +
    scale_color_manual(values=cb1) +
    theme(legend.position= "bottom") 
  
  ggsave(file=paste0(figF,ctryF,scname,file = "_cases_inc.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)
  
  
  # Death population - all population age groups
  
  tmp <- deathC_scen %>% mutate(modeltime = as.Date(td + (modeltime-startyear) * (tf-td)))
  
  g_tmp <- ggplot() +
    theme_classic() +
    geom_line(aes(x = tmp$modeltime , y = tmp$value ,col = tmp$scenario),size=1.5) +
    geom_point(aes(x=obsF$date,y=obsF$deaths_inc),shape = 23 , fill = "grey50" ) +
    scale_color_manual(values=cb1) +
    labs(x = "Time",  y = "deats per 1M", title = "Deaths", col = element_blank()) +
    scale_x_date(labels = date_format("%m/%y"),date_breaks = datebreak) +
    theme(legend.position= "bottom") 
  ggsave(file=paste0(figF,ctryF,scname,file = "_death_inc.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)
  

