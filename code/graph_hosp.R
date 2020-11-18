
neff = length(scenv)
nscen = length(scens_lbl)

t_deb  = round((per_g[1]-1)/dt/12) + 2 # number of time steps period 1
t_fin  = round(per_g[2]/dt/12) + 1  # number of time steps period 1
days = seq(td,td+t_fin-t_deb,1)
col = c(2:(N+1))
colc = c(2:(N+2))
incid = 1e6


# all compartments
ctry_par = read.csv(paste0(scenF,ctryref,".csv"), stringsAsFactors = F)

ctryF = paste0(ctry_par[1,3],"/")
obsF = read.csv(paste0(inputF,ctryF,ctry_par[1,2],".csv"), stringsAsFactors = F)
obsF$date = as.Date(obsF$date)

ctryNPI = read.csv(paste0(inputF,ctry_par[10,3],"/",ctry_par[10,2],".csv"), stringsAsFactors = F)
thresh_g[1] = ctryNPI$value[3]*incid

state <-read.csv(paste0(resF,ctryF,scenref,file = "_state.csv"))
modeltime<-startyear - dt +state[,1]

state <- state[t_deb:t_fin,]
pop_tot <- rowSums(state) / incid

# hospit
hosp <- read.csv(file = paste0(resF,ctryF,scenref,"_hospitalisations.csv"))
hosp <-data.frame(cbind(modeltime,hosp))
colnames(hosp) = c("time",catname)
hosp <- hosp[t_deb:t_fin,]
hosp_tot <-data.frame(scenario = scens_lbl[1],cbind(modeltime = modeltime[t_deb:t_fin],value=rowSums(hosp[,col])/pop_tot))
tmp = list(hosp_tot)

# 1st scenario 
nb = 1
while (nb <= length(scens[1,]))
{
  # hospit
  eff = c()
  for (i in 1:neff)
  {
    ef <- read.csv(file = paste0(resF,ctryF,scens[i,nb],"_hospitalisations.csv"))
    eff = cbind(eff,rowSums(ef[t_deb:t_fin,1:N])/pop_tot)
  }
  
  high = apply(eff,1,max)
  low = apply(eff,1,min)
  med = apply(eff,1,mean)

  hosp_tot <-data.frame(scenario = scens_lbl[nb+1],cbind(modeltime = modeltime[t_deb:t_fin],high=high,low=low,med=med))
  tmp[[nb+1]] = hosp_tot
  nb = nb +1 
}


##########################################################################
# plot
##########################################################################


# hospitalisations - all population age groups

line_col = c("darkblue","darkgreen")
area_col = c("tomato","lightblue")

g_tmp <- ggplot() +
  theme_classic() +
  scale_x_date(labels = date_format("%m/%y"),date_breaks = datebreak) +
  theme(legend.position= "bottom") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))  +
  scale_fill_manual(values=area_col) +
  labs(x = "Time",  y = "hospitalisations per 1M", title = "Hospitalisations", col = element_blank()) 


for(i in 2:nscen)
{
  g1 <- data.frame(tmp[[i]]) %>% mutate(modeltime = as.Date(td + (modeltime-startyear) * (tf-td)))
  g_tmp <- g_tmp +
  geom_ribbon(data=g1, aes(x = modeltime , ymin = low, ymax = high, fill = scenario), inherit.aes = FALSE) + 
  geom_line(data=g1, aes(x = modeltime , y = med),size = 0.5, inherit.aes = FALSE,col = "white") 
}
g1 <- data.frame(tmp[[1]]) %>% mutate(modeltime = as.Date(td + (modeltime-startyear) * (tf-td)))
g_tmp <- g_tmp + geom_line(data=g1, aes(x = modeltime , y = value),size = 1, inherit.aes = FALSE,col = "darkblue") +
   geom_point(aes(x=obsF$date,y=obsF$hosp_inc),shape = 23 , fill = "grey50" ) +
   geom_line(aes(x = g1$modeltime , y = thresh_g[1]),size = 0.5, linetype = "dashed", color = "black") 
  
ggsave(file=paste0(figF,ctryF,scname,file = "_graph_hosp.tiff"), plot=g_tmp, width=10, height=6, device ='tiff', scale = 0.75)
