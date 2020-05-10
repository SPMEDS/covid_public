library(deSolve)
library(tidyverse)

#########################################
#
# Parameters
#
#########################################

S0 = 1
  
inputF = "input/"
resF = "results/"
obsF = "observed/"

startyear=2020 # starting year of simulation

dir.create(file.path(paste0("./",resF), scen, "/"), showWarnings = FALSE, recursive = TRUE)

cbPalette <- c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
td = as.Date("2020-01-01")

##########################################
#
#  Model functions  
#
############################################

N<-9   # number of age classes
B<-8   # number of variables per class
A<-31  # number of transitions per class
V<-N*B # total number of variables
L<-N*A # total number of transitions
yd <- 360 # number of days per year
md <-30  # number of days per month
dd <- 1 #  timestep in days
pd <- 5 # period for distancing and imported
dt<-dd/yd # timestep in years
tsteps<-tmonths*md # number of time steps
tyears <- tmonths/12
time<-startyear+seq(0,tyears,dt) # time vector

t_deb  = round((per_deb-1)/12/dt) + 2 # number of time steps period 1
t_fin  = round(per_fin/12/dt) + 1  # number of time steps period 1
days = seq(td,td+t_fin-t_deb,1)


EQ2<-function(L, N, oldeq, transit,trans)
{
  eq = oldeq
  for(i in 1:L){
    
    iu1 = trans[i,1];
    iu2 = trans[i,3];
    iv1 = trans[i,2];
    iv2 = trans[i,4];
    
    eq[iu1] = eq[iu1] + transit[i]*iv1;
    eq[iu2] = eq[iu2] + transit[i]*iv2;
    
  }
  eq
}


# ************************************************************************************* #
# define variables
# ************************************************************************************* #

# 1=S:  susceptible
# 2=E:  exposed
# 3=I1: infectious -  asymptomatic
# 4=I2: infectious -  mild
# 5=I3: infectious -  severe
# 6=R:  recovered 
# 7=W:  Susceptible to secondary infection 
# 8=WE: Exposed to secondary infection

# ************************************************************************************* #
# define indices
# ************************************************************************************* #
varind<-matrix(0,nrow=B,ncol=N)
traind<-matrix(0,nrow=A,ncol=N)
for (n in 1:N){
  for (b in 1:B){
    varind[b,n]<-(n-1)*B+b
  }
  for (a in 1:A){
    traind[a,n]<-(n-1)*A+a
  }
}

# ************************************************************************************* #
# Transitions
# ************************************************************************************* #

trans <-matrix(0,nrow=L,ncol=4)
for (n in 1:N){
  nxt=(n<N)*(n+1)+(n==N)*(N)
  trans[traind[1,n],]<-c(varind[1,n],0,varind[1,1],1) # birth 
  trans[traind[2,n],]<-c(varind[1,n], -1,varind[1,n],0) # death S=1 
  trans[traind[3,n],]<-c(varind[2,n], -1,varind[1,n],0) # death E 
  trans[traind[4,n],]<-c(varind[3,n], -1,varind[1,n],0) # death I1 
  trans[traind[5,n],]<-c(varind[4,n], -1,varind[1,n],0) # death I2 
  trans[traind[6,n],]<-c(varind[5,n], -1,varind[1,n],0) # death I3 
  trans[traind[7,n],]<-c(varind[6,n], -1,varind[1,n],0) # death R 
  trans[traind[8,n],]<-c(varind[7,n], -1,varind[1,n],0) # death W 
  trans[traind[9,n],]<-c(varind[8,n], -1,varind[1,n],0) # death WE 

  trans[traind[10,n],]<-c(varind[1,n], -1,varind[1,nxt],+1) # age next level S 
  trans[traind[11,n],]<-c(varind[2,n], -1,varind[2,nxt],+1) # age next level E 
  trans[traind[12,n],]<-c(varind[3,n], -1,varind[3,nxt],+1) # age next level I1 
  trans[traind[13,n],]<-c(varind[4,n], -1,varind[4,nxt],+1) # age next level I2 
  trans[traind[14,n],]<-c(varind[5,n], -1,varind[5,nxt],+1) # age next level I3 
  trans[traind[15,n],]<-c(varind[6,n], -1,varind[6,nxt],+1) # age next level R 
  trans[traind[16,n],]<-c(varind[7,n], -1,varind[7,nxt],+1) # age next level W 
  trans[traind[17,n],]<-c(varind[8,n], -1,varind[8,nxt],+1) # age next level WE 

  trans[traind[18,n],]<-c(varind[1,n], -1,varind[2,n],+1) # inc S -> E
  trans[traind[19,n],]<-c(varind[2,n], -1,varind[3,n],+1) # inc E -> I1
  trans[traind[20,n],]<-c(varind[2,n], -1,varind[4,n],+1) # inc E -> I2
  trans[traind[21,n],]<-c(varind[2,n], -1,varind[5,n],+1) # inc E -> I3
  trans[traind[22,n],]<-c(varind[3,n], -1,varind[6,n],+1) # recovery I1 -> R
  trans[traind[23,n],]<-c(varind[4,n], -1,varind[6,n],+1) # recovery I2 -> R
  trans[traind[24,n],]<-c(varind[5,n], -1,varind[6,n],+1) # recovery I3 -> R
  trans[traind[25,n],]<-c(varind[5,n], -1,varind[6,n],0) # specific deaths due to the disease for I3
  trans[traind[26,n],]<-c(varind[6,n], -1,varind[7,n],+1) # waning natural R -> W
  trans[traind[27,n],]<-c(varind[7,n], -1,varind[8,n],+1) # inc W -> WE
  trans[traind[28,n],]<-c(varind[8,n], -1,varind[3,n],+1) # inc WE -> I1
  trans[traind[29,n],]<-c(varind[8,n], -1,varind[4,n],+1) # inc WE -> I2
  trans[traind[30,n],]<-c(varind[8,n], -1,varind[5,n],+1) # inc WE -> I3
  
  trans[traind[31,n],]<-c(varind[1,n], 0,varind[1,n],0) # epidemic period
  
}

# ************************************************************************************* #
# Function to calculate transition rates, given variables and parameters
# ************************************************************************************* #
trates <- function(x, parmal, t,ti) {
  with(as.list(c( parmal)),
       {
         t_internal<-as.integer(round(t*yd)) + ti
         p_int = (floor(t_internal/pd)) %% (md/pd) + 1 
         m_int = min(floor(t_internal/md) + 1, max_month)
         y_int = min(floor(t_internal*pd/yd) + 1, max_yr)
         m_vac = 1+ (m_int-1) %% 12 
         
         #Set up matrices
         pop<-c(rep(0,N))    # population sizes  
         foif<-c(rep(0,N))     # forces of infection
         tranrate<-matrix(0,nrow=N,ncol=A)   # transition rate matrix
         
         #population counter
         for (n in 1:N){pop[n]<-sum(x[varind[,n]])}   #  human population per age class
         poptot = sum(x)
         
         # Social distancing
         outbreak = sum(x[varind[5,]]) 
         dist = distancing[m_int,p_int]
         
         # imported cases
         ext = imported[m_int,p_int]
         
         # Social distancing threshold
         
         if((m_int < thresh_start) |  (m_int > thresh_stop))
         {
           dist = distancing[m_int,p_int]
           if (dist < 1)
           {
             if ((thresh_per == 0) & (wave == 0))
             {
               thresh_day <<- t_internal
               wave  <<-  1
             }
             
             if((t_internal - thresh_day) < thresh_dur3) 
             {
               thresh_per <<- wave
             }
             else
             {
               thresh_per <<- wave
             }
             
           }
         }
         else if (outbreak > thresh)
         {
           ext = 0
           if (thresh_per < 2)
           {
             thresh_day <<- t_internal
             wave  <<- wave + 1
           }
           thresh_per <<- wave
           if ((t_internal - thresh_day) < thresh_dur1)
           {
             dist = min(dist,thresh_dist1)
           }
           else
           {
             dist = min(dist,thresh_dist2)
           }
         }
         else if((t_internal - thresh_day) < thresh_dur2) 
         {
           dist = min(dist,thresh_dist2)
           thresh_per <<- wave
         }
         else if((t_internal - thresh_day) < thresh_dur3) 
         {
           thresh_per <<- 0
           dist = min(dist,thresh_dist3)
         }
         else 
         {
           thresh_per <<- 0
         }
         
         #imported cases
         #ext = ifelse(dist<1,0,imported[m_int,p_int])
         
         # seasonality
         
         t_seas = t + ti/yd
         seas = 1 + seas_amp*cospi(2*(t_seas - seas_month/md))
         

         for (n in 1:N){
           
           #Force of infection
           foif[n]<-dist*seas*p_trans*sum(contact[,n]*(f1*x[varind[3,]]+f2*x[varind[4,]]+f3*x[varind[5,]]))/poptot + ext 
           
           mortrate<-mortmat[n,y_int]
           births<-ifelse(n==1,1,0) * birthmat[1,y_int]

           
           tranrate[n,]<-c(      
             births,   # rate of birth -> R0   1
             mortrate*x[varind[1,n]],       # rate of death of S   2
             mortrate*x[varind[2,n]],       # rate of death of E   3
             mortrate*x[varind[3,n]],       # rate of death of I1  4
             mortrate*x[varind[4,n]],       # rate of death of I2  5
             mortrate*x[varind[5,n]],       # rate of death of I3  6
             mortrate*x[varind[6,n]],       # rate of death of R   7
             mortrate*x[varind[7,n]],       # rate of death of W   8
             mortrate*x[varind[8,n]],       # rate of death of WE  9

             (1-mortrate)*agerate[n]*x[varind[1,n]],        # aging s              10
             (1-mortrate)*agerate[n]*x[varind[2,n]],        # aging E              11
             (1-mortrate)*agerate[n]*x[varind[3,n]],        # aging I1              12
             (1-mortrate)*agerate[n]*x[varind[4,n]],        # aging I2              13
             (1-mortrate)*agerate[n]*x[varind[5,n]],        # aging I3              14
             (1-mortrate)*agerate[n]*x[varind[6,n]],        # aging R              15
             (1-mortrate)*agerate[n]*x[varind[7,n]],        # aging W              16
             (1-mortrate)*agerate[n]*x[varind[8,n]],        # aging WE              17

             foif[n]*x[varind[1,n]],     # inc S -> E     18
             poiE[n,1]*epsilon*x[varind[2,n]],     # inc E -> I1   19
             poiE[n,2]*epsilon*x[varind[2,n]],     # inc E -> I2   20
             poiE[n,3]*epsilon*x[varind[2,n]],     # inc E -> I3   21
             delta*x[varind[3,n]],  # recovery I1 -> R     22
             delta*x[varind[4,n]],  # recovery I2 -> R     23
             (1-poiE[n,4])*delta*x[varind[5,n]],   # recovery I3 -> R     24
             poiE[n,4]*delta*x[varind[5,n]], # death I3      25
             tau*x[varind[6,n]],               # natural waning R -> W  26
             foif[n]*x[varind[7,n]],     # inc W -> WE     27
             poiWE[n,1]*epsilon*x[varind[8,n]],     # inc WE -> I1   28
             poiWE[n,2]*epsilon*x[varind[8,n]],     # inc WE -> I2   29
             poiWE[n,3]*epsilon*x[varind[8,n]],     # inc WE -> I3   30
             
             thresh_per  # outbreak 31
           )
         }
         return(c(t(tranrate)))
       })
}

# ************************************************************************************* #
# Post-processing function
# ************************************************************************************* #
postproc <- function(parpro,out,tran) {
  with(as.list(c( parpro)),
       { 
         # ************************************************************************************* #
         # for outputting the  time series for each patch
         # ************************************************************************************* #
         #Incidence       
         outbreak<-sev<-deathC<-deathN<-birth<-symp<-inf<-matrix(0,nrow=length(out[,1]),ncol=N)
         for (n in 1:N){
           inf[,n]<-rowSums(tran[,c(traind[c(c(19:21),c(28:30)),n])])*dt
           symp[,n]<-rowSums(tran[,c(traind[c(c(20:21),c(29:30)),n])])*dt
           birth[,n]<-tran[,c(traind[1,n])]*dt
           deathN[,n]<-rowSums(tran[,c(traind[c(2:9),n])])*dt
           outbreak[,n]<-tran[,c(traind[31,n])]
           deathC[,n]<-tran[,c(traind[25,n])]*dt
           sev[,n]<-rowSums(tran[,c(traind[c(21,30),n])])*dt
         }
         
         return(cbind(inf,    #1
                      symp,#2
                      birth,      #3
                      deathN,       #4
                      outbreak,           #5
                      deathC,           #6
                      sev     #7  
         ))
         
       })
}
ti<-0
# ************************************************************************************* #
# ************************************************************************************* #
# ************************************************************************************* #
# ODE SOLVER
# ************************************************************************************* #
# ************************************************************************************* #

epiModel<-function(t,state, parode) 
{
  with(as.list(c(state, parode)),
       {
         
         #   # ************************************************************************************* #
         #   # define variables
         #   # ************************************************************************************* #
         
         Z=state
         
         # rates of change
         ti<-0
         transit<-trates(Z[1:V],parode,Z[V+1],ti)  
         
         if (sum(is.na(transit))>0)  {
           stop("transit NA   ",Z[V+1], "                                      ", 
                as.data.frame(transit))
         }
         
         eq<-rep(0.0, V)
         
         #eq<-EQ(L, N, eq, transit,trsu1,trsu2,trsv1,trsv2)
         eq<-EQ2(L, N, eq, transit,trans)
         
         eq[V+1]<-1
         
         dZ<-eq
         
         # return the rate of change
         list(c(dZ))
       }
  ) 
  # end with(as.list ...
}

runmod<-function(parrun){
  initcondrun<-NULL
  for(n in 1:N){
    initcondrun<-c(initcondrun,c(
      S=S0*parrun$population[n], 
      E=0.0*parrun$population[n], 
      I1=0*parrun$population[n], 
      I2=0*parrun$population[n], 
      I3=0*parrun$population[n], 
      R=(1-S0)*parrun$population[n], 
      W=0*parrun$population[n], 
      WE=0*parrun$population[n] 
    ))
  }
  #initoderun<-round(initcondrun)
  initoderun<-initcondrun
  staterun <- c(initcondrun,0)
  ti<-0  
  
  transitrun <- trates(initoderun,parrun,0,ti)
  
  
  # SOLVE THE ODEs and get output
  timesrun <- seq(0, tyears, by = dt) # Model run time
  #Solve ODE
  outoderun <- ode(y = staterun, times = timesrun, func = epiModel, parms = parrun, method  = "euler", hini=0.01)
  # Compute transitions at each time step
  tranoderun<-matrix(0,nrow=length(outoderun[,1]),ncol=L)
  
  wave <<- 0
  thresh_day <<- - yd
  thresh_per <<- 0
  
  for (ti in 0:tsteps){
    tranoderun[ti+1,]<-t(trates(outoderun[ti+1,2:(1+V)],parrun,0,ti))
  }
  #Compute outputs
  ppoutrun<-postproc(parrun,outoderun,tranoderun)
  modeltimes<-outoderun[,1]+startyear
  wave <<- 0

  inf_ode<-ppoutrun[,1:N]
  symp_ode<-ppoutrun[,(N+1):(2*N)]
  birth_ode<-ppoutrun[,(2*N+1):(3*N)]
  deathN_ode<-ppoutrun[,(3*N+1):(4*N)]
  out_ode<-ppoutrun[,(4*N+1):(5*N)]
  deathC_ode<-ppoutrun[,(5*N+1):(6*N)]
  hosp_ode<-ppoutrun[,(6*N+1):(7*N)]
  
  return(list( inf_ode, #1
               symp_ode,    #2
               outoderun,    #3
               tranoderun,    #4
               birth_ode, #5
               deathN_ode,  #6
               out_ode,  #7               
               deathC_ode,  #8
               hosp_ode
  ) )
}

runscen<-function(scen) 
{
  
  # ************************************************************************************* #
  # import Scenario data
  # ************************************************************************************* #
  
  # read input files
  
  transmission = read.csv(paste0(inputF,scen,"/transmission.csv"))
  demographics = read.csv(paste0(inputF,scen,"/demographics.csv"))
  mortality = read.csv(paste0(inputF,scen,"/mortality.csv"))
  birth = read.csv(paste0(inputF,scen,"/birth.csv"))
  contact= read.csv(paste0(inputF,scen,"/contact.csv"))
  distancing = read.csv(paste0(inputF,scen,"/distancing.csv"))
  imported = read.csv(paste0(inputF,scen,"/imported.csv"))
  distancing_threshold = read.csv(paste0(inputF,scen,"/distancing_threshold.csv"))
  generic = read.csv(paste0(inputF,scen,"/generic.csv"))
  poi = read.csv(paste0(inputF,scen,"/poi.csv"))
  
  catname <<-demographics[,1]

  thresh_day <<- - yd
  thresh_per <<- 0
  wave <<-0

  pars=list(
    p_trans=transmission$value[1]*yd,
    f1=generic$value[1], 
    f2=generic$value[2], 
    f3=generic$value[3], 
    epsilon= yd/generic$value[4], 
    delta= yd/generic$value[5], 
    tau=1/generic$value[6], 
    seas_month = generic$value[7], 
    seas_amp = generic$value[8], 
    max_yr = dim(mortality)[2], 
    max_month = dim(distancing)[1], 
    thresh_start = distancing_threshold$value[1], 
    thresh_stop = distancing_threshold$value[2], 
    thresh = distancing_threshold$value[3], 
    thresh_dur1 = distancing_threshold$value[4], 
    thresh_dur2 = distancing_threshold$value[5], 
    thresh_dur3 = distancing_threshold$value[6], 
    thresh_dist1 = distancing_threshold$value[7], 
    thresh_dist2 = distancing_threshold$value[8], 
    thresh_dist3 = distancing_threshold$value[9], 
    poiE = poi[1:N,2:5], 
    poiWE = poi[(N+1):(2*N),2:5], 
    contact = contact ,
    distancing = distancing[,-1], 
    imported  = imported[,-1], 
    catname = catname ,#0k
    agerate = demographics[,2] , #0k
    population = demographics[,3] ,
    mortmat = mortality ,
    birthmat  = birth
    
  ) #ok
  

  # ************************************************************************************* #
  # Run model
  # ************************************************************************************* #
  
  st<-Sys.time()
  run1<-runmod(pars)
  ed<-Sys.time()
  ed-st
  return(run1)
}

#####################################################
#
#    scenarios
#
#####################################################



################ 
#  observed data
################    

obs = read.csv(paste0(obsF,fobs,".csv"), stringsAsFactors = F)
obs$date = as.Date(obs$date)
obsd  = obs[c(1,6)] %>% mutate(mcases = stats::filter(cases_base, rep(1/pd,pd), sides=2))
nobs <- dim(obs)[1]
obsd <- obsd[-1,]
obsd$mcases[1] <- 0
obsd = obsd[is.na(obsd$mcases) != TRUE,]


################ 
#  run scenario
################    

st<-Sys.time()
run1 <- runscen(scen) 
ed<-Sys.time()
ed-st



# ************************************************************************************* #
#Plot outputs
# [[1]] infections 
# [[2]] symptomatic cases
# [[3]] Model compartments
# [[4]] transition rate
# [[5]] births
# [[6]] death all causes
# [[7]] outbreak
# [[8]] Death covid
# [[9]] Hospitalisation COVID

# ************************************************************************************* #

modeltime<-startyear - dt +run1[[3]][,1]

timeA<-time 
col = c(2:(N+1))

#### CREATE TABLES

# infections
inf <-data.frame(cbind(modeltime,run1[[1]]))
colnames(inf) = c("time",catname)
write.csv(inf,file = paste0(resF,scen,"/infections.csv"), row.names = FALSE)

# cases
cases <-data.frame(cbind(modeltime,run1[[2]]))
colnames(cases) = c("time",catname)
write.csv(cases,file = paste0(resF,scen,"/cases.csv"), row.names = FALSE)

# hospit
hosp <-data.frame(cbind(modeltime,run1[[9]]))
colnames(hosp) = c("time",catname)
write.csv(hosp,file = paste0(resF,scen,"/hospitalisations.csv"), row.names = FALSE)

# covid deaths
deathC <-data.frame(cbind(modeltime,run1[[8]]))
colnames(deathC) = c("time",catname)
write.csv(deathC,file = paste0(resF,scen,"/deaths.csv"), row.names = FALSE)

# state compartments
state <-data.frame(run1[[3]])
write.csv(state,file = paste0(resF,scen,"/state.csv"), row.names = FALSE)

# outbreak
outb <- data.frame(run1[[7]])
write.csv(outb,file = paste0(resF,scen,"/outbreak.csv"), row.names = FALSE)


#### CREATE TABLES

# infections
inf <- inf[t_deb:t_fin,]
inf_tot <-data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],inf=rowSums(inf[,col])))
# cases
cases <- cases[t_deb:t_fin,]
cases_tot <-data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],cases=rowSums(cases[,col])))
# hospit
hosp <- hosp[t_deb:t_fin,]
hosp_tot <-data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],hosp=rowSums(hosp[,col])))

# covid deaths
deathC <- deathC[t_deb:t_fin,]
deathC_tot <-data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],deathC=rowSums(deathC[,col])))

# outbreak
outb <- outb[t_deb:t_fin,]
outb_tot <-data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],outb=outb[,1]))


# all compartments
state <- state[t_deb:t_fin,]
susc_states <- c(varind[1,],varind[7,])
susc_tot <- data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],susceptibles = rowSums(state[,1+c(varind[1,],varind[7,])])))
pop_tot <- data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],pop = rowSums(state)))
susc_tot$susceptibles = susc_tot$susceptibles /pop_tot$pop 

hospit <- data.frame(scenario = scen,cbind(modeltime = modeltime[t_deb:t_fin],hosp = rowSums(state[,1+varind[5,]])))


##########################################################################
# plot
##########################################################################


# infections - all population age groups

tmp <- inf_tot %>% mutate(modeltime = days)
g_tmp <- ggplot() +
  theme_classic() +
  geom_line(aes(x = tmp$modeltime , y = tmp$inf,col = tmp$scenario),size = 1.5) +
  labs(x = "Time",  y = "Number of events", title = "Infections",col= scen) +
  scale_x_date(labels = scales::date_format("%m/%y"),date_breaks = datebreak) +
  theme(legend.position= "bottom") 
ggsave(file=paste0(resF,scen,"/infections.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)

# cases - all population age groups

tmp <- cases_tot %>% mutate(modeltime = days)
g_tmp <- ggplot() +
  theme_classic() +
  geom_line(aes(x = tmp$modeltime , y = tmp$cases,col = tmp$scenario),size = 1.5) +
  geom_point(aes(x=obsd$date,y=obsd$mcases)) +
  labs(x = "Time",  y = "Number of events", title = "Symptomatic cases",col= scen) +
  scale_x_date(labels = scales::date_format("%m/%y"),date_breaks = datebreak) +
  theme(legend.position= "bottom") 
ggsave(file=paste0(resF,scen,"/cases.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)

# hospitalization  - all population age groups

tmp <- hosp_tot %>% mutate(modeltime = days)
g_tmp <- ggplot(tmp,aes(x = modeltime , y = hosp , col = scenario)) +
  theme_classic() +
  geom_line(size = 1.5) +
  labs(x = "Time",  y = "Number of events", title = "Hospitalisations",col= scen) +
  scale_x_date(labels = scales::date_format("%m/%y"),date_breaks = datebreak) +
  theme(legend.position= "bottom") 
g_tmp
ggsave(file=paste0(resF,scen,"/hospitalisations.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)

# deaths - all population age groups

tmp <- deathC_tot %>% mutate(modeltime = days)
g_tmp <- ggplot(tmp,aes(x = modeltime , y = deathC , col = scenario)) +
  theme_classic() +
  geom_line(size = 1.5) +
  labs(x = "Time",  y = "Number of events", title = "Deaths",col= scen) +
  scale_x_date(labels = scales::date_format("%m/%y"),date_breaks = datebreak) +
  theme(legend.position= "bottom") 
ggsave(file=paste0(resF,scen,"/deaths.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)

# Susceptibles - all population age groups

tmp <- susc_tot %>% mutate(modeltime = days)
g_tmp <- ggplot(tmp,aes(x = modeltime , y = susceptibles)) +
  theme_classic() +
  geom_line(size = 1.5) +
  labs(x = "Time",  y = "% of susceptibles",title="Susceptibles",col= scen) +
  scale_x_date(labels = scales::date_format("%m/%y"),date_breaks = datebreak) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position= "bottom") 
ggsave(file=paste0(resF,scen,"/susceptible.png"), plot=g_tmp, width=10, height=6, device ='png', scale = 0.75)

