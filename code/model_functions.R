##########################################
#
#  Model functions  
#
############################################


library(deSolve)
library(tidyverse)
library(scales)


inputF = "input/"
resF = "result/"
codeF = "code/"
scenF = 'input/scenarios/'
figF = 'figures/'


# ************************************************************************************* #
# model outputs
# [[1]] Infections
# [[2]] Symptomatics cases
# [[3]] Model compartments
# [[4]] transition rate
# [[5]] births
# [[6]] death all causes
# [[7]] vaccinated subjects
# [[8]] Death covid
# [[9]] Hospitalisation COVID

# ************************************************************************************* #

cbPalette <- c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
cb1 <- c("#000000", "darkblue", "red", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7") 

td = as.Date("2020-01-01")
tf = as.Date("2020-12-31")

S0 = 1
wave <- 0
startyear=2020 # starting year of simulation

N <<- 9   # number of age classes
NC <<- N+1
B<-11   # number of variables per class
A<-49  # number of transitions per class
V<-NC*B # total number of variables
L<-NC*A # total number of transitions
yd <- 360 # number of days per year
md <-30  # number of days per month
dd <- 1 #  timestep in days
pd <- 5 # period for distancing and imported
dt<-dd/yd # timestep in years
per_m = md/pd
SCnovac = 'SC0'
incid = 1e6
grp<-seq(0,5*yd,pd) 

#function to read in all sheets in a workbook
read_input <- function(filelist)
{
  counter  <- 1:length(filelist$Category)
  x <- lapply(counter, function(n) read.csv(paste0(inputF,filelist$Folder[n],"/",filelist$File[n],".csv"), stringsAsFactors = F) )
  names(x) <- filelist$Category
  x
}


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
# 7=X:  Susceptible to secondary infection 
# 8=XE: Exposed to secondary infection
# 9=V:  Vaccinated - Protected
# 10=VS: Vaccinated - Susceptible 
# 11=VE: Vaccinated - exposed 

# ************************************************************************************* #
# define indices
# ************************************************************************************* #
varind <<- matrix(0,nrow=B,ncol=NC)
traind<-matrix(0,nrow=A,ncol=NC)
for (n in 1:NC){
  for (b in 1:B){
    varind[b,n]<-(n-1)*B+b
  }
  for (a in 1:A){
    traind[a,n]<-(n-1)*A+a
  }
}

# ************************************************************************************* #
# define transitions
# ************************************************************************************* #

trans <-matrix(0,nrow=L,ncol=4)
for (n in 1:NC){
  nxt=(n<N)*(n+1)+(n>=N)*(n)
  trans[traind[1,n],]<-c(varind[1,n],0,varind[1,1],1) # birth 
  trans[traind[2,n],]<-c(varind[1,n], -1,varind[1,n],0) # death S=1 
  trans[traind[3,n],]<-c(varind[2,n], -1,varind[1,n],0) # death E 
  trans[traind[4,n],]<-c(varind[3,n], -1,varind[1,n],0) # death I1 
  trans[traind[5,n],]<-c(varind[4,n], -1,varind[1,n],0) # death I2 
  trans[traind[6,n],]<-c(varind[5,n], -1,varind[1,n],0) # death I3 
  trans[traind[7,n],]<-c(varind[6,n], -1,varind[1,n],0) # death R 
  trans[traind[8,n],]<-c(varind[7,n], -1,varind[1,n],0) # death W 
  trans[traind[9,n],]<-c(varind[8,n], -1,varind[1,n],0) # death WE 
  trans[traind[10,n],]<-c(varind[9,n], -1,varind[1,n],0) # death V 
  trans[traind[11,n],]<-c(varind[10,n], -1,varind[1,n],0) # death VS 
  trans[traind[12,n],]<-c(varind[11,n], -1,varind[1,n],0) # death VE 
  
  trans[traind[13,n],]<-c(varind[1,n], -1,varind[1,nxt],+1) # age next level S 
  trans[traind[14,n],]<-c(varind[2,n], -1,varind[2,nxt],+1) # age next level E 
  trans[traind[15,n],]<-c(varind[3,n], -1,varind[3,nxt],+1) # age next level I1 
  trans[traind[16,n],]<-c(varind[4,n], -1,varind[4,nxt],+1) # age next level I2 
  trans[traind[17,n],]<-c(varind[5,n], -1,varind[5,nxt],+1) # age next level I3 
  trans[traind[18,n],]<-c(varind[6,n], -1,varind[6,nxt],+1) # age next level R 
  trans[traind[19,n],]<-c(varind[7,n], -1,varind[7,nxt],+1) # age next level W 
  trans[traind[20,n],]<-c(varind[8,n], -1,varind[8,nxt],+1) # age next level WE 
  trans[traind[21,n],]<-c(varind[9,n], -1,varind[9,nxt],+1) # age next level V 
  trans[traind[22,n],]<-c(varind[10,n], -1,varind[10,nxt],+1) # age next level VS 
  trans[traind[23,n],]<-c(varind[11,n], -1,varind[11,nxt],+1) # age next level VE 
  
  trans[traind[24,n],]<-c(varind[1,n], -1,varind[2,n],+1) # inc S -> E
  trans[traind[25,n],]<-c(varind[2,n], -1,varind[3,n],+1) # inc E -> I1
  trans[traind[26,n],]<-c(varind[2,n], -1,varind[4,n],+1) # inc E -> I2
  trans[traind[27,n],]<-c(varind[2,n], -1,varind[5,n],+1) # inc E -> I3
  trans[traind[28,n],]<-c(varind[3,n], -1,varind[6,n],+1) # recovery I1 -> R
  trans[traind[29,n],]<-c(varind[4,n], -1,varind[6,n],+1) # recovery I2 -> R
  trans[traind[30,n],]<-c(varind[5,n], -1,varind[6,n],+1) # recovery I3 -> R
  trans[traind[31,n],]<-c(varind[5,n], -1,varind[6,n],0) # specific deaths due to the disease for I3
  trans[traind[32,n],]<-c(varind[6,n], -1,varind[7,n],+1) # waning natural R -> W
  trans[traind[33,n],]<-c(varind[7,n], -1,varind[8,n],+1) # inc W -> WE
  trans[traind[34,n],]<-c(varind[8,n], -1,varind[3,n],+1) # inc WE -> I1
  trans[traind[35,n],]<-c(varind[8,n], -1,varind[4,n],+1) # inc WE -> I2
  trans[traind[36,n],]<-c(varind[8,n], -1,varind[5,n],+1) # inc WE -> I3
  
  trans[traind[37,n],]<-c(varind[1,n], -1,varind[9,n],+1) # inc S -> V
  trans[traind[38,n],]<-c(varind[1,n], -1,varind[10,n],+1) # inc S -> VS
  trans[traind[39,n],]<-c(varind[6,n], -1,varind[9,n],+1) # inc R -> V
  trans[traind[40,n],]<-c(varind[7,n], -1,varind[9,n],+1) # inc W -> V
  trans[traind[41,n],]<-c(varind[7,n], -1,varind[10,n],+1) # inc W -> VS

  trans[traind[42,n],]<-c(varind[9,n], -1,varind[7,n],+1) # inc V -> W
  trans[traind[43,n],]<-c(varind[10,n], -1,varind[11,n],+1) # inc VS -> VE
  trans[traind[44,n],]<-c(varind[11,n], -1,varind[3,n],+1) # inc VE -> I1
  trans[traind[45,n],]<-c(varind[11,n], -1,varind[4,n],+1) # inc VE -> I2
  trans[traind[46,n],]<-c(varind[11,n], -1,varind[5,n],+1) # inc VE -> I3
  trans[traind[47,n],]<-c(varind[10,n], -1,varind[7,n],+1) # inc Vs -> W
  
  trans[traind[48,n],]<-c(varind[1,n], 0,varind[1,n],0) # number vaccinated
  trans[traind[49,n],]<-c(varind[1,n], 0,varind[1,n],0) # number vaccinated
  
}

# ************************************************************************************* #
# Function to calculate transition rates, given variables and parameters
# ************************************************************************************* #
trates <- function(x, parmal, t,ti) {
  with(as.list(c( parmal)),
       {
         t_internal<-as.integer(round(t*yd)) + ti
         p_int = (floor(t_internal/pd)) %% per_m + 1 
         m_int = min(floor(t_internal/md) + 1, max_month)
         y_int = min(floor(t_internal/yd) + 1, max_yr)
         m_vac = 1+ (m_int-1) %% 12 

         #Set up matrices
         pop<-c(rep(0,N))    # population sizes  
         foif<-c(rep(0,N))     # forces of infection
         tranrate<-matrix(0,nrow=NC,ncol=A)   # transition rate matrix
         
         #population counter
         for (n in 1:N){pop[n]<-sum(x[varind[,n]])}   #  human population per age class
         poptot = sum(x)
         
         # Social distancing
         outbreak = hosp_cur/poptot
         
         age_dist <- rep(1,N)
         
         # imported cases
         ext = imported[m_int,p_int]

         # Risk-based distancing
         
        if (m_int < thrage_start |  (m_int > thrage_stop))
         # if (m_int < thresh_age) 
         {
           age_dist <- rep(1,N)
         }
         else
         {
           age_dist <- dist_age
         }

         # Threshold-based distancing
         
         if((m_int < thresh_start) |  (m_int > thresh_stop))
         {
           sdist <<- distancing[m_int,p_int]
           ext = ifelse(sdist<1,0,ext)
         }
         else if (outbreak > thresh)
         {
           ext = 0
           wave <<- 1
           if (Rt> thresh_Rt)
           {
             sdist <<- sdist * sdist_reduc
           }
         }
         else if ((outbreak > thresh * thresh_low) & (wave == 1))
         {
           ext = 0
         }
         else if ((outbreak < thresh * thresh_low) & (wave == 1))
         {
           wave <<- 0
           sdist <<- min(distancing[m_int,p_int],sdist*sdist_relax)
         }
         else if (wave == 0)
         {
           sdist <<- min(distancing[m_int,p_int],sdist*sdist_relax)
         }
         
         # seasonality
         
         t_seas = t + ti/yd
         seas = 1 + seas_amp*cospi(2*(t - seas_month*md/yd))
        # seas = 1 + seas_amp*cospi(2*(t_seas - (seas_month-0.5)/md))
         
# population                  
         
         for (n in 1:N){
           
           #Force of infection
           foif[n]<-sdist*age_dist[n]*seas*p_trans*sum(contact[,n]*(f1*x[varind[3,1:N]]+f2*x[varind[4,1:N]]+f3*x[varind[5,1:N]]))/poptot + ext 

           mortrate<-mortmat[n,y_int]
           births<-ifelse(n==1,1,0) * birthmat[1,y_int]
           pvac = ifelse(m_int > avac_month,as.numeric(vcov_y[m_vac,n]), as.numeric(vcov[m_int,n]) )
           pvac2 = ifelse(m_int > avac_month,as.numeric(vcov2_y[m_vac,n]), as.numeric(vcov2[m_int,n]) )
           
           tranrate[n,]<-c(      
             births,   # rate of birth -> R0   1
             mortrate*x[varind[1,n]],       # rate of death of S   2
             mortrate*x[varind[2,n]],       # rate of death of E   3
             mortrate*x[varind[3,n]],       # rate of death of I1  4
             mortrate*x[varind[4,n]],       # rate of death of I2  5
             mortrate*x[varind[5,n]],       # rate of death of I3  6
             mortrate*x[varind[6,n]],       # rate of death of R   7
             mortrate*x[varind[7,n]],       # rate of death of X   8
             mortrate*x[varind[8,n]],       # rate of death of XE  9
             mortrate*x[varind[9,n]],       # rate of death of V   10
             mortrate*x[varind[10,n]],       # rate of death of VS  11
             mortrate*x[varind[11,n]],       # rate of death of VE  12
             
             (1-mortrate)*agerate[n]*x[varind[1,n]],        # aging s              13
             (1-mortrate)*agerate[n]*x[varind[2,n]],        # aging E              14
             (1-mortrate)*agerate[n]*x[varind[3,n]],        # aging I1              15
             (1-mortrate)*agerate[n]*x[varind[4,n]],        # aging I2              16
             (1-mortrate)*agerate[n]*x[varind[5,n]],        # aging I3              17
             (1-mortrate)*agerate[n]*x[varind[6,n]],        # aging R              18
             (1-mortrate)*agerate[n]*x[varind[7,n]],        # aging X              19
             (1-mortrate)*agerate[n]*x[varind[8,n]],        # aging XE              20
             (1-mortrate)*agerate[n]*x[varind[9,n]],        # aging V              21
             (1-mortrate)*agerate[n]*x[varind[10,n]],        # aging VS              22
             (1-mortrate)*agerate[n]*x[varind[11,n]],        # aging VE              23
             
             foif[n]*x[varind[1,n]],     # inc S -> E     24
             poiE[n,1]*epsilon*x[varind[2,n]],     # inc E -> I1   25
             poiE[n,2]*epsilon*x[varind[2,n]],     # inc E -> I2   26
             poiE[n,3]*epsilon*x[varind[2,n]],     # inc E -> I3   27
             delta*x[varind[3,n]],  # recovery I1 -> R     28
             delta*x[varind[4,n]],  # recovery I2 -> R     29
             (1-poiE[n,4])*delta*x[varind[5,n]],   # recovery I3 -> R     30
             poiE[n,4]*delta*x[varind[5,n]], # death I3      31
             tau*x[varind[6,n]],               # natural waning R -> X  32
             foif[n]*x[varind[7,n]],     # inc X -> XE     33
             poiWE[n,1]*epsilon*x[varind[8,n]],     # inc XE -> I1   34
             poiWE[n,2]*epsilon*x[varind[8,n]],     # inc XE -> I2   35
             poiWE[n,3]*epsilon*x[varind[8,n]],     # inc XE -> I3   36
             
             es[vn[m_int],1]*pvac*x[varind[1,n]], # vacc S->V    37                                  
             es[vn[m_int],2]*pvac*x[varind[1,n]], # vacc S->VS    38                                  
             er[vn[m_int],1]*pvac*x[varind[6,n]], # vacc R->V    39
             ew[vn[m_int],1]*pvac*x[varind[7,n]], # vacc  W->V  40                                  
             ew[vn[m_int],2]*pvac*x[varind[7,n]], # vacc  W->VS  41
             
             wan[vn[m_int]]*x[varind[9,n]],               # vaccine waning V -> W  42
             foi_v[vn[m_int]] * foif[n]*x[varind[10,n]],     # inc VS -> VE     43
             poiVE[vn[m_int],n,1]*epsilon*x[varind[11,n]],     # inc VE -> I1   44
             poiVE[vn[m_int],n,2]*epsilon*x[varind[11,n]],     # inc VE -> I2   45
             poiVE[vn[m_int],n,3]*epsilon*x[varind[11,n]],     # inc VE -> I3   46
             wan[vn[m_int]]*x[varind[10,n]],               # vaccine waning VS -> W  47
             pvac2*pop[n],  # number vaccinated 48
             sdist*age_dist[n]  # distancing 49
           )
         }
# cohort          
          if(m_int < mdeb)
          {
            tranrate[NC,]<-rep(0,A)      
          }
          else
          {
            mortrate<-mortmat[ncoh,y_int]
            births<-0
            pvac = as.numeric(vcov_c[m_int])
            
            tranrate[NC,]<-c(      
              births,   # rate of birth -> R0   1
              mortrate*x[varind[1,NC]],       # rate of death of S   2
              mortrate*x[varind[2,NC]],       # rate of death of E   3
              mortrate*x[varind[3,NC]],       # rate of death of I1  4
              mortrate*x[varind[4,NC]],       # rate of death of I2  5
              mortrate*x[varind[5,NC]],       # rate of death of I3  6
              mortrate*x[varind[6,NC]],       # rate of death of R   7
              mortrate*x[varind[7,NC]],       # rate of death of X   8
              mortrate*x[varind[8,NC]],       # rate of death of XE  9
              mortrate*x[varind[9,NC]],       # rate of death of V   10
              mortrate*x[varind[10,NC]],       # rate of death of VS  11
              mortrate*x[varind[11,NC]],       # rate of death of VE  12
              
              (1-mortrate)*agerate[ncoh]*x[varind[1,NC]],        # aging s              13
              (1-mortrate)*agerate[ncoh]*x[varind[2,NC]],        # aging E              14
              (1-mortrate)*agerate[ncoh]*x[varind[3,NC]],        # aging I1              15
              (1-mortrate)*agerate[ncoh]*x[varind[4,NC]],        # aging I2              16
              (1-mortrate)*agerate[ncoh]*x[varind[5,NC]],        # aging I3              17
              (1-mortrate)*agerate[ncoh]*x[varind[6,NC]],        # aging R              18
              (1-mortrate)*agerate[ncoh]*x[varind[7,NC]],        # aging X              19
              (1-mortrate)*agerate[ncoh]*x[varind[8,NC]],        # aging XE              20
              (1-mortrate)*agerate[ncoh]*x[varind[9,NC]],        # aging V              21
              (1-mortrate)*agerate[ncoh]*x[varind[10,NC]],        # aging VS              22
              (1-mortrate)*agerate[ncoh]*x[varind[11,NC]],        # aging VE              23
              
              foif[ncoh]*x[varind[1,NC]],     # inc S -> E     24
              poiE[ncoh,1]*epsilon*x[varind[2,NC]],     # inc E -> I1   25
              poiE[ncoh,2]*epsilon*x[varind[2,NC]],     # inc E -> I2   26
              poiE[ncoh,3]*epsilon*x[varind[2,NC]],     # inc E -> I3   27
              delta*x[varind[3,NC]],  # recovery I1 -> R     28
              delta*x[varind[4,NC]],  # recovery I2 -> R     29
              (1-poiE[ncoh,4])*delta*x[varind[5,NC]],   # recovery I3 -> R     30
              poiE[ncoh,4]*delta*x[varind[5,NC]], # death I3      31
              tau*x[varind[6,NC]],               # natural waning R -> X  32
              foif[ncoh]*x[varind[7,NC]],     # inc X -> XE     33
              poiWE[ncoh,1]*epsilon*x[varind[8,NC]],     # inc XE -> I1   34
              poiWE[ncoh,2]*epsilon*x[varind[8,NC]],     # inc XE -> I2   35
              poiWE[ncoh,3]*epsilon*x[varind[8,NC]],     # inc XE -> I3   36
              
              es[vn[m_int],1]*pvac*x[varind[1,NC]], # vacc S->V    37                                  
              es[vn[m_int],2]*pvac*x[varind[1,NC]], # vacc S->VS    38                                  
              er[vn[m_int],1]*pvac*x[varind[6,NC]], # vacc R->V    39
              ew[vn[m_int],1]*pvac*x[varind[7,NC]], # vacc  w->V  40                                  
              ew[vn[m_int],2]*pvac*x[varind[7,NC]], # vacc  w->VS  41
              
              wan[vn[m_int]]*x[varind[9,NC]],               # vaccine waning V -> VS  42
              foi_v[vn[m_int]] * foif[ncoh]*x[varind[10,NC]],     # inc VS -> VE     43
              poiVE[vn[m_int],ncoh,1]*epsilon*x[varind[11,NC]],     # inc VE -> I1   44
              poiVE[vn[m_int],ncoh,2]*epsilon*x[varind[11,NC]],     # inc VE -> I2   45
              poiVE[vn[m_int],ncoh,3]*epsilon*x[varind[11,NC]],     # inc VE -> I3   46
              wan[vn[m_int]] *x[varind[10,NC]],               # vaccine waning Vs -> W  47
              pvac,  # number vaccinated 48
              sdist*age_dist[ncoh]  # distancing 49
            )
            
          }
          rate <- c(t(tranrate))
          # Rt calculation
          hosp_prec <<- hosp_cur
          hosp_cur <<-sum (rate[c(traind[c(27,36,46),])])*dt
          Rho <- ifelse(hosp_cur>0,log(hosp_cur),0) - ifelse(hosp_prec>0,log(hosp_prec),0)
          Rt <<- (Rho/dt+epsilon)*(Rho/dt+delta)/delta/epsilon

         return(rate)
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
         vacc<-sev<-deathC<-deathN<-birth<-symp<-inf<-matrix(0,nrow=length(out[,1]),ncol=NC)
         for (n in 1:NC){
           inf[,n]<-rowSums(tran[,c(traind[c(c(25:27),c(34:36),c(44:46)),n])])*dt
           symp[,n]<-rowSums(tran[,c(traind[c(c(26:27),c(35:36),c(45:46)),n])])*dt
           birth[,n]<-tran[,c(traind[49,n])]
           deathN[,n]<-rowSums(tran[,c(traind[c(2:12),n])])*dt
           vacc[,n]<-tran[,c(traind[48,n])]*dt
           deathC[,n]<-tran[,c(traind[31,n])]*dt
           sev[,n]<-rowSums(tran[,c(traind[c(27,36,46),n])])*dt
         }
         
         return(cbind(inf,    #1
                      symp,#2
                      birth,      #3
                      deathN,       #4
                      vacc,           #5
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
         
         # Rt calculation
         # hosp_prec <<- hosp_cur
         # hosp_cur <<-sum (transit[c(traind[c(27,36,46),])])*dt
         # Rho <- ifelse(hosp_cur>0,log(hosp_cur),0) - ifelse(hosp_prec>0,log(hosp_prec),0)
         # Rt <<- (Rho/dt+epsilon)*(Rho/dt+delta)/delta/epsilon
         # cat(paste(hosp_cur, Rt),'\n')
         
         eq<-rep(0.0, V)
         
         eq<-EQ2(L, NC, eq, transit,trans)
         
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
      WE=0*parrun$population[n], 
      VP=0*parrun$population[n], 
      VS=0*parrun$population[n], 
      VE=0*parrun$population[n] 
    ))
   }
    initcondrun<-c(initcondrun,c(
      S=parrun$sc[1], 
      E=parrun$sc[2], 
      I1=parrun$sc[3], 
      I2=parrun$sc[4], 
      I3=parrun$sc[5], 
      R=parrun$sc[6], 
      W=parrun$sc[7], 
      WE=parrun$sc[8], 
      VP=parrun$sc[9], 
      VS=parrun$sc[10], 
      VE=parrun$sc[11] 
    ))
  #initoderun<-round(initcondrun)
  initoderun<-initcondrun
  staterun <- c(initcondrun,0)
  ti<-0  
  tyears <- parrun$tmonth/12 
  tsteps<-round(tyears/dt) # number of time steps
  
  transitrun <- trates(initoderun,parrun,0,ti)
  
  
  # SOLVE THE ODEs and get output
  
  timesrun <- seq(0, tyears, by = dt) # Model run time
  #Solve ODE
  outoderun <- ode(y = staterun, times = timesrun, func = epiModel, parms = parrun, method  = "euler", hini=0.01)
  # Compute transitions at each time step
  tranoderun<-matrix(0,nrow=length(outoderun[,1]),ncol=L)
  

  for (ti in 0:tsteps){
    tranoderun[ti+1,]<-t(trates(outoderun[ti+1,2:(1+V)],parrun,0,ti))
  }
  #Compute outputs
  ppoutrun<-postproc(parrun,outoderun,tranoderun)
  modeltimes<-outoderun[,1]+startyear

  inf_ode<-ppoutrun[,1:NC]
  symp_ode<-ppoutrun[,(NC+1):(2*NC)]
  dist_ode<-ppoutrun[,(2*NC+1):(3*NC)]
  deathN_ode<-ppoutrun[,(3*NC+1):(4*NC)]
  out_ode<-ppoutrun[,(4*NC+1):(5*NC)]
  deathC_ode<-ppoutrun[,(5*NC+1):(6*NC)]
  hosp_ode<-ppoutrun[,(6*NC+1):(7*NC)]
  
  return(list( inf_ode, #1
               symp_ode,    #2
               outoderun,    #3
               tranoderun,    #4
               dist_ode, #5
               deathN_ode,  #6
               out_ode,  #7               
               deathC_ode,  #8
               hosp_ode #9
  ) )
}



inputdata <-function(ctry,vacscen,tmonth) 
{
  
  # ************************************************************************************* #
  # import Scenario data
  # ************************************************************************************* #
  
  tmp1 <-read.csv(paste0(scenF,ctry,".csv"))
  tmp2 <-read.csv(paste0(scenF,vacscen,".csv"))
  scenario = rbind(tmp1,tmp2)
  resFolder = tmp1[1,3]

  alldata = read_input(scenario)
  
  # demographic parameters 
  
  catname <- alldata$demo[,2]
  catname <<- c(catname,"coh")
  
  agestart<-alldata$demo[,3] #0k
  ageend<-alldata$demo[,4] #0k
  agerate<-alldata$demo[,5] #0k
  population<-alldata$demo[,6] #0k
  mortmat<-alldata$mortality[,3:ncol(alldata$mortality)] #0k
  birthmat<-alldata$birth[,2:ncol(alldata$birth)] #0k
  
  # agerate<- mortmat[,1]/(exp(mortmat[,1]/agerate)-1)

  # dem_tf<-as.numeric(names(birthmat)) #ok
  
  # disease parameters 
  
  parameters<-alldata$disease #ok
  epsilon=yd/parameters$value[4] 
  delta=yd/parameters$value[5] 
  tau=log(2)/parameters$value[6]
  poiWE <- array(0,dim=c(N,4))
  poiWE = as.matrix(alldata$poi[1:N,2:5])
  poiWE[,2:3] = poiWE[,2:3] * (1-parameters$value[7])
  poiWE[,1] = 1 -poiWE[,2] - poiWE[,3]
  

  # npi parameters 
  
  npi<-alldata$npi #ok
  
  hosp_cur <<- 1e-10
  hosp_prec <<- 1e-10
  sdist <<- 1
  wave <<- 0
  
  # vaccine parameters
  
  vaccine <- alldata$vaccine
  
  poiVE <- array(0,dim=c(2,N,4))
  es  <- array(0,dim=c(2,2))
  ew  <- array(0,dim=c(2,2))
  er  <- array(0,dim=c(2,2))
  foi_v <- rep(0,2)
  wan <- rep(0,2)
  
  # vaccine #1
  
  poiVE[1,,] = as.matrix(alldata$poi[1:N,2:5])
  poiVE[1,,2:3] = poiVE[1,,2:3] * (1-vaccine$vaccine1[7])
  poiVE[1,,1] = 1 -poiVE[1,,2] - poiVE[1,,3]
  es[1,1] = vaccine$vaccine1[1] 
  es[1,2] = vaccine$vaccine1[2] 
  ew[1,1] = vaccine$vaccine1[1] 
  ew[1,2] = vaccine$vaccine1[3] 
  er[1,1] = vaccine$vaccine1[4] 
  foi_v[1] = 1 - vaccine$vaccine1[5]  
  wan[1] = log(2)/vaccine$vaccine1[6]  
  
  # vaccine #2
  
  poiVE[2,,] = as.matrix(alldata$poi[1:N,2:5])
  poiVE[2,,2:3] = poiVE[2,,2:3] * (1-vaccine$vaccine2[7])
  poiVE[2,,1] = 1 -poiVE[2,,2] - poiVE[2,,3]
  es[2,1] = vaccine$vaccine2[1] 
  es[2,2] = vaccine$vaccine2[2] 
  ew[2,1] = vaccine$vaccine2[1] 
  ew[2,2] = vaccine$vaccine2[3] 
  er[2,1] = vaccine$vaccine2[4] 
  foi_v[2] = 1 - vaccine$vaccine2[5]  
  wan[2] = log(2)/vaccine$vaccine2[6]  
  
  vn <-alldata$program[,12] #ok
  coverage<-alldata$program[,3:11] #ok
  max_month = dim(coverage)[1] 
  
  vcov<-coverage
  vcov[1,] <- 0 
  for (m in 2:dim(coverage)[1]){vcov[m,] <- 12 * log((1-coverage[m-1,])/(1-coverage[m,]))  }   
  vcov2<-coverage*12
  for (m in 2:dim(coverage)[1]){vcov2[m,] <-(coverage[m,]-coverage[m-1,])*12  }   
  
  vcov <- vcov %>% mutate_all(~lag(.x,vaccine$vaccine1[8]))
  vcov[is.na(vcov)] <- 0
  
  coverage_y <- alldata$routine[,2:10] #ok
  vcov_y<-coverage_y
  vcov_y[1,] <- log(1/(1-coverage_y[1,]))
  for (m in 2:12){vcov_y[m,] <- 12 * log((1-coverage_y[m-1,])/(1-coverage_y[m,]))  }
  vcov2_y <-coverage_y*12
  for (m in 2:dim(coverage_y)[1]){vcov2_y[m,] <-(coverage_y[m,]-coverage_y[m-1,])*12  }   
  
  
  # cohort data
  
  popcoh <- alldata$cohort$value[2]
  population[NC] <- popcoh
  vcoh <- alldata$cohort$value[3]
  mcoh <- alldata$cohort$value[4]
  mdeb = alldata$cohort$value[5]
  
  vcov_c <- rep(0,dim(coverage)[1])
  
  sc <- rep(0,B)
  sc[alldata$cohort$value[6]] <- popcoh
  
  if (vcoh==1)
  {
    if (mdeb < mcoh)
    {
      vcov_c[mcoh] <- 200
    }
    else
    {
      sc[9] = es*sc[1] + er*sc[6] + er*sc[7]
      sc[10] = (1-es)*sc[1] + (1-ew)*sc[7]
      sc[1] = 0
      sc[7] = 0
      sc[6] = (1-er)*sc[6]
      
    }
  }
  
  # calibrated parameters 
  
  calval <- alldata$calibrated$value
  calper <- alldata$calibrated$period
  calval <- calval[calper>0]
  calper <- calper[calper>0]
  scale_imp <- as.numeric(alldata$calibrated[6,4])
  lag_death <- as.numeric(alldata$calibrated[9,4])
  lag_hosp <- as.numeric(alldata$calibrated[10,4])

  calibrated <-cbind(period = calper, value = calval)
  imported <- array(calval[2],dim=c(max_month,per_m))

  distancing <- array(1,dim=c(max_month,per_m))
  relax <- 4
  calf <- calper[length(calper)]
  thresh_start = ceiling(calf/per_m) + 1
  
  tmp <- calf:(calf+relax*per_m)
  cper = c(calper[-1],tmp[-1],max_month*per_m)
  tmp <- seq (calval[length(calval)],1,length.out = relax*per_m)
  cval = c(1,calval[-(1:2)],tmp,1)
  num = 1
  for (m in 1:max_month) {
    for (p in 1:per_m) {
      val = (m-1)*per_m + p
      if(val > cper[num]) {num = num + 1}
      distancing[m,p] <- cval[num]
    }
  }
  
  pars=list(
    resFolder = resFolder,
    p_trans=calval[1]*yd,
    f1=parameters$value[1], 
    f2=parameters$value[2], 
    f3=parameters$value[3], 
    epsilon=epsilon, 
    delta=delta, 
    tau=tau, 
    seas_month = parameters$value[8], 
    seas_amp = parameters$value[9], 
    es=es, 
    ew=ew, 
    er=er, 
    foi_v = foi_v,  
    wan= wan,  
    avac_month = alldata$routine_start$value[1],
    max_month = dim(vcov)[1], 
    max_yr = dim(mortmat)[2], 
    thresh_start = thresh_start,  
    thrage_start = npi$value[1],  
    thrage_stop = npi$value[2],  
    thresh_stop = npi$value[8], 
    thresh = npi$value[3], 
    thresh_low = npi$value[4], 
    sdist_reduc = npi$value[5], 
    sdist_relax = npi$value[6], 
    thresh_Rt = npi$value[7], 
    poiE = alldata$poi[1:N,2:5], 
    poiWE = poiWE, 
    poiVE = poiVE, 
    contact = alldata$contact[1:N,3:(N+2)], 
    distancing = distancing, 
    imported  = imported, 
    vcov = vcov,
    vcov_y = vcov_y ,
    vcov2 = vcov2,
    vcov2_y = vcov2_y ,
    population = population ,
    agerate = agerate,
    mortmat = mortmat ,
    birthmat  = birthmat ,
    dist_age = alldata$distancing_age$value[1:N],
    ncoh = alldata$cohort$value[1],
    mdeb = mdeb,
    sc = sc,
    vcov_c = vcov_c,
    tmonth = tmonth,
    vn = vn,
    calibrated = calibrated,
    lag_death = lag_death,
    lag_hosp = lag_hosp
    
  ) #ok
  
  return(pars)
}

run_scenario <- function (ctry,vacscen,tmonth)
{

  scen = paste0(ctry,vacscen)

  pars = inputdata(ctry,vacscen,tmonth)
  rescenF = paste0(pars$resFolder,"/")
  run1<-runmod(pars)

  # infections
  inf <-data.frame(run1[[1]])
  colnames(inf) = c(catname)
  write.csv(inf,file = paste0(resF,rescenF,scen,"_infections.csv"), row.names = FALSE)

  # cases
  cases <-data.frame(run1[[2]])
  colnames(cases) = c(catname)
  write.csv(cases,file = paste0(resF,rescenF,scen,"_cases.csv"), row.names = FALSE)

  # hospit
  hosp <-data.frame(run1[[9]])
  colnames(hosp) = c(catname)
  hosp <- hosp %>% mutate_all(~ lag(.x,pars$lag_hosp))
  hosp[1:pars$lag_hosp,] = 0 
  write.csv(hosp,file = paste0(resF,rescenF,scen,"_hospitalisations.csv"), row.names = FALSE)

  # covid deaths
  deathC <-data.frame(run1[[8]])
  colnames(deathC) = c(catname)
  deathC <- deathC %>% mutate_all(~ lag(.x,pars$lag_death))
  deathC[1:pars$lag_death,] = 0 
  write.csv(deathC,file = paste0(resF,rescenF,scen,"_death COVID.csv"), row.names = FALSE)

  # all compartments
  state <-data.frame(run1[[3]])
  write.csv(state,paste0(resF,rescenF,scen,file = "_state.csv"), row.names = FALSE)
  
  # number vaccinated
  nvac <-data.frame(run1[[7]])
  colnames(nvac) = c(catname)
  write.csv(nvac,paste0(resF,rescenF,scen,file = "_vaccinated.csv"), row.names = FALSE)
  
  # number vaccinated
  distage <-data.frame(run1[[5]])
  colnames(distage) = c(catname)
  write.csv(distage,paste0(resF,rescenF,scen,file = "_distancing.csv"), row.names = FALSE)

  return(scen)

}

write_results <- function (run1)
{
  
  # infections
  inf <-data.frame(run1[[1]])
  colnames(inf) = c(catname)
  write.csv(inf,file = paste0(resF,rescenF,scen,"_infections.csv"), row.names = FALSE)
  
  # cases
  cases <-data.frame(run1[[2]])
  colnames(cases) = c(catname)
  write.csv(cases,file = paste0(resF,rescenF,scen,"_cases.csv"), row.names = FALSE)
  
  # hospit
  hosp <-data.frame(run1[[9]])
  colnames(hosp) = c(catname)
  hosp <- hosp %>% mutate_all(~ lag(.x,pars$lag_hosp))
  hosp[1:pars$lag_hosp,] = 0 
  write.csv(hosp,file = paste0(resF,rescenF,scen,"_hospitalisations.csv"), row.names = FALSE)
  
  # covid deaths
  deathC <-data.frame(run1[[8]])
  colnames(deathC) = c(catname)
  deathC <- deathC %>% mutate_all(~ lag(.x,pars$lag_death))
  deathC[1:pars$lag_death,] = 0 
  write.csv(deathC,file = paste0(resF,rescenF,scen,"_death COVID.csv"), row.names = FALSE)
  
  # all compartments
  state <-data.frame(run1[[3]])
  write.csv(state,paste0(resF,rescenF,scen,file = "_state.csv"), row.names = FALSE)
  
  # number vaccinated
  nvac <-data.frame(run1[[7]])
  colnames(nvac) = c(catname)
  write.csv(nvac,paste0(resF,rescenF,scen,file = "_vaccinated.csv"), row.names = FALSE)
  
  # number vaccinated
  distage <-data.frame(run1[[5]])
  colnames(distage) = c(catname)
  write.csv(distage,paste0(resF,rescenF,scen,file = "_distancing.csv"), row.names = FALSE)
  
  return(scen)
  
}

# run_scenario(ctry,scenv,tmonth)
# inputdata(ctry,vacscen,tmonth)
# run1<-runmod(pars)

#test <- run_optim(scen,svac,scoh,scenL,tmonth)
