read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

scen <-"IT1_wan5y_seas0_relax1"  
ctry <- "IT1/"
inputF = "input/"

dir.create(file.path(paste0("./",inputF), ctry), showWarnings = FALSE, recursive = TRUE)


alldata = read_excel_allsheets(paste0(inputF,scen,".xlsm"))


trans <- alldata$transmission[c(1),c(1,3)]
write.csv(trans,file = paste0(inputF,ctry,"transmission.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"transmission.csv")))

demo <- alldata$dem[,c(2,5,6)]
write.csv(demo,file = paste0(inputF,ctry,"demographics.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"demographics.csv")))

mort <- alldata$mortality[,c(3:7)]
write.csv(mort,file = paste0(inputF,ctry,"mortality.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"mortality.csv")))

birth <- alldata$birth[,c(2:6)]
write.csv(birth,file = paste0(inputF,ctry,"birth.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"birth.csv")))

contact = alldata$contact[1:N,3:(N+2)] 
write.csv(contact,file = paste0(inputF,ctry,"contact.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"contact.csv")))

distancing = alldata$distancing[,2:8]
write.csv(distancing,file = paste0(inputF,ctry,"distancing.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"distancing.csv")))

imported = alldata$imported[,2:8]
write.csv(imported,file = paste0(inputF,ctry,"imported.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"imported.csv")))

distancing_threshold = alldata$distancing_threshold[c(1:9),c(1,3)]
write.csv(distancing_threshold,file = paste0(inputF,ctry,"distancing_threshold.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"distancing_threshold.csv")))

generic = alldata$generic[c(1:6,8:9),c(1,3)]
write.csv(generic,file = paste0(inputF,ctry,"generic.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"generic.csv")))

poi = alldata$poi[c(1:18),c(1:5)]
write.csv(poi,file = paste0(inputF,ctry,"poi.csv"), row.names = FALSE)
(tmp = read.csv(paste0(inputF,ctry,"poi.csv")))

sprintf("%.20f",distancing[3,])