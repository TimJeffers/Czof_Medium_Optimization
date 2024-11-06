# Media Optimization Script


# Set up ----
#Set the working directory to file path you downloaded this code 
setwd("~/Documents/Czof_Media_Optimization/")
dir.create("OptimizedMediaFigs")

## load ---- 
library(tidyverse)
library(plotly); library(pastecs); library(tweenr); library(directlabels)
library(ggsignif)
library(lubridate)

## Functions ---- 

## standard error
standard.error<- function(x){
  sd(x)/sqrt(length(x))
}
## CoulterCount counts all the data
CellCounter2 = function(M3Folder){
  # read in file directories 
  filelist <- list.files(M3Folder) ; 
  filelist<- filelist[!grepl("av" , filelist)] #This is used to ignore any "av"  files
  
  #set up dataframe objects () 
  counterdata <- NULL
  culturedata<-data.frame(ID=0,Date=0, Rep=0, Day=0, CellDensity=0, Dilution=0,  Mean=0, P1=0, P2=0,
                          P3=0, TotVol=0, filename=0, Experiment = 0) 
  
  # For loop over files: data extracted from M3 files, including raw diameter
  #and number measurements and calculates 
  #the cell density, totalvol, and give approximations of the peaks according 
  #to cellular distribution 
  
  for(i in 1:length(filelist)){
    culturedata[i, "filename"]= filelist[i]
    #file is saved into the document even if there is an annotation error so we
    #can  see where the error is for now 
    
    
    #Extract useful Date and Dilution information from file name 
    s <- strsplit(filelist[i], "_" )
    ID <- s[[1]][2]
    Dilution <- as.numeric(s[[1]][3])
    Date <- s[[1]][1]
    Rep <- gsub(".#m3", "",  s[[1]][4])
    
    #  skip the file if there is any missing information;
    #the any function is looking for any NAs in the list 
    #we could add a next function here if ANY of the 
    if(any(is.na(c(Date, Dilution))))  next 
    
    #Extract  BinHeight and BinDiam from each file added 
    Practice <- read.csv(file=paste(M3Folder,"/",filelist[i], sep=""),
                         row.names=NULL)
    Practice <- Practice %>% mutate(RowNum = as.numeric(row.names(Practice))) 
    
    #The following steps extract the data related to diameter of cells BinDIAM and the number in those (BinHeight)
    Points <- filter(Practice, row.names=="[#Bindiam]"| row.names=="[#Binheight]"|row.names=="[Binunits]" |row.names=="[SizeStats]")$RowNum  #This lists the values where the numeric values start and stop.
    
    BinDiam = as.numeric(Practice$row.names[(Points[1]+1):(Points[2]-1)])
    BinHeight = as.numeric(Practice$row.names[(Points[3]+1):(Points[4]-1)]) *Dilution 
    
    #Let's make a file with these two things together 
    cell.distrib <- data_frame(ID,  Date, BinDiam, BinHeight) %>% filter(BinDiam>2) %>% mutate(Experiment= M3Folder, Rep= Rep) 
    counterdata <- rbind(counterdata, cell.distrib)
    #counterdata will build the whole thing
    cell.distrib.mean <- cell.distrib %>% 
      filter(BinDiam > 2 & BinDiam  < 20)  %>% 
      #We almost always look in between these two datapoints for this 
      summarise(Count = sum(BinHeight),
                Mean=sum(BinHeight*BinDiam)/sum(BinHeight),
                Vol = sum(4*pi*(BinDiam/2)^3*BinHeight/3))
    
    #Input values into culture data dataframe 
    culturedata[i,"Dilution"] <- Dilution
    culturedata[i,"ID"] <-ID
    culturedata[i, "filename"]= filelist[i]
    culturedata[i,"Date"] <-Date
    culturedata[i,"CellDensity"]<- cell.distrib.mean$Count[1]
    culturedata[i,"Mean"]<- cell.distrib.mean$Mean[1]
    culturedata[i,"TotVol"]<- cell.distrib.mean$Vol[1]
    culturedata[i,"Rep"] <- Rep
    culturedata[i,"Experiment"] <- M3Folder
    
    ## Modeling distributions by polynomial distribution to find local peaks
    ## (Not essential for this publications)
    poly <- lm( BinHeight ~ poly(BinDiam, 50, raw=TRUE),cell.distrib %>% 
                  filter(BinDiam > 2 & BinDiam  < 20) )
    
    #Fitting that data with the fitted function:  
    CellandFit <- cell.distrib %>% filter(BinDiam > 2 & BinDiam  < 20)  %>% 
      mutate(fitted=fitted(poly))
    tp2 <- turnpoints(CellandFit$fitted)
    
    #tp2$peaks identifies all peaks in the regression. Tp$tppos identifies all turning points.Returns a true false schematic. 
    CellandFit$BinDiam[tp2$peaks] #
    pd <- data.frame(Diam=CellandFit$BinDiam[tp2$peaks], 
                     Height  =CellandFit$BinHeight[tp2$peaks] )
    pd <- pd %>% arrange(desc(Height))
    
    
    #Calls local peacs from the PD diam and puts these in the culture data
    culturedata[i,"P1"]<- pd$Diam[1]
    culturedata[i,"P2"]<- if(pd$Height[2] != 0) {pd$Diam[2]} else{0}
    culturedata[i,"P3"]<- if(pd$Height[3] != 0) {pd$Diam[3]} else{0}
    
    
  }
  return(list(culturedata, counterdata))
}


## This allows all variables to be represented with zero 
ZeroProcess_Data <- function(data, group_col) {
  # Filter data for Day == 0 and remove ID and group_col
  zero <- data %>%
    filter(Day == 0) %>%
    ungroup() %>%
    select(-ID, -{{group_col}})
  
  # Create zero_examp data frame
  zero_examp <- data %>%
    ungroup() %>%
    distinct({{group_col}}) %>%
    # Set the variable as 0 days, the same format as subrtracting two dates together. Pick the same date so it is 0
    mutate(Day = as.Date(as.character("20230101"), "%Y%m%d")-as.Date(as.character("20230101"), "%Y%m%d")) %>%
    ## I believe this part is the same 
    inner_join(zero, by = c( "Day")) %>%
    bind_rows(data %>% filter(Day != 0)) 
  
  return(zero_examp)
}


## LoopCount is updated here to add trycatch in case there is a problem with that 
LoopCount = function(path, x, inoc, keyfile){
  
  
  # try-catch block to handle errors
  result <- tryCatch({
    
    setwd(path)
    
    # key file join
    key <- read.csv(keyfile, header=TRUE)
    key$ID <- as.character(key$ID)
    
    
    CountLoop = CellCounter2(x)
    
    # join the inoc file and convert to dates in the third
    CountLoop[[3]] <- CountLoop[[1]] %>% 
      mutate(Date = as.Date(as.character(Date), "%Y%m%d")) %>% 
      group_by(ID, Date) %>% 
      dplyr::summarise(CellDensity = mean(CellDensity), 
                       Mean = mean(Mean), 
                       P1 = mean(P1), 
                       P2 = mean(P2), 
                       P3 = mean(P3), 
                       TotVol = mean(TotVol)) %>% 
      mutate(Day = Date - as.Date(as.character(inoc), "%Y%m%d")) %>%
      inner_join(key)
    return(CountLoop)
  }, error = function(e) {
    return(paste("Error:", e$message))
  })
  
  return(result)
}

#  Long Term Glucose -----
cult.dat <- read.csv(file= "20190318_ADJMgS_GlcCons_GrowthData.csv", header=T)
keyfile="20190225_ADJMgS_key.csv"
key <- read.csv(keyfile, header=TRUE)
key$ID <- as.character(key$ID)
inoc <- "20190225"
## Plot the cult dat ----
cult.dat %>% mutate(Glc = ifelse(Glc == "0glc", "Autotrophic", "Glucose")) %>% 
  filter(!(Glc=="Glucose" & Day %in% c(9,10) & TotVol < 1.2E10)) %>%  #Here are the measure sample outliers 
  #filter(Media=="adj2x") 
  ggplot(aes(x=Day, y=TotVol, fill=Glc, linetype=Glc, shape=Glc)) + 
  stat_summary(geom="errorbar", color="black", linetype="solid") +
  stat_summary(geom="line", color="dark green") +
  geom_jitter( alpha=0.8, size=1.5, width=0.2) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24) )+
  theme_classic(base_size=10)+ # geom_point( color="dark green") +
  scale_fill_manual(values = c("Autotrophic"="#4ef542", "Glucose"="dark green") )+
  scale_y_log10(breaks=c(1E5, 1E6,1E7,1E8,1E9,1E10,1E11,1E12,1E13), 
                labels=c(1,10,100,1000,10000,100000,1000000,10000000,100000000 )/10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9," ", mu*m^{3}, " ", mL^{-1} ))) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=45))

ggsave("OptimizedMediaFigs/GlcucoseAdditionADD.pdf", device= "pdf", unit="in", width=4,height=3)

###Glc Doubling time:: 
start  <- (cult.dat %>% filter(Day==1) %>% group_by(Day) %>% dplyr::summarise(mean=mean(TotVol)))$mean
photo.td <- cult.dat %>% filter(Glc=="0glc")  %>% filter(Day==1| Day==6) %>%   group_by(ID,  Glc) %>%
  dplyr::mutate(mu=(log(TotVol)-log(start))/5) %>% mutate(td = (log(2)/mu)*24) %>%filter(Day==6) %>%
  select( Glc, td)# %>% group_by(Media, Glc) %>% dplyr::summarise(mean=mean(td), se=standard.error(td))

#Calculating doubling time. Use the inner_join function so that you can substract each by the average of the TimePoint 0 point
td.key <- cult.dat %>% filter(Glc=="Glc")  %>% filter(Day==5) %>% group_by(Glc) %>% dplyr::summarize(start=mean(TotVol)) 
mixo.td <- cult.dat %>% filter(Glc=="Glc")  %>% filter(Day==8) %>% inner_join(td.key) %>%   group_by(ID, Glc) %>%
  dplyr::mutate(  mu=(log(TotVol)-log(start))/3) %>% mutate(td = (log(2)/mu)*24) %>%
  select( Glc, td) #%>%group_by(Media, Glc) %>% dplyr::summarise(mean=mean(td), se=standard.error(td))


bind_rows(photo.td, mixo.td) %>% 
  mutate(Glc = ifelse(grepl("0",Glc), "Autotrophic", "Glucose")) %>% 
  #filter(Glc == "Glucose") %>% 
  ggplot() + 
  stat_summary(aes(x=Glc, y= td, fill=Glc), color="black", geom = "bar") +
  stat_summary(aes(x=Glc, y= td, fill=Glc), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Glc, y= td, fill=Glc, shape=Glc), color="black", width=0.1) +
  #facet_wrap(~Glc) +
  scale_fill_manual(values = c("Autotrophic"="#4ef542", "Glucose"="dark green") )+
  theme_classic(base_size=10) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  scale_fill_manual(values = c("Autotrophic"="#4ef542", "Glucose"="dark green") )+
  scale_y_continuous(breaks=c(0,6,12, 18,24)) +
  ylab("Doubling time (h)") + xlab(NULL) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("OptimizedMediaFigs/GlcucoseAddition_Doubling.pdf", device= "pdf", unit="in", width=2,height=3)


## glucose measurements -----
Glc <- read.csv("20190303_MgS_GlcKey.csv", header=T)
Glc$ID <- as.character(Glc$ID)

Glc <- Glc %>% mutate(Date=as.Date(as.character(Date), "%Y%m%d")) %>%
  inner_join(key)  %>% mutate(Day= Date- as.Date(as.character(inoc), "%Y%m%d"))


## Figure
ggplot(data=Glc %>% filter(Time == ""), aes(x=Day, y=mMGLC)) + 
  geom_point(aes(x=Day, y= CumulativeAddedGlc), size=3, color="purple")+
   stat_summary(aes(x=Day, y= mMGLC), color="black", geom = "errorbar", width=0.5) +
  
  geom_jitter( color="black",fill="dark green", width=0.1, shape=24) +
  theme_classic(base_size=10) +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("mM Glc"))) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=45))
ggsave("OptimizedMediaFigs/Glucose_Consumption.pdf", device= "pdf", unit="in", width=4,height=3)

Glc.mean <- Glc  %>% mutate(Date=as.Date(as.character(Date), "%Y%m%d")) %>% group_by(ID, Date,CumulativeAddedGlc,Time) %>%
  dplyr::summarize(mMGlc = mean(mMGLC)) %>% inner_join(key)  %>% mutate(Day= Date- as.Date(as.character(inoc), "%Y%m%d"))


# ICP Gradient Boost----

## Per cell data----
icp.by.cell <-read.csv("20190319_ICP_bycell_MgS.csv", header=T)
icp.by.cell$ID <- as.character(icp.by.cell$ID)
icp.by.cell$ID[1:3] <- c("02", "03", "04")

nutr_levels = icp.by.cell %>% 
  gather(element, concentration, -ID, -Date ) %>%
  group_by(element) %>% dplyr::summarise(concentration = mean(concentration)) %>% 
  arrange(desc(concentration))

icp.by.cell %>% 
gather(element, concentration, -ID, -Date ) %>%
  mutate(Day = as.Date(as.character(Date), "%Y%m%d") - as.Date(as.character(inoc), "%Y%m%d")) %>% 
  inner_join(key) %>% 
 # filter(Day==14) %>% 
  mutate(element =factor(element, nutr_levels$element)) %>% 
  filter(element != "Na") %>% 
  mutate(Glc = ifelse(Glc == "0glc", "Autotrophic", "Glucose")) -> plot_icp

pval = plot_icp %>% select(element, concentration, Glc) %>% 
  group_by(element, Glc) %>% 
  mutate(concentration =list(concentration)) %>% 
  distinct() %>% 
  spread(Glc, concentration) %>% 
  mutate(Pval = t.test(unlist(Autotrophic), unlist(Glucose))$p.value,
         Diff =  t.test(unlist(Autotrophic), unlist(Glucose))$estimate[1]) %>% 
  mutate(Stars = case_when(Pval < 0.001 ~ "***",  Pval < 0.01 ~"**", Pval < 0.05~ "*")) %>% 
  select(-Autotrophic, -Glucose) %>% 
  inner_join(plot_icp %>% group_by(element) %>% filter(concentration== max(concentration)) %>% select(element, concentration) %>% mutate(y=concentration*2))
  
plot_icp %>% 
  ggplot(aes(x=element, y=concentration )) +
  #stat_summary(aes(fill=Glc, shape=Glc, group=interaction(element, Glc)), geom="errorbar", color="black", linetype="solid") +
  geom_jitter(aes(fill=Glc, shape=Glc), alpha=0.6, size=1, width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24) )+
  scale_fill_manual(values = c("Autotrophic"="#4ef542", "Glucose"="dark green") )+
  scale_y_log10()+
  scale_y_log10(breaks = c(1e-2, 1e-1,1e0,1e1,1e2,1e3,1e4),labels= c(0.01, 0.1,1,10,100,1000, 10000) ) +
  annotation_logticks(sides="l") +
  ylab(expression(paste("x", 10^7," atoms", " ",cell^-1 ))) +
  theme_classic(base_size=10)  +
  theme(legend.position = "none", axis.text.x = element_text(angle=45)) +
  geom_text(data=pval, mapping =aes(x=element, y, label=Stars)) +
  xlab(NULL)

ggsave("OptimizedMediaFigs/CellICP.pdf", device= "pdf", unit="in", width=4,height=3)

## Per volume data----

icp.by.vol <-read.csv("20190319_ICP_byvol_MgS.csv", header=T)
icp.by.vol$ID <- as.character(icp.by.vol$ID)
icp.by.vol$ID[1:3] <- c("02", "03", "04")
keyfile="20190225_ADJMgS_key.csv"
key <- read.csv(keyfile, header=TRUE)
key$ID <- as.character(key$ID)

library(ggsignif)
icp.by.vol %>% 
  gather(element, concentration, -ID, -Date ) %>%
  mutate(Day = as.Date(as.character(Date), "%Y%m%d") - as.Date(as.character(inoc), "%Y%m%d")) %>% 
  inner_join(key) %>% 
  mutate(element =factor(element, nutr_levels$element)) %>% 
  filter(element != "Na") %>% 
 # filter(element == "Fe") %>% 
  mutate(Glc = ifelse(Glc == "0glc", "Autotrophic", "Glucose")) -> plot_icp

#What a nice plotter! 
pval = plot_icp %>% select(element, concentration, Glc) %>% 
  group_by(element, Glc) %>% 
  mutate(concentration =list(concentration)) %>% 
  distinct() %>% 
  spread(Glc, concentration) %>% 
  mutate(Pval = t.test(unlist(Autotrophic), unlist(Glucose))$p.value,
         Diff =  t.test(unlist(Autotrophic), unlist(Glucose))$estimate[1]) %>% 
  mutate(Stars = case_when(Pval < 0.001 ~ "***",  Pval < 0.01 ~"**", Pval < 0.05~ "*")) %>% 
  select(-Autotrophic, -Glucose) %>% 
  inner_join(plot_icp %>% group_by(element) %>% filter(concentration== max(concentration)) %>% select(element, concentration) %>% mutate(y=concentration*2))
 
plot_icp %>% 
  ggplot(aes(x=element, y=concentration )) +
  #stat_summary(aes(fill=Glc, shape=Glc, group=interaction(element, Glc)), geom="errorbar", color="black", linetype="solid") +
  geom_jitter(aes(fill=Glc, shape=Glc), alpha=0.8, size=2, width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24) )+
  scale_fill_manual(values = c("Autotrophic"="#4ef542", "Glucose"="dark green") )+
  scale_y_log10(breaks = c(1e-4, 1e-3,1e-2,1e-1,1e0,1e1,1e2),labels= c(0.01, 0.1,1,10,100,1000, 10000) ) +
  annotation_logticks(sides="l") +
  ylab(expression(paste("x ", 10^5," atoms", " ",µm^-3 ))) +
  theme_classic(base_size=12)  +
  xlab(NULL) +
  theme(legend.position = "none", axis.text.x = element_text(angle=45)) +
  geom_text(data=pval, mapping =aes(x=element, y, label=Stars))
ggsave("OptimizedMediaFigs/Vol_ICP.pdf", device= "pdf", unit="in", width=8,height=3.5)


# Elemental Ratio-----
# This dataset creates different media compositions based on the scalar ratio 
cult.dat <- read.csv(file= "20190326_RatioGradient_GrowthData.csv", header=T)

media <- read.csv("20190415_MediaConcentrations.csv", header=T)
med.ord =  c("Kropat","ADJ.2xMgS" , "P14"  ,      "M10" ,   "M10.Ca" ,     "M14"  )     

media =media %>% select(-Bristols, -NewMedia)
#This isn'necessary 
cult.dat$Media <-factor(cult.dat$Media, levels = med.ord)
man.col2 <- c("dark green", "orange", "light blue", "blue", "red"  )

## Plot the MEDIA comp ----
media %>% 
  gather(Media, concentration, -element) %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>%
  filter(!element %in% c("Na", "Mo","Se")) %>% 
  ggplot(aes(x=element, y= concentration, color=Media)) +
  geom_point(alpha=0.9,  size=2.5, shape=15) +
  scale_color_manual(values= man.col2) +
  scale_y_log10() +
  annotation_logticks(sides="l") +
  ylab(expression(paste("Adapted media concentrations(µM)"))) +
  theme_classic(base_size=10)  +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle=45), legend.position = "none")

ggsave("OptimizedMediaFigs/Media Ratio Gradient_NoLegend.pdf", device= "pdf", unit="in", width=5,height=3)

##  Growth ----
### Cell Density ----
cult.dat %>% group_by(Day,Media, Glc ) %>% dplyr::summarise(CellDensity.mn=mean(CellDensity), CellDensity.SE=standard.error(CellDensity),
                                                            TotVol.mn=mean(TotVol), TotVol.SE=standard.error(TotVol)) %>%  
  ggplot(aes(x=Day, y=CellDensity.mn, color=Media, linetype=Glc))  +
  geom_errorbar(aes(ymin=CellDensity.mn-CellDensity.SE, ymax=CellDensity.mn+CellDensity.SE), width=0.1, size=0.5) +geom_line(size=1.2)  + theme_classic()+geom_line() + geom_point() +
  scale_y_log10(breaks=c(1E5, 1E6,1E7,1E8)) + annotation_logticks(sides="l") +
  ylab(expression(paste("Cells ", mL^{-1}, " ", log[10])))+
  scale_x_continuous(breaks=0:40, labels=0:40) +
  ylab(expression(paste("Cells ", mL^{-1}, " "))) +scale_color_manual(values = man.col2)

### Volumetric Biomass ---- 
cult.dat %>% 
  mutate(Glc = ifelse(grepl("0",Glc), "Autotrophic", "Glucose")) %>% 
  group_by(Day,Media, Glc) %>% dplyr::summarise(TotVol.mn=mean(TotVol), TotVol.SE=standard.error(TotVol),
                                                           TotVol.mn=mean(TotVol), TotVol.SE=standard.error(TotVol)) %>% 
  ggplot(aes(x=Day, y=TotVol.mn, fill=Media,color=Media, shape=Glc, linetype=Glc)) + 
  geom_errorbar(aes(ymin=TotVol.mn-TotVol.SE, ymax=TotVol.mn+TotVol.SE), width=0.5, color="black", linetype="solid")+
  geom_line()  + 
  theme_classic()+
  geom_point(color= "black") +
  scale_y_log10(breaks=c(1E5, 1E6,1E7,1E8,1E9,1E10,1E11,1E12,1E13), 
                 labels=c(1,10,100,1000,10000,100000,1000000,10000000,100000000 )/10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x", 10^9," ", mu*m^{3}, " ", mL^{-1} ))) +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  scale_color_manual(values = man.col2)+
  scale_fill_manual(values = man.col2)+
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
   theme_classic(base_size=10) +
   theme( legend.position = "none",
         axis.text.x = element_text(angle=45))
 
ggsave("OptimizedMediaFigs/RatioGradient GrowthCurve_NoLegend.pdf", device= "pdf", unit="in", width=4,height=3)

## Doubling time calculation ----
start  <- (cult.dat %>% filter(Day==0) %>% group_by(Day) %>% dplyr::summarise(mean=mean(TotVol)))$mean
photo.td <- cult.dat %>% filter(Glc=="0Glc")  %>% filter(Day==0 | Day==5) %>% #  group_by(ID, Media, Glc, Nutrient.Boost) %>%
  dplyr::mutate(mu=(log(TotVol)-log(start))/5) %>% mutate(td = (log(2)/mu)*24) %>%filter(Day==5) %>%
  select(Media, Glc, td)# %>% group_by(Media, Glc) %>% dplyr::summarise(mean=mean(td), se=standard.error(td))

#Calculating doubling time. Use the inner_join function so that you can substract each by the average of the TimePoint 0 point
td.key <- cult.dat %>% filter(Glc=="Glc")  %>% filter(Day==4) %>% group_by( Media) %>% dplyr::summarize(start=mean(TotVol)) 
mixo.td <- cult.dat %>% filter(Glc=="Glc")  %>% filter(Day==7) %>% inner_join(td.key) %>% #  group_by(ID, Media, Glc, Nutrient.Boost) %>%
  dplyr::mutate(  mu=(log(TotVol)-log(start))/3) %>% mutate(td = (log(2)/mu)*24) %>%
  select(Media, Glc, td) #%>%group_by(Media, Glc) %>% dplyr::summarise(mean=mean(td), se=standard.error(td))


bind_rows(photo.td, mixo.td) %>% 
  mutate(Glc = ifelse(grepl("0",Glc), "Autotrophic", "Glucose")) %>% 
 #filter(Glc == "Glucose") %>% 
  ggplot() + 
  stat_summary(aes(x=Media, y= td, fill=Media), color="black", geom = "bar") +
  stat_summary(aes(x=Media, y= td, fill=Media), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Media, y= td, fill=Media, shape=Glc), color="black", width=0.1) +
  facet_wrap(~Glc) + theme_classic(base_size=10) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  scale_y_continuous(breaks=c(0,6,12, 18,24)) + scale_fill_manual(values=man.col2) +
  ylab("Doubling time (h)") + xlab(NULL) +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
ggsave("OptimizedMediaFigs/RatioGradientDoublingTime_NoLegend.pdf", device= "pdf", unit="in", width=4,height=3)


# ICP-MS : Day Ratio -----
##ICP-MS data of the Ratio gradient 

#load all files 
setwd("~/Documents/Czof_Media_Optimization/")
icp.key <- read.csv( "RatioGradientICPMS_key.csv", header=T)
icp.pellet <-read.csv( "20190326_ICP_Pellet.csv", header=T)
icp.spent <- read.csv( "20190326_RatioGradient_ICP_Spent.csv", header=T)
icp.s.norm <- read.csv("20190326_RatioGradient_SulfNormICP.csv" , header=T)
media <- read.csv("~/Documents/20190415_MediaConcentrations.csv", header=T)


##Spent ----
icp.spent %>% gather(element, value, -ID) %>% 
  inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID") %>% 
  filter(grepl("Fresh", ID)) %>% rename(Fresh=value) %>% select(-Day, -Glc, -ID) %>%
  inner_join(icp.spent %>% gather(element, value, -ID) %>% 
               inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID")) %>% 
  #inner_join(media2 %>% rename(Media=media)) %>% 
  mutate(Percentage = 100*(value/Fresh)) %>% 
  #filter(Day==4) %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>%
  filter(!element %in% c("Na", "Mo","Se")) %>% 
  filter(element %in%  c("S", "P")) %>% 
  filter(Day>4, Glc=="Glc") %>% 
  mutate(Glc = ifelse(grepl("0",Glc),"Autotrophic","Glucose")) %>% 
  ggplot() +
  geom_hline(yintercept = 100, color="blue")+
  stat_summary(aes(x=Media, y=Percentage, fill=Media), color="black", geom = "bar") +
  stat_summary(aes(x=Media, y=Percentage, fill=Media), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Media, y= Percentage, fill=Media, shape=Glc), color="black", width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  theme_classic(base_size=10) +
  facet_grid(rows=vars(element), cols = vars(Glc)) +
  ylab("% of fresh media") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100, 110,120), labels=c(0,10,20,30,40,50,60,70,80,90,100,110,120)) +
  scale_fill_manual(values = man.col2) +xlab(NULL)+
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

ggsave("~/Documents/OptimizedMediaFigs/PandS_Spent.pdf", device= "pdf", unit="in", width=3,height=6)

#For reporting: 
icp.spent %>% gather(element, value, -ID) %>% 
  inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID") %>% 
  filter(grepl("Fresh", ID)) %>% rename(Fresh=value) %>% select(-Day, -Glc, -ID) %>%
  inner_join(icp.spent %>% gather(element, value, -ID) %>% 
               inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID")) %>% 
  #inner_join(media2 %>% rename(Media=media)) %>% 
  mutate(Percentage = 100*(value/Fresh)) %>% filter(Day==7, Glc=="Glc") %>%
  group_by(element, Media) %>% dplyr::summarise(Percentage=mean(Percentage)) %>% 
  spread(element, Percentage)


icp.pellet%>% 
  gather(element, value, -ID) %>%  
  inner_join(icp.key %>% select(ID, Media, Day, Glc,PelletVol,TotVol), by="ID") %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>%
  mutate(value=value/PelletVol) %>% 
  mutate(intMol = (value*TotVol*1E7*1000*1E6)/6E23) %>% 
  #inner_join(media2 ) %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>% 
  filter(Day>4, Glc=="Glc") %>% 
  filter(element %in%  c("S", "P")) %>% 
  mutate(Glc = ifelse(grepl("0",Glc),"Autotrophic","Glucose")) %>% 
  ggplot() +
  stat_summary(aes(x=Media, y=value*100, fill=Media), color="black", geom = "bar") +
  stat_summary(aes(x=Media, y=value*100, fill=Media), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Media, y=value*100, fill=Media, shape=Glc), color="black", width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  theme_classic(base_size=10) +
  facet_grid(rows=vars(element), cols = vars(Glc), scale="free_y") +
  scale_fill_manual(values = man.col2) + xlab(NULL)+
 # scale_y_continuous(breaks = c(0,1,2,3,4,5),labels= c(0,1,2,3,4,5)*100 ) +
  ylab(expression(paste("x ", 10^5," atoms", " ",µm^-3 )))  +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))


ggsave("OptimizedMediaFigs/PandS_Pellet.pdf", device= "pdf", unit="in", width=3,height=6)


  ##Mg and Mn separately-----
icp.spent %>% gather(element, value, -ID) %>% 
  inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID") %>% 
  filter(grepl("Fresh", ID)) %>% rename(Fresh=value) %>% select(-Day, -Glc, -ID) %>%
  inner_join(icp.spent %>% gather(element, value, -ID) %>% 
               inner_join(icp.key %>% select(ID, Media, Day, Glc), by="ID")) %>% 
  #inner_join(media2 %>% rename(Media=media)) %>% 
  mutate(Percentage = 100*(value/Fresh)) %>% 
  #filter(Day==4) %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>%
  filter(!element %in% c("Na", "Mo","Se")) %>% 
  filter(element %in%  c("Mg", "Mn")) %>% 
  filter(Day>4, Glc=="Glc") %>% 
  mutate(Glc = ifelse(grepl("0",Glc),"Autotrophic","Glucose")) %>% 
  ggplot() +
  geom_hline(yintercept = 100, color="blue")+
  stat_summary(aes(x=Media, y=Percentage, fill=Media), color="black", geom = "bar") +
  stat_summary(aes(x=Media, y=Percentage, fill=Media), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Media, y= Percentage, fill=Media, shape=Glc), color="black", width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  theme_classic(base_size=10) +
  facet_grid(rows=vars(element), cols = vars(Glc)) +
  ylab("% of fresh media") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100, 110,120), labels=c(0,10,20,30,40,50,60,70,80,90,100,110,120)) +
  scale_fill_manual(values = man.col2) +xlab(NULL)+
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))
ggsave("OptimizedMediaFigs/MnandMg_Spent.pdf", device= "pdf", unit="in", width=3,height=6)

icp.pellet%>% 
  gather(element, value, -ID) %>%  
  inner_join(icp.key %>% select(ID, Media, Day, Glc,PelletVol,TotVol), by="ID") %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>%
  mutate(value=value/PelletVol) %>% 
  mutate(intMol = (value*TotVol*1E7*1000*1E6)/6E23) %>% 
  #inner_join(media2 ) %>% 
  mutate(Media = factor(Media, levels= med.ord)) %>% 
  mutate(element =factor(element, c("N", nutr_levels$element))) %>% 
  filter(Day>4) %>% 
  filter(element %in%  c("Mg", "Mn")) %>% 
  mutate(Glc = ifelse(grepl("0",Glc),"Autotrophic","Glucose")) %>% 
  mutate(value = value*100) %>%  #make 10 to the 5
  ggplot() +
  stat_summary(aes(x=Media, y=value, fill=Media), color="black", geom = "bar") +
  stat_summary(aes(x=Media, y=value, fill=Media), color="black", geom = "errorbar", width=0.5) +
  geom_jitter(aes(x=Media, y= value, fill=Media, shape=Glc), color="black", width=0.1) +
  scale_shape_manual(values = c("Autotrophic"=21, "Glucose"=24)) +
  theme_classic(base_size=10) +
  facet_grid(rows=vars(element), cols = vars(Glc), scale="free_y") +
  scale_fill_manual(values = man.col2) + xlab(NULL)+
  ylab(expression(paste("x ", 10^5," atoms", " ",µm^-3 )))  +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))


ggsave("OptimizedMediaFigs/MnMg_Pellet.pdf", device= "pdf", unit="in", width=5,height=6)





# Concentration IMPACT-----
  #A) collect all day ratio samples 
#Ratio Gradient
  icp.key <- read.csv( "RatioGradientICPMS_key.csv", header=T)
  icp.s.norm <- read.csv("20190326_RatioGradient_SulfNormICP.csv" , header=T)
#####plotting all the media adjustments-> pick a color scale! 
media <- read.csv("20190415_MediaConcentrations.csv", header=T)
med.ord = colnames(media)[2:length(colnames(media))]
media2 <- media %>% gather(Media, "concentration", -element) 
#media.dr = read.csv("~/Documents/DAYRATIO_CONCENTRATION.csv",header=T) %>%   gather(element, concentration, -Media) 
  
#MEDIA_TOT = full_join(media.dr,media2 )

MEDIA_TOT = media2
#ADJ  
icp.by.cell <-read.csv("~/Documents/20190319_ICP_bycell_MgS.csv", header=T)
icp.by.cell$ID <- as.character(icp.by.cell$ID)
icp.by.cell$ID[1:3] <- c("02", "03", "04")
keyfile="20190225_ADJMgS_key.csv"
icp.adj.key <- read.csv(keyfile, header=TRUE)
icp.adj.key$ID <- as.character(icp.adj.key$ID)
  inoc <- "20190225" 
#Get this normalized
icp.adj.snorm =  icp.by.cell %>% 
  gather(element, concentration,-ID, -Date, -S) %>% 
  mutate(Snorm = (concentration/S)*1000) %>% 
  inner_join(icp.adj.key) %>% 
  filter(Glc == "0glc") %>% 
  mutate(Day = as.Date(as.character(Date), "%Y%m%d") - as.Date(inoc, "%Y%m%d"))


#MediaRatio = 
icp.rg.key <- read.csv( "RatioGradientICPMS_key.csv", header=T)
icp.rg.norm1 <- read.csv("20190326_RatioGradient_SulfNormICP.csv" , header=T)

icp.rg.norm = icp.rg.norm1 %>% 
  gather(element, Snorm, -ID, -S) %>% 
  inner_join(icp.rg.key) %>% 
  filter(Glc == "0Glc") 

## Day Ratio
icp.dr.snorm1 = read.csv("20190515_DR_ICPSULFNORM.csv", header=T)
icp.dr.key = read.csv("20190514_DR_Key.csv")

icp.dr.snorm = icp.dr.snorm1 %>% 
  gather(element, Snorm, -ID, -S) %>% 
  inner_join(icp.dr.key) %>% 
  mutate(Glc = "0Glc") 

o

COMBINED_ICP= full_join(icp.rg.norm %>% filter(Glc == "0Glc"), icp.dr.snorm ) %>% 
    
    select(element, Snorm,Media, Day) %>% 
  full_join(icp.adj.snorm  %>% 
  mutate(Media = case_when(Media == "adj" ~"ADJ",
                           Media == "adj2x" ~"ADJ.2xMgS")  
            ) %>% select(element, Snorm,Media, Day) %>% mutate(Day=as.numeric(Day))) %>% 
  inner_join(MEDIA_TOT)

## Observe Results 
COMBINED_ICP %>% 
  mutate(element =factor(element, nutr_levels$element)) %>% 
  filter(!element %in% c("Na", "Mo", "Se")) %>% 
  filter(element %in% c("Fe", "Mn","Zn", "Cu")) %>%
  ggplot(aes(x=concentration, y=Snorm, color=Day)) +
  geom_point(alpha=0.5)+
  facet_wrap(~element, scales="free") +
  facet_grid(cols =vars(element), scales = "free_x") +
  scale_y_log10() + scale_x_log10() +
  theme_classic(base_size=10) +
  stat_smooth(aes(x=concentration, y=Snorm), method="lm", color="black")+
  scale_color_viridis_c()
 
COMBINED_ICP %>% 
  mutate(element =factor(element, nutr_levels$element)) %>% 
  filter(!element %in% c("Na", "Mo", "Se")) %>% 
  #filter(element %in% c("Fe", "Mn","Zn", "Cu")) %>%
  #filter(element %in% c( "Fe")) %>%
  
  ggplot(aes(x=concentration, y=Snorm, color=Day)) +
  geom_point(alpha=0.5)+
  facet_wrap(~element, scales="free") +
  facet_grid(cols =vars(element), scales = "free_x") +
  scale_y_log10() + scale_x_log10() +
  theme_classic(base_size=10) +
  stat_smooth(aes(x=concentration, y=Snorm), method="lm", color="black")+
  scale_color_viridis_c()



# GET EQUATION AND R-SQUARED AS STRING
# Inspiration: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

##Set up a means to find the pearson 
setup = COMBINED_ICP %>% filter(element == "Fe", Day==7 | Day ==4 )
cor(setup$Snorm, setup$concentration)

m = lm(Snorm ~ concentration, COMBINED_ICP %>% filter(element == "Fe") %>% filter(Day ==14))
sm = summary(m) 
sm$r.squared
  day_order = COMBINED_ICP %>% distinct(Day) %>% pull() %>% as.character()

library(viridis)
COMBINED_ICP %>% 
  mutate(Day = factor(as.character(Day), levels = day_order)) %>% 
  group_by(element) %>% mutate(RelativeCon = concentration/max(concentration)) %>% 
  mutate(element =factor(element, nutr_levels$element)) %>% 
  filter(!element %in% c("Na", "Mo", "Se")) -> COMBINED_ICP_1

COMBINED_ICP_2 = COMBINED_ICP_1 %>% filter(!(element == "Mn"& Media=="P14" ))
 
#Test run, not ordered byfactor  
COMBINED_ICP %>% 
  #filter(Day %in% c(7)) %>% 
  ggplot() +
  stat_smooth(data=COMBINED_ICP_2, aes(x=concentration, y=Snorm), method="lm", color="black")+
  geom_jitter(aes(x=concentration, y=Snorm, fill=as.character(Day)),alpha=0.5, size=1.5,pch=21, color="black", width=0.1)+
  #facet_wrap(~element, scales="free") +
  facet_grid(cols =vars(element), scales="free_x" ) +
  scale_y_log10() + 
  scale_x_log10(breaks =c(1,10,100,1000,10000,100000),  labels = c(1,10,100,1000,10000,100000)) +
  theme_bw(base_size=10) +
  scale_fill_manual(values = viridis(5) ) +
  ylab(expression(paste("mmol ", mol^-1, " S"))) +
  annotation_logticks(sides = "bl") +
  xlab("Media concentration (µM)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1), strip.background = element_rect(fill="white"))

### Final Ordered For Figure ----
#use factors
factor_order = c("P", "K", "Mg","Fe", "Ca", "Zn", "Mn", "Cu")
COMBINED_ICP$element = factor(COMBINED_ICP$element , levels = factor_order)
COMBINED_ICP_2$element = factor(COMBINED_ICP_2$element , levels = factor_order)

##Plot elements of importance
COMBINED_ICP %>% 
  filter(element %in% factor_order) %>% 
  #filter(Day %in% c(7)) %>% 
  ggplot() +
  stat_smooth(data=COMBINED_ICP_2 %>%   filter(element %in% factor_order), aes(x=concentration, y=Snorm), method="lm", color="black")+
  geom_jitter(aes(x=concentration, y=Snorm, fill=as.character(Day)),alpha=0.5, size=1.5,pch=21, color="black", width=0.1)+
  #facet_wrap(~element, scales="free") +
  facet_grid(cols =vars(element), scales="free_x" ) +
  scale_y_log10() + 
  scale_x_log10(breaks =c(1,10,100,1000,10000,100000),  labels = c(1,10,100,1000,10000,100000)) +
  theme_bw(base_size=10) +
  scale_fill_manual(values = viridis(5) ) +
  ylab(expression(paste("mmol ", mol^-1, " S"))) +
  annotation_logticks(sides = "bl") +
  xlab("Media concentration (µM)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1), strip.background = element_rect(fill="white"))
ggsave("OptimizedMediaFigs/SNORM_Concentration_Day_nolegend.pdf", device= "pdf", unit="in", width=8,height=4)




# N source ----

#Create a directory in your computer of the experiment name in your working directoyr  
#This has to be reset due to Loop have a setwd function internally
setwd("~/Documents/Czof_Media_Optimization/")

path="20190805_NSource/";x="20190805_NSourceFinal";inoc="20190805";keyfile="20190805_NSource_key.csv"   #the working directoy you are working with
nsource = LoopCount(path, x, inoc, keyfile)

## Zero connector 
nsource[[3]] %>% 
  mutate(int =interaction(as.character(Concentration),N) ) %>% pull(int) %>% unique()

nsource[[3]] %>% #distinct(Concentration)
 # filter(N == "NH4") %>% 
ggplot(aes(x=Day, y=TotVol, group=interaction(as.character(Concentration),N),
           fill=interaction(as.character(Concentration),N),color=interaction(as.character(Concentration),N),
        shape = as.character(Concentration))) + 
  stat_summary(geom="errorbar", linetype="solid") +
  stat_summary(geom="line") +
  geom_jitter( alpha=0.8, size=1.5, width=0.2, color="black") +
  scale_fill_manual(values = c("#FFE178", "#FFA800", "#BF7E00", "#78CEFF", "#0075FF", "#001AA1"))  +
  scale_color_manual(values = c("#FFE178", "#FFA800", "#BF7E00", "#78CEFF", "#0075FF", "#001AA1"))+

#  scale_fill_manual(values = c("NO3"="blue", "NH4"="orange") )+
  #scale_color_manual(values = c("NO3"="blue", "NH4"="orange"))+
  scale_shape_manual(values = c("11.25"=21,"45"=23, "22.5"=24) )+
  theme_classic(base_size=10)+ # geom_point( color="dark green") +
  scale_y_log10(breaks=c(1E5, 1E6,1E7,1E8,1E9,1E10,1E11,1E12,1E13), 
                labels=c(1,10,100,1000,10000,100000,1000000,10000000,100000000 )/10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9," ", mu*m^{3}, " ", mL^{-1} ))) +
  theme(#legend.position = "none", 
        axis.text.x = element_text(angle=45))

## Nsource code for the above, we can move along smoothly for this  
nsource[[3]] %>% 
  ggplot(aes(x=Day, y=TotVol, group=interaction(Concentration, N, sep = " mM "),
             fill=interaction(Concentration, N, sep = " mM "),
             color=interaction(Concentration, N, sep = " mM "),
             shape=as.character(Concentration))) + 
  stat_summary(geom="errorbar", linetype="solid") +
  stat_summary(geom="line") +
  geom_jitter(alpha=0.8, size=1.5, width=0.2, color="black") +
  scale_fill_manual(name = "N treatment",
                    values = c("#FFE178", "#FFA800", "#BF7E00", "#78CEFF", "#0075FF", "#001AA1"),
                    guide = guide_legend(override.aes = list(shape = c(21, 24, 23)))) +
  scale_color_manual(name = "N treatment",
                     values = c("#FFE178", "#FFA800", "#BF7E00", "#78CEFF", "#0075FF", "#001AA1"),
                     guide = guide_legend(override.aes = list(shape = c(21, 24, 23)))) +
  scale_shape_manual(guide = "none",
                     values = c("11.25" = 21, "45" = 23, "22.5" = 24)) +
  theme_classic(base_size=10) +
  scale_y_log10(breaks=c(1E5, 1E6, 1E7, 1E8, 1E9, 1E10, 1E11, 1E12, 1E13), 
                labels=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000) / 10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
  theme(axis.text.x = element_text(angle=45))

setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/NsourceTreatment.pdf", device= "pdf", unit="in", width=5,height=3)


# pH ----
## Check ALL loops to run everything ----
setwd("~/Documents/Czof_Media_Optimization/")
path="pH_Curves";x="20190513 pH grad II";inoc="20190513";keyfile="20190513_pHIIgrad.csv"
ph1 = LoopCount(path, x, inoc, keyfile)

setwd("~/Documents/Czof_Media_Optimization/")
path="pH_Curves/";x="20190506_pHGradient";inoc="20190506";keyfile="20190506_phGrad_key.csv"
ph2 =  LoopCount(path, x, inoc, keyfile)
## pH Gradient Growth ---- 

### Create the Zero data representatives for all values ----
## Then Bind 2 experiments 
BothPh_Reps = ZeroProcess_Data(ph1[[3]], pH) %>% 
  mutate(Exp = "Exp. Rep 1") %>% 
  bind_rows(ZeroProcess_Data(ph2[[3]], pH) %>% 
            mutate(Exp="Exp. Rep 2"))
  
## Log Plotting of data  
BothPh_Reps %>% 
  ggplot(aes(x=Day, y=TotVol, group=interaction(as.character(pH),Exp),
             fill=as.character(pH),
             color=as.character(pH))) + 
  stat_summary(geom="errorbar", linetype="solid") +
  stat_summary(geom="line") +
  geom_jitter(alpha=0.8, size=1.5, width=0.2, color="black", shape=21) +
  theme_classic(base_size=10) +
  scale_y_log10(breaks=c(1E5, 1E6, 1E7, 1E8, 1E9, 1E10, 1E11, 1E12, 1E13), 
                labels=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000) / 10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
 labs(color = "pH", fill = "pH")+
  theme(axis.text.x = element_text(angle=45)) +
  facet_grid(cols = vars(Exp)) +
  theme(legend.position = "none")
setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/LogpHGradient.pdf", device= "pdf", unit="in", width=3.5,height=3.5)

##
p_both = BothPh_Reps %>% 
  ggplot(aes(x=Day, y=TotVol, group=interaction(as.character(pH),Exp),
             fill=as.character(pH),
             color=as.character(pH))) + 
  stat_summary(geom="errorbar", linetype="solid") +
  stat_summary(geom="line") +
  geom_jitter(alpha=0.8, size=1.5, width=0.2, color="black", shape=21) +
  theme_classic(base_size=10) +
  scale_y_continuous(breaks=c(0, 3E8,6E8, 9E8, 1.2E9), 
                labels=c(0,0.3,0.6,0.9,1.2 )) +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
  labs(color = "pH", fill = "pH")+
  theme(axis.text.x = element_text(angle=45)) +
  facet_grid(cols = vars(Exp)) 
 #Plot the values 
p_both

setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/UnlogpHGradient.pdf", device= "pdf", unit="in", width=4.5,height=3.5)


# Buffer Information ---- 

# input data  
pH.3 <- read.csv("20190528_buffers_culturedata.csv", header=T)%>% mutate(Experiment = "pH.3") %>% filter(!ID %in% c("buff10", "buff13"))


## Volumetric Biomass ----- 
pH.3 %>% 
  mutate(pH = paste("pH ", pH)) %>% 
  ggplot(aes(x = Day, y = TotVol, 
             group = interaction(as.character(pH), Buffer, sep =" "), 
             fill = interaction(as.character(pH), Buffer, sep =" "), 
             color = interaction(as.character(pH), Buffer, sep =" "), 
             shape = interaction(as.character(pH), Buffer, sep =" "), 
             linetype = interaction(as.character(pH), Buffer, sep =" "))) + 
  theme_classic(base_size = 10) +
  
  stat_summary(geom = "errorbar", linetype = "solid") +
  stat_summary(geom = "line") +
  geom_jitter(alpha = 0.8, size = 1.5, width = 0.2, color = "black") +
  
  # Set scales for shape, fill, color, and linetype, and integrate them into one legend
  scale_shape_manual(name = "Treatment",
                     values = c(21, 22, 21, 22, 21, 22)) +
  
  scale_fill_manual(name = "Treatment",
                    values = c("#C80815", "#C80815", "#800080", "#800080", "#CCAC00", "#CCAC00")) +
  
  scale_color_manual(name = "Treatment",
                     values = c("#C80815", "#C80815", "#800080", "#800080", "#CCAC00", "#CCAC00")) +
  
  scale_linetype_manual(name = "Treatment",
                        values = c("solid", "dashed", "solid", "dashed", "solid", "dashed")) +
  
  scale_y_log10(breaks=c(1E5, 1E6, 1E7, 1E8, 1E9, 1E10, 1E11, 1E12, 1E13),
                labels=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000) / 10000) +
  annotation_logticks(sides="l") +
  scale_x_continuous(breaks=0:40, labels=0:40) + 
  ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
  theme(axis.text.x = element_text(angle=45)) +
  theme(legend.position = "none") #There is no legend position, as it will be used in the next set 

#output
setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/Buffer_Tot_Volume.pdf", device= "pdf", unit="in", width=3.5,height=3.5)

## Cell density ---- 
## Tris buffer creates a growth depression

pH.3 %>% 
  mutate(pH = paste("pH ", pH)) %>% 
  ggplot(aes(x = Day, y = CellDensity, 
             group = interaction(as.character(pH), Buffer, sep =" "), 
             fill = interaction(as.character(pH), Buffer, sep =" "), 
             color = interaction(as.character(pH), Buffer, sep =" "), 
             shape = interaction(as.character(pH), Buffer, sep =" "), 
             linetype = interaction(as.character(pH), Buffer, sep =" "))) + 
  theme_classic(base_size = 10) +
  
  stat_summary(geom = "errorbar", linetype = "solid") +
  stat_summary(geom = "line") +
  geom_jitter(alpha = 0.8, size = 1.5, width = 0.2, color = "black") +
  
  # Set scales for shape, fill, color, and linetype, and integrate them into one legend
  scale_shape_manual(name = "Treatment",
                     values = c(21, 22, 21, 22, 21, 22)) +
  
  scale_fill_manual(name = "Treatment",
                    values = c("#C80815", "#C80815", "#800080", "#800080", "#CCAC00", "#CCAC00")) +
  
  scale_color_manual(name = "Treatment",
                     values = c("#C80815", "#C80815", "#800080", "#800080", "#CCAC00", "#CCAC00")) +
  
  scale_linetype_manual(name = "Treatment",
                        values = c("solid", "dashed", "solid", "dashed", "solid", "dashed")) +
  
  scale_y_log10(breaks = c(1E5, 1E6, 1E7, 1E8, 1E9), labels = c(1E5, 1E6, 1E7, 1E8, 1E9)/10^6) + 
  annotation_logticks(sides = "l") +
  scale_x_continuous(breaks = 0:40, labels = 0:40) + 
  ylab(expression(paste("x " , 10^6 , " cells ", mL^{-1}))) +
  theme(axis.text.x = element_text(angle = 45))

setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/Buffer_Cell_Density.pdf", device= "pdf", unit="in", width=4.5,height=3.5)


## Density Day Compare ----

pH.3 %>% 
  filter(Day == 7) -> data_day7

data_day7$Buffer -> x
data_day7$CellDensity -> y

# Perform pairwise t-tests
PairwiseT_Process = function(x, y){

pairwise_t_test_results <- pairwise.t.test(y, x, p.adjust.method = "bonferroni")
# Extract p-values and format them for plotting
 broom::tidy(pairwise_t_test_results) %>%
  mutate(sci.p.value = format(p.value, scientific = TRUE, digits = 3)) %>% 
   return()
}
## how to process the three results 
PairwiseT_Process(data_day7$Buffer,data_day7$TotVol )
PairwiseT_Process(data_day7$Buffer,data_day7$CellDensity )
PairwiseT_Process(data_day7$Buffer,data_day7$Mean )

## Significance 

# Each Day 7 signficance plot 
data_day7 %>% 
  ggplot(aes(x=Buffer, fill=Buffer, y=CellDensity))+
  stat_summary(geom = "bar", color="black") +
  geom_jitter(aes(shape=as.character(pH)), width=0.2) +
  
  scale_shape_manual(name = "pH",
                     values = c(21, 22)) +
  scale_fill_manual(name = "Buffer",
                    values = c("#C80815", "#800080", "#CCAC00")) +
  theme_classic(base_size = 10) +
    geom_signif(
      comparisons = list(c("control", "HEPES"), c("control", "TRIS"), c("HEPES", "TRIS")),
      map_signif_level = TRUE,
      step_increase = 0.1,
      tip_length = 0.02
    ) +
  xlab(NULL) +
  scale_y_continuous(breaks = c(0, 1E7, 2E7), labels =  c(0, 1E7, 2E7)/10^6)+
  ylab(expression(paste("x " , 10^6 , " cells ", mL^{-1}))) +
  theme(legend.position = "none")

ggsave("OptimizedMediaFigs/Buffer_Day7_Cell_Density.pdf", device= "pdf", unit="in", width=2,height=1.75)


data_day7 %>% 
  ggplot(aes(x=Buffer, fill=Buffer, y=Mean))+
  stat_summary(geom = "bar", color="black") +
  geom_jitter(aes(shape=as.character(pH)), width=0.2) +
  
  scale_shape_manual(name = "pH",
                     values = c(21, 22)) +
  scale_fill_manual(name = "Buffer",
                    values = c("#C80815", "#800080", "#CCAC00")) +
  theme_classic(base_size = 10) +
  geom_signif(
    comparisons = list(c("control", "HEPES"), c("control", "TRIS"), c("HEPES", "TRIS")),
    map_signif_level = TRUE,
    step_increase = 0.1,
    tip_length = 0.02
  ) +
  ylab(expression(paste(mu, "m"))) +
  xlab(NULL)  + theme(legend.position = "none")

ggsave("OptimizedMediaFigs/Buffer_Day7_CellSize.pdf", device= "pdf", unit="in", width=2,height=1.75)


data_day7 %>% 
  ggplot(aes(x=Buffer, fill=Buffer, y=TotVol))+
  stat_summary(geom = "bar", color="black") +
  geom_jitter(aes(shape=as.character(pH)), width=0.2) +
  
  scale_shape_manual(name = "pH",
                     values = c(21, 22)) +
  scale_fill_manual(name = "Buffer",
                    values = c("#C80815", "#800080", "#CCAC00")) +
  theme_classic(base_size = 10) +
  geom_signif(
    comparisons = list(c("control", "HEPES"), c("control", "TRIS"), c("HEPES", "TRIS")),
    map_signif_level = TRUE,
    step_increase = 0.1,
    tip_length = 0.02
  ) +
  scale_y_continuous(breaks = c(0, 5E8, 1e9), labels  =c(0, 5E8, 1e9)/10^9 ) +
  ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
  theme(legend.position = "none") + xlab(NULL)



ggsave("OptimizedMediaFigs/Buffer_Day7_TotVOL .pdf", device= "pdf", unit="in", width=2,height=1.75)

### pH in spent media----
#pH over time
setwd("~/Documents/Czof_Media_Optimization/")
pH <- read.csv("20190528_buffer_pH.csv", header=T)

##find the key conversion
key = pH.3 %>% distinct(ID, pH,Buffer)

#Gatherandfindthekeyfile 
ph_spent = pH %>% gather(Date, pH.spent, -ID) %>% mutate(Date= gsub("X", "", Date)) %>% filter(!ID %in% c("buff10", "buff13")) %>%
  mutate(Date= as.Date(as.character(Date), "%Y%m%d")) %>% mutate(Day= Date-as.Date("20190528", "%Y%m%d")) %>% 
  inner_join(key)

#Plot the value   
ph_spent %>% 
  filter(pH == 7.5) %>% 
  mutate(pH = paste("pH ", pH)) %>% 
  ggplot(aes(x = Day, y = pH.spent, 
             group =  Buffer, 
             fill = Buffer, 
             color = Buffer)) + 
  theme_classic(base_size = 10) +
  stat_summary(geom = "errorbar", linetype = "solid") +
  stat_summary(geom = "line") +
  geom_jitter(alpha = 0.8, size = 1.5, width = 0.2, color = "black",shape=21) +
  
  scale_fill_manual(name = "Buffer",
                    values = c("#C80815",  "#800080","#CCAC00")) +
  
  scale_color_manual(name = "Buffer",
                     values = c("#C80815",  "#800080","#CCAC00"))+ 
  scale_y_continuous(breaks =c(7.5, 8.0, 8.5), labels =c(7.5, 8.0, 8.5) ) +
  ylab("pH") +
  xlab("Time (d)") +
  scale_x_continuous(breaks = 0:40, labels = 0:40) + 
  theme(axis.text.x = element_text(angle = 45))
  
ggsave("OptimizedMediaFigs/ph_spent.pdf", device= "pdf", unit="in", width=4,height=3.5)

  
# K to NA -----
setwd("~/Documents/Czof_Media_Optimization/")
path="KtoNA";x="20190605 ratios";inoc="20190605";keyfile="20190605_KNaRatio_key.csv"
KNA = LoopCount(path, x, inoc, keyfile)
KNA[[3]] %>% #ungroup %>% distinct(pH, K.to.Na)
   filter(pH == 7.5) %>%
   mutate(K.to.Na = ifelse(K.to.Na== "1 to 1.6(M10)", "1 to 1.6 (CORE)", K.to.Na)) %>% 
   mutate(K.to.Na = gsub(" to ", ":", K.to.Na)) %>% 
     ggplot(aes(x=Day, y=TotVol, group= K.to.Na,
                fill= K.to.Na,
                color= K.to.Na)) + 
     stat_summary(geom="errorbar", linetype="solid") +
     stat_summary(geom="line") +
     geom_jitter(alpha=0.8, size=1.5, width=0.2, color="black", shape=21) +
     scale_fill_viridis_d(name = "K:Na Ratio") +
      scale_color_viridis_d(name = "K:Na Ratio") +
     theme_classic(base_size=10) +
     scale_y_log10(breaks=c(1E5, 1E6, 1E7, 1E8, 1E9, 1E10, 1E11, 1E12, 1E13),
                   labels=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000) / 10000) +
     annotation_logticks(sides="l") +
     scale_x_continuous(breaks=0:40, labels=0:40) + 
     ylab(expression(paste("x ", 10^9, " ", mu*m^{3}, " ", mL^{-1} ))) +
     theme(axis.text.x = element_text(angle=45)) +
     xlab("Time (d)")

setwd("~/Documents/Czof_Media_Optimization/")
ggsave("OptimizedMediaFigs/K_to_NA_Ratio.pdf", device= "pdf", unit="in", width=4,height=3.5)
 
 