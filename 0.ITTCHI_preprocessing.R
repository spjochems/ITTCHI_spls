library(xlsx)
library(ggplot2)
library(reshape2)


rm(list=ls()) #clear workspace

setwd('R://Para-CIH//CIH Group member folders//Simon//Results//ITTCHI mixomics//')

######################################################################
#Load in Clinical
######################################################################
Clinical <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 3)

#change the names in clinical
Clinical$Group[Clinical$Group == 1] <- 'Placebo'
Clinical$Group[Clinical$Group == 2] <- 'Intervention'
Clinical$Sex[Clinical$Sex == 1] <- 'Male'
Clinical$Sex[Clinical$Sex == 2] <- 'Female'

#check conversion
table(Clinical$Group)
table(Clinical$Sex)


######################################################################
#Load in the AE data
######################################################################
AE <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 4)
AE2 <- AE[,1:9]
colnames(AE2) <- gsub('..days.', '', colnames(AE2))
colnames(AE2)[9] <- 'Gr3_abdominal_AE_duration'
colnames(AE2) <- gsub('[.]', '_', colnames(AE2))

#check colnames
colnames(AE2)

#Convert the legends
AE2$Skin_AEs[AE2$Skin_AEs == 1] <- 'No'
AE2$Skin_AEs[AE2$Skin_AEs == 2] <- 'Yes'
AE2$Severe_abdominal_AE[AE2$Severe_abdominal_AE == 1] <- 'No'
AE2$Severe_abdominal_AE[AE2$Severe_abdominal_AE == 2] <- 'Yes'

table(AE2$Skin_AEs)
table(AE2$Severe_abdominal_AE)

#remove NA
sum(AE2 == -99)
AE2[AE2 == -99] <- NA
sum(AE2 == -99, na.rm=T)
sum(is.na(AE2))

#check the data
plot(as.factor(AE2$Skin_AEs), AE2$Rash_A0)
plot(as.factor(AE2$Skin_AEs), AE2$Rash_A03)
plot(as.factor(AE2$Skin_AEs), AE2$Rash_A06)
plot(as.factor(AE2$Skin_AEs), AE2$Rash_B0)
plot(as.factor(AE2$Skin_AEs), AE2$Rash_B02)
plot(as.factor(AE2$Severe_abdominal_AE), AE2$Gr3_abdominal_AE_duration)


#merge with clinical
All <- merge(Clinical, AE2, by = 'VolID')
nrow(All) == 19 #check all matched

All.melt <- melt(All[,1:11], id = colnames(All)[1:6])
ggplot(All.melt, aes(x=variable, y=value, group = VolID, colour = Group)) + 
    geom_line() + 
    stat_summary(aes(group = Group), fun = mean, geom = "line", size = 3)

######################################################################
#Egg data
######################################################################
Eggs <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 5)

#remove -99s
sum(Eggs == -99)
Eggs[Eggs == -99] <- NA
sum(Eggs == -99, na.rm=T)
sum(is.na(Eggs))

#check colnames
colnames(Eggs)

#Convert the legends
Eggs$Protected[Eggs$Protected == 1] <- 'No'
Eggs$Protected[Eggs$Protected == 2] <- 'Yes'
table(Eggs$Protected)

#check the data
plot(as.factor(Eggs$Protected), Eggs$Mean_egg_load)
plot(as.factor(Eggs$Protected), Eggs$Hatching_B12)
plot(as.factor(Eggs$Protected), Eggs$Hatching_B16)
plot(as.factor(Eggs$Protected), Eggs$Mean_Ct_wk_12.16)


#merge with all
All2 <- merge(All, Eggs, by = 'VolID')
nrow(All2) == 19 #check all matched

plot(as.factor(All2$Group), All2$Mean_egg_load)
plot(as.factor(All2$Group), All2$Gr3_abdominal_AE_duration)
plot(as.factor(All2$Protected), All2$Rash_A06)
plot(as.factor(All2$Protected), All2$Rash_A03)
plot(as.factor(All2$Protected), All2$Rash_A0)
plot(as.factor(All2$Protected), All2$Rash_B0)
plot(as.factor(All2$Protected), All2$Rash_B02)
ggplot(All2, aes(x=Gr3_abdominal_AE_duration, y = Mean_egg_load, colour = Skin_AEs))+
    geom_point() + facet_grid(.~Group)



######################################################################
#Load in Eosinophils
######################################################################
Eosinophils <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 6)

#remove -99s
sum(Eosinophils == -99)
Eosinophils[Eosinophils == -99] <- NA
sum(Eosinophils == -99, na.rm=T)
sum(is.na(Eosinophils))

Eosinophils[,3] <- as.numeric(Eosinophils[,3])
hist(Eosinophils[,3])

#make to long format
Eosinophils2 <- dcast(Eosinophils, VolID ~ Timepoint)
colnames(Eosinophils2)[2:ncol(Eosinophils2)] <- paste0('Eosinophils_', 
                      colnames(Eosinophils2)[2:ncol(Eosinophils2)])

colnames(Eosinophils2)


#Merge
All3 <- merge(All2, Eosinophils2, by = 'VolID')
nrow(All3) == 19 #check all matched

#checks
All3.melt <- melt(All3, id = colnames(All3)[1:17])
ggplot(All3.melt, aes(x=variable, y=value, group = VolID, colour = Group)) + 
    geom_line()

ggplot(All3.melt, aes(x=variable, y=value, group = VolID, 
                      colour = Protected)) +  geom_line()  +
    facet_grid(.~Group)+ geom_vline(xintercept = 8, linetype = 'dashed')+
    stat_summary(aes(group = Protected), fun = mean, geom = "line", size = 2)

ggplot(All3.melt, aes(x=variable, y=value, group = VolID, 
                      colour = Severe_abdominal_AE)) +  geom_line()  +
    facet_grid(.~Group)+ geom_vline(xintercept = 8, linetype = 'dashed')+
    stat_summary(aes(group = Severe_abdominal_AE), fun = mean, geom = "line", size = 2)

ggplot(All3.melt, aes(x=variable, y=value, group = VolID, 
                      colour = Skin_AEs)) +  geom_line()  +
    facet_grid(.~Group) + geom_vline(xintercept = 8, linetype = 'dashed')+
    stat_summary(aes(group = Skin_AEs), fun = mean, geom = "line", size = 2)




######################################################################
#Load in Luminex
######################################################################
Luminex <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 7)

#fix colnames
names <- strsplit(colnames(Luminex), split = '[..]')
colnames(Luminex) <- do.call(c, lapply(names, function(x) x[1]))

#remove cytokines that are not really found
co = 0.75 #cutoff to remove cytokines, all cytokines with less than this number measurements are discarded
keep <- colSums(Luminex == 'OOR <', na.rm = T) / nrow(Luminex) < co
Luminex2 <- Luminex[,keep]

colSums(Luminex2 == 'OOR <', na.rm = T) / nrow(Luminex2) #checl it worked

#remove the sample with missing values
#Luminex3 <- Luminex2[rowSums(is.na(Luminex2)) < 5,]
Luminex3 <- Luminex2

#make the data good for use
for(i in 3:ncol(Luminex3)){
    Luminex3[,i] <- gsub(',', '.', Luminex3[,i])
    Luminex3[,i] <- gsub('[*]', '', Luminex3[,i])
    low <- min(as.numeric(Luminex3[,i]), na.rm=T)
    high <- max(as.numeric(Luminex3[,i]), na.rm=T)
    Luminex3[,i][Luminex3[,i] == 'OOR <'] <- low
    Luminex3[,i][Luminex3[,i] == 'OOR >'] <- high
    Luminex3[,i]  <-   as.numeric(Luminex3[,i])
    print(colnames(Luminex3)[i])
    print(sum(Luminex3[,i] == low, na.rm=T))
    print(sum(Luminex3[,i] < low, na.rm=T))
    print(sum(Luminex3[,i] == high, na.rm=T))
    print( sum(Luminex3[,i] > high, na.rm=T))
}

#change timepoint names to match others
Luminex3$Timepoint <- gsub('00', '0', Luminex3$Timepoint)
Luminex3$Timepoint <- gsub('-', '', Luminex3$Timepoint)
Luminex3$Timepoint 

#exploratory plots
Lum_clin <- merge(All3, Luminex3, by = 'VolID')

for(i in 45:ncol(Lum_clin)){
    Lum_clin2 <- Lum_clin[,c(1:44, i)]
    name <-  colnames(Lum_clin2)[45]
    colnames(Lum_clin2)[45] <- 'value'
    
    plot1 <- ggplot(Lum_clin2, aes(x=Timepoint, y=value, group = VolID, 
                          colour = Protected)) +  geom_line()  +
        facet_grid(.~Group)+ ggtitle(name)+ 
        geom_vline(xintercept = 5, linetype = 'dashed')+
        stat_summary(aes(group = Protected), fun = mean, geom = "line", size = 2)
    print(plot1)
    
}

#make to long format
Luminex4 <- split(Luminex3, f = Luminex3$Timepoint)
lapply(Luminex4, function(x) print(nrow(x)))
Luminex5 <- lapply(Luminex4, function(x){
    a <- x[,-1]
    colnames(a)[2:ncol(a)] <- paste(colnames(a)[2:ncol(a)], x[1,1], sep = '_')
    colnames(a) <- gsub('Serum_', '',colnames(a))
    return(a)})

#divide all by the baselines and do a log transform
Luminex5.FC <- lapply(Luminex5, function(x) log2(x/Luminex5[[1]]))

Luminex6 <- do.call(cbind, Luminex5.FC)

Luminex6[,grepl('VolID', colnames(Luminex6))]
del <- which(grepl('VolID', colnames(Luminex6)))
Luminex6$VolID <- Luminex5[[1]]$VolID

Luminex7 <- Luminex6[,-del] 
colnames(Luminex7) <- do.call(c, lapply(strsplit(colnames(Luminex7), '[.]'), function(i) i[2]))
colnames(Luminex7)[ncol(Luminex7)] <- 'VolID'


#check it worked
Luminex7.melt <- melt(Luminex7, id = 'VolID')
Luminex7.melt$Timepoint <- do.call(c, lapply(strsplit(as.character(Luminex7.melt$variable), '_'), function(i) i[length(i)]))
Luminex7.melt$Cytokine <- as.character(Luminex7.melt$variable)  
Luminex7.melt$Cytokine <- do.call(c, lapply(strsplit(as.character(Luminex7.melt$Cytokine), '_A0'), function(i) i[1]))
Luminex7.melt$Cytokine <- do.call(c, lapply(strsplit(as.character(Luminex7.melt$Cytokine), '_B'), function(i) i[1]))
table(Luminex7.melt$Cytokine)
Luminex7.cast <- dcast(Luminex7.melt, VolID + Timepoint  ~ Cytokine, value.var = 'value')
LumFC_clin <- merge(All3, Luminex7.cast, by = 'VolID')

for(i in 45:ncol(LumFC_clin)){
    Lum_clin2 <- LumFC_clin[,c(1:44, i)]
    name <-  colnames(Lum_clin2)[45]
    colnames(Lum_clin2)[45] <- 'value'
    
    plot1 <- ggplot(Lum_clin2, aes(x=Timepoint, y=value, group = VolID, 
                                   colour = Protected)) +  geom_line()  +
        facet_grid(.~Group)+ ggtitle(name)+ 
        geom_vline(xintercept = 5, linetype = 'dashed')+
        stat_summary(aes(group = Protected), fun = mean, geom = "line", size = 2)
    print(plot1)
    
}
   

#remove the day 0 samples
Luminex8 <- Luminex7[,-(1:22)]
colnames(Luminex7)
colnames(Luminex8)

#Merge
All4 <- merge(All3, Luminex8, by = 'VolID', all.x = T)
nrow(All4) == 19 #check all matched
All3$VolID[!All3$VolID %in% Luminex7$VolID]





######################################################################
#Load in IgG1 FC
######################################################################
IgG1 <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 9)


table(IgG1$Volunteer)
colnames(IgG1) <- c('Timepoint', 'VolID', 'IgG1_FC') #rename

#make to numeric
IgG1$IgG1_FC[IgG1$IgG1_FC == '-99'] <- NA
IgG1$IgG1_FC <- gsub('>', '', IgG1$IgG1_FC)
IgG1$IgG1_FC <- as.numeric(IgG1$IgG1_FC)
IgG1$IgG1_FC <- log2(IgG1$IgG1_FC)
hist(IgG1$IgG1_FC)

#exploratory plots
IgG1_clin <- merge(All4, IgG1, by = 'VolID')

ggplot(IgG1_clin, aes(x=Timepoint, y=IgG1_FC, group = VolID, 
                                   colour = Protected)) +  geom_line()  +
        facet_grid(.~Group)+  
        geom_vline(xintercept = 5, linetype = 'dashed')+
        stat_summary(aes(group = Protected), fun = mean, geom = "line", size = 2)


IgG1_2 <- dcast(IgG1, VolID ~ Timepoint)
colnames(IgG1_2)  <- gsub('-', '', colnames(IgG1_2))   
colnames(IgG1_2)  <- gsub('00', '0', colnames(IgG1_2))  
IgG1_3 <- IgG1_2[,-2]
colnames(IgG1_3)[2:ncol(IgG1_3)] <- paste0('IgG1_', colnames(IgG1_3)[2:ncol(IgG1_3)])


#Merge
All5 <- merge(All4, IgG1_3, by = 'VolID', all.x = T)
nrow(All5) == 19 #check all matched
All4$VolID[!All4$VolID %in% IgG1_3$VolID]




######################################################################
#Load in IgG4 FC
######################################################################
IgG4 <- read.xlsx('Raw/Integrated data file_05082021.xlsx', 11)


table(IgG4$Volunteer)
colnames(IgG4) <- c('Timepoint', 'VolID', 'IgG4_FC') #rename

#make to numeric
IgG4$IgG4_FC[IgG4$IgG4_FC == '-99'] <- NA
IgG4$IgG4_FC <- log2(IgG4$IgG4_FC)
hist(IgG4$IgG4_FC)

#exploratory plots
IgG4_clin <- merge(All4, IgG4, by = 'VolID')

ggplot(IgG4_clin, aes(x=Timepoint, y=IgG4_FC, group = VolID, 
                      colour = Protected)) +  geom_line()  +
    facet_grid(.~Group)+  
    geom_vline(xintercept = 5, linetype = 'dashed')+
    stat_summary(aes(group = Protected), fun = mean, geom = "line", size = 2)


IgG4_2 <- dcast(IgG4, VolID ~ Timepoint)
colnames(IgG4_2)  <- gsub('-', '', colnames(IgG4_2))   
colnames(IgG4_2)  <- gsub('00', '0', colnames(IgG4_2))  
IgG4_3 <- IgG4_2[,-2]
colnames(IgG4_3)[2:ncol(IgG4_3)] <- paste0('IgG4_', colnames(IgG4_3)[2:ncol(IgG4_3)])


#Merge
All6 <- merge(All5, IgG4_3, by = 'VolID', all.x = T)
nrow(All6) == 19 #check all matched
All5$VolID[!All5$VolID %in% IgG4_3$VolID]


#write table
write.xlsx(All6, 'All_ITTCHI_LumFC.xlsx', row.names=F)
write.table(All6, 'All_ITTCHI_LumFC.txt', sep='\t', quote=F, row.names=F)


print(Sys.time())
print(sessionInfo())
