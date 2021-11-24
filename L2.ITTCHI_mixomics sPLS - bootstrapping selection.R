rm(list=ls()) #clear workspace

#load in packages
library(plyr)
library(dplyr)
library(corrplot)
library(pROC)
library(tidyverse)
library(tidyr)
library(mixOmics)
library(NMF)
library(RColorBrewer)
library(cowplot)
library(scales)
library(knitr)
library(UpSetR)
library(grid)
library(igraph)
library(reshape2)
library(pheatmap)
library(corrplot)

#setwd
setwd('R:/Para-CIH/CIH Group member folders/Simon/Results/ITTCHI mixomics/')

## load data
data <- read.delim('All_ITTCHI_LumFC.txt')
#data_raw <- read.delim('All_ITTCHI.txt')

write_dir <- './Results/L2.sPLS_LOO_bootstrap/'
name <- 'sPLSC_LuminexFC_boot_'
  
  
dir.create(write_dir)
setwd(write_dir)

table(colSums(is.na(data))<5) #check for missing data 
table(apply(data, 2, sd, na.rm=T) == 0) #check data without any changes


rownames(data) <- data$VolID


#remove the one without protection data
data <- data[!is.na(data$Protected),]


colnames(data)

#mixomics data format per datatype

AE <- data[,3:12]
AE$Sex[AE$Sex == 'Male'] <- 0
AE$Sex[AE$Sex == 'Female'] <- 1
AE$Sex <- as.numeric(AE$Sex)
AE$Skin_AEs[AE$Skin_AEs == 'No'] <- 0
AE$Skin_AEs[AE$Skin_AEs == 'Yes'] <- 1
AE$Skin_AEs <- as.numeric(AE$Skin_AEs)
AE$Severe_abdominal_AE[AE$Severe_abdominal_AE == 'No'] <- 0
AE$Severe_abdominal_AE[AE$Severe_abdominal_AE == 'Yes'] <- 1
AE$Severe_abdominal_AE <- as.numeric(AE$Severe_abdominal_AE)

Eosino <- data[,18:43]
Luminex <- data[,44:263] #get Luminex from the raw data so not normalized
Antibodies <- data[,264:283]


#X <- list(Eosino = Eosino, AE = AE, Antibodies = Antibodies, Luminex = Luminex)
X <- cbind(Eosino, AE, Antibodies, Luminex)
#X <- Luminex
Y <- data$Mean_egg_load

X <- X[rownames(X) != '875',]
X <- X[,colSums(apply(X, 2, is.na))==0]
X <- X[,which(!colnames(X) == 'IFN_g_B24')]#remove the IFNg_b24 that is only in the outlier different
X.scale <- data.frame(scale(X)) #if you scale manually it will come with 12 features instead of 26

Y <- Y[rownames(Eosino) != '875']

data2 <- data[data$VolID != '875',]
data2 <- data2[,which(!colnames(data2) == 'IFN_g_B24')]#remove the IFNg_b24 that is only in the outlier different


tuning <- tune.spls(X, Y, validation = 'loo', scale = T,
                    ncomp = 1, test.keepX = 1:50, measure =  'MAE' )
pdf(paste0(name, '_tuning.pdf'), width =3, height = 3)  
plot(tuning)
dev.off()

tuning$choice.keepX

sample <- data2$VolID



important <- list()
for(i in 1:100){

  rm <- sample(1:17, 15, replace = F) #SAMPLE 14 donors at random

    trainX <- X[rm,] 
    trainY <- Y[rm]
    
    sp <- spls(trainX, trainY, scale = T, mode = 'regression', 
               keepX = tuning$choice.keepX, keepY = 1, ncomp = 1)

    
    important[[i]] <- names(sp$loadings$X[,1][abs(sp$loadings$X[,1])>0])
    
}

  important2 <- as.character(do.call(c, important))
  write.table(important2, paste0('features_in_all_bootstraps.txt'), sep = '\t', row.names = F)
  imp <- data.frame(table(important2)/100*100)
  imp2 <- imp[imp$Freq>10,]
  imp2$important2 <- factor(imp2$important2, 
                            levels = imp2$important2[order(imp2$Freq)])
  ggplot(imp2, aes(x=important2, y=Freq)) + 
    geom_bar(stat='identity') + 
    xlab('Variable') +
    ylab('Importance') + 
    coord_flip()  + theme_bw()
  ggsave(paste0('keep_', keep, 'importance.pdf'), width=8, height=12)
  write.table(imp2, paste0('keep_', keep, 'freq_feats_in_all_LOO.txt'), sep = '\t', row.names = F)
  


#build the full model
sp <- spls(X, Y, mode = 'regression', keepX = tuning$choice.keepX, keepY = 1, ncomp = 1)


pdf(paste0(name, '_plotLoadings.pdf'), width =6, height = 6)  
plotLoadings(sp, block = 'X')
dev.off()

pos = data.frame(Name = rownames(loadings(sp)$X)[loadings(sp)$X>0],
                 Loading = loadings(sp)$X[loadings(sp)$X>0])
neg = data.frame(Name = rownames(loadings(sp)$X)[loadings(sp)$X<0],
                 Loading = loadings(sp)$X[loadings(sp)$X<0])

all <- rbind(pos, neg)
all$Name <- factor(all$Name, 
                          levels = all$Name[order(abs(all$Loading), decreasing = T)])

all$Consensus <- 'No'
cons <- as.character(imp2$important2[imp2$Freq>50])
all$Consensus[all$Name %in% cons] <- 'Yes' 
all$Direction <- 'Protect'
all$Direction[all$Loading > 0 ] <- 'Suscept'  
  
ggplot(all, aes(x=Loading, y = Name, fill = Direction, alpha=Consensus)) + 
  geom_bar(stat='identity')  + theme_bw() + 
  scale_alpha_manual(values = c(0.2, 1)) + 
  scale_fill_manual(values = c('red', 'blue'))
ggsave('features_w_consensus.pdf', width=12, height=12, units = 'cm')



keep <- as.character(all$Name[all$Consensus == 'Yes'])
data_keep <- data2[,c('VolID', 'Protected', 'Group', 'Mean_egg_load', keep)]  
write.table(all, 'loadings_X.txt', sep = '\t', row.names = F)
data_keep$Skin_AEs[data_keep$Skin_AEs == 'No'] <- 0
data_keep$Skin_AEs[data_keep$Skin_AEs == 'Yes'] <- 1
data_keep$Skin_AEs <- as.numeric(data_keep$Skin_AEs)

data.melt <- melt(data_keep, id = c(colnames(data_keep)[1:4], 'Skin_AEs'))
data.melt <- melt(data_keep, id = c(colnames(data_keep)[1:4]))
ggplot(data.melt, aes(x=Mean_egg_load, y=value)) + 
  facet_wrap(~variable, scale = 'free') + 
  geom_point(aes(colour = Group))  + 
  stat_smooth(method = 'lm', colour = 'black', size = 2) + 
  theme_bw() + theme(aspect.ratio = 1) +
  scale_colour_manual(values = c('red', 'blue'))
ggsave('Features_vs_load.pdf', width = 10, height = 10)


ann <- data.frame(group <- data$Group,
                  Load <- data$Mean_egg_load)
rownames(ann) <- rownames(data)

pheat <- scale(data_keep[,5:ncol(data_keep)])
pheat[pheat>2] <- 2
pheat[pheat < -2] <- -2

pdf('heatmap_selected.pdf', width = 10, height = 7)
pheatmap(t(pheat), annotation_col = ann)
dev.off()

pheat2 <- pheat[order(data_keep$Mean_egg_load, decreasing = T),]
pdf('heatmap_selected_arranged.pdf', width = 10, height = 7)
pheatmap(pheat2, annotation_row = ann, cluster_rows = F)
dev.off()


corMat <- cor(data_keep[,4:ncol(data_keep)], use = 'na.or.complete')  
corMat[lower.tri(corMat)] <- 0
diag(corMat) <- 0
color.blocks = brewer.pal(n = 12, name = "Paired")[seq(2, 10, by = 2)]
links <- corMat %>% 
  as.data.frame() %>% 
  mutate(to = rownames(.)) %>% 
  gather(from, cor, -to) %>% 
  filter(abs(cor) > 0.5) %>% 
  mutate(Color = ifelse(cor > 0, "red", "blue"))
nodes <- data.frame(id = unique(c(links$to, links$from)))
nodes$datasets <- 'empty' 
nodes$datasets[grepl('_B', as.character(nodes$id))] <- 'Challenge' 
nodes$datasets[!grepl('_B', as.character(nodes$id))] <- 'Immunization' 
nodes$datasets[grepl('egg', as.character(nodes$id))] <- 'Outcome' 

net <- graph_from_data_frame(d=links, vertices=nodes, directed=FALSE) 
E(net)$color <- links$Color
V(net)$color <- color.blocks[as.numeric(factor(nodes$datasets))]
E(net)$weight <- abs(links$cor)
weight <- E(net)$weight
weight[weight == 0.5] <- 0.3

pdf(paste0(name, '_network manual.pdf'), width = 8, height = 8)
plot(net, edge.curved=.2, vertex.label.cex=1, vertex.size=6, 
      vertex.label.color="black")
dev.off()



corrvl <- Hmisc::rcorr(as.matrix(data_keep[,4:ncol(data_keep)]), type = 'pearson')
colnames(corrvl$r)
corr2 <- corrvl$r
ann <- data.frame(type = rep('Clinical', ncol(corr2)))
row.names(ann) = colnames(corr2)
ann$type[grepl("Nose", row.names(ann))] <- 'Nose'
ann$type[grepl("CyTOF", row.names(ann))] <- 'WB_CyTOF'
ann$type[grepl("IM", row.names(ann))] <- 'WB_IM'

corrplot(corr2)

stats <- corrvl$P
stats[corrvl$P > 0.05] <- ''
stats[corrvl$P < 0.05] <- '*'
stats[corrvl$P < 0.01] <- '**'
stats[corrvl$P < 0.001] <- '***'
diag(stats) <- ''


dev.off()
pdf('heatmap.pdf', width = 10, height = 10)
pheatmap(corr2, border_color = NA, display_numbers = stats)
dev.off()


save.image(paste0(name, 'spls.RData'))

