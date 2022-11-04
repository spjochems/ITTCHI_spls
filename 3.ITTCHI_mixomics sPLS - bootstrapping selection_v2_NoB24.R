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
setwd('R:/Para-CIH/CIH Group member folders/Simon/Results/ITTCHI mixomics/Final/')

## load data
data <- read.delim('All_ITTCHI_LumFC.txt')
#data_raw <- read.delim('All_ITTCHI.txt')

write_dir <- '../Results/L4.sPLS_LOO_bootstrap_noB24/'
name <- 'sPLSC_LuminexFC_boot_'
  
  
dir.create(write_dir)
setwd(write_dir)

table(colSums(is.na(data))<5) #check for missing data 
table(apply(data, 2, var, na.rm=T) == 0) #check data without any changes

#filter data with too low variance
data0 <- data
hist(log10(apply(data0, 2, var, na.rm=T)))
data <- data0[,unique(c(1:6,17,which(apply(data0, 2, var, na.rm=T) > 0.1)))]


plot(x= log10(apply(data, 2, var, na.rm=T)),
     y=log10(apply(data, 2, sd, na.rm=T)  ))
rownames(data) <- data$VolID


#remove the one without protection data
data <- data[!is.na(data$Protected),]


colnames(data)

#mixomics data format per datatype

AE <- data[,c(3:6, 8:13)]
AE$Sex[AE$Sex == 'Male'] <- 0
AE$Sex[AE$Sex == 'Female'] <- 1
AE$Sex <- as.numeric(AE$Sex)
AE$Skin_AEs[AE$Skin_AEs == 'No'] <- 0
AE$Skin_AEs[AE$Skin_AEs == 'Yes'] <- 1
AE$Skin_AEs <- as.numeric(AE$Skin_AEs)
AE$Severe_abdominal_AE[AE$Severe_abdominal_AE == 'No'] <- 0
AE$Severe_abdominal_AE[AE$Severe_abdominal_AE == 'Yes'] <- 1
AE$Severe_abdominal_AE <- as.numeric(AE$Severe_abdominal_AE)

Eosino <- data[,18:31]
Luminex <- data[,32:185] 
Antibodies <- data[,186:205]


#X <- list(Eosino = Eosino, AE = AE, Antibodies = Antibodies, Luminex = Luminex)
X <- cbind(Eosino, AE, Antibodies, Luminex)


#only keep immunization timepoints
X <- X[,!(grepl('_B24', colnames(X)) | grepl('_B17', colnames(X)) )]



#X <- Luminex
Y <- data$Mean_egg_load

X <- X[rownames(X) != '875',] #exclude the donor without complete datasets
X <- X[,colSums(apply(X, 2, is.na))==0]
X.scale <- data.frame(scale(X)) #scale manually
Y <- Y[rownames(Eosino) != '875']

data2 <- data[data$VolID != '875',]


#tuning <- tune.spls(X, Y, validation = 'loo', scale = T,
 #                   ncomp = 1, test.keepX = 1:50, measure =  'MAE' )

tuning <- tune.spls(X.scale, Y, validation = 'loo', scale = F,
                    ncomp = 1, test.keepX = 1:50, measure =  'MAE' )

pdf(paste0(name, '_tuning.pdf'), width =4, height = 4)  
plot(tuning)
dev.off()

tuning$choice.keepX

sample <- data2$VolID



important <- list()
for(i in 1:100){

  rm <- sample(1:17, 15, replace = F) #SAMPLE 15 donors at random

    trainX <- X.scale[rm,] 
    trainY <- Y[rm]
    
    sp <- spls(trainX, trainY, scale = F, mode = 'regression', 
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
sp <- spls(X.scale, Y, mode = 'regression', keepX = tuning$choice.keepX, keepY = 1, ncomp = 1)


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
cons <- as.character(imp2$important2[imp2$Freq>25])
all$Consensus[all$Name %in% cons] <- 'Yes' 
all$Direction <- 'Protect'
all$Direction[all$Loading > 0 ] <- 'Suscept'  

all2 <- all

all$Name <- gsub('IFN_g', 'IFNg', all$Name)
all$Name <- gsub('IL_13', 'IL13', all$Name)
all$Name <- gsub('MIP_1b', 'MIP1b', all$Name)
all$Name <- gsub('PDGF_bb', 'PDGFbb', all$Name)
all$Name <- gsub('MCP_1', 'MCP1', all$Name)
all$Name <- gsub('IL_9', 'IL9', all$Name)
all$Name <- gsub('_A', '_I', all$Name)
all$Name <- gsub('_B', '_C', all$Name)
all$Name <- gsub('TNF_a', 'TNFa', all$Name)
all$Name <- gsub('Skin_IEs', 'Skin_AEs', all$Name)

all$Name <- factor(all$Name, 
                   levels = all$Name[order(abs(all$Loading), decreasing = T)])

ggplot(all, aes(x=Loading, y = Name, fill = Direction)) + 
  geom_bar(stat='identity')  + theme_bw() + 
  scale_alpha_manual(values = c(0.2, 1)) + 
  scale_fill_manual(values = c('red', 'blue'))
ggsave('features_w_consensus2.pdf', width=12, height=12, units = 'cm')


write.table(all, 'loadings_X.txt', sep = '\t', row.names = F)


keep <- as.character(all2$Name[all2$Consensus == 'Yes'])
data_keep <- data2[,c('VolID', 'Protected', 'Group', 'Mean_egg_load', keep)]  

data_keep$Skin_AEs[data_keep$Skin_AEs == 'No'] <- 0
data_keep$Skin_AEs[data_keep$Skin_AEs == 'Yes'] <- 1
data_keep$Skin_AEs <- as.numeric(data_keep$Skin_AEs)

colnames(data_keep) <- gsub('IFN_g', 'IFNg', colnames(data_keep))
colnames(data_keep) <- gsub('IL_13', 'IL13', colnames(data_keep))
colnames(data_keep) <- gsub('MIP_1b', 'MIP1b', colnames(data_keep))
colnames(data_keep) <- gsub('PDGF_bb', 'PDGFbb', colnames(data_keep))
colnames(data_keep) <- gsub('MCP_1', 'MCP1', colnames(data_keep))
colnames(data_keep) <- gsub('IL_9', 'IL9', colnames(data_keep))
colnames(data_keep) <- gsub('_A', '_I', colnames(data_keep))
colnames(data_keep) <- gsub('_B', '_C', colnames(data_keep))
colnames(data_keep) <- gsub('Skin_IEs', 'Skin_AEs', colnames(data_keep))

data.melt <- melt(data_keep, id = c(colnames(data_keep)[1:4], 'Skin_AEs'))
data.melt <- melt(data_keep, id = c(colnames(data_keep)[1:4]))
ggplot(data.melt, aes(x=Mean_egg_load, y=value)) + 
  facet_wrap(~variable, scale = 'free', ncol = 4) + 
  geom_point(aes(colour = Group, shape=Group))  + 
  stat_smooth(method = 'lm', colour = 'black', size = 2) + 
  theme_bw() + theme(aspect.ratio = 1) +
  scale_colour_manual(values = c('#008000', '#0000C0')) + 
  scale_shape_manual(values = c(15,16))
ggsave('Features_vs_load2.pdf', width = 18, height = 20, units= 'cm')


ann <- data.frame(group <- data2$Group,
                  Load <- data2$Mean_egg_load)
rownames(ann) <- rownames(data2)

pheat <- scale(data_keep[,5:ncol(data_keep)])
pheat[pheat>2] <- 2
pheat[pheat < -2] <- -2

mat <- pheat[order(ann$Load....data2.Mean_egg_load),]
mat1 <- mat[rownames(mat) %in% rownames(ann)[ann$group....data2.Group == 'Intervention'],]
mat2 <- mat[rownames(mat) %in% rownames(ann)[ann$group....data2.Group != 'Intervention'],]

mat <- pheat[order(ann$group....data2.Group),]
colnames(mat)
pdf('heatmap_selected_manualordered.pdf', width = 10, height = 7)
pheatmap(t(rbind(mat1, mat2)), annotation_col = ann, cluster_cols = F)
dev.off()


pdf('heatmap_selected2.pdf', width = 10, height = 7)
pheatmap(t(pheat), annotation_col = ann)
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
nodes$datasets[grepl('_C', as.character(nodes$id))] <- 'Challenge' 
nodes$datasets[!grepl('_I', as.character(nodes$id))] <- 'Immunization' 
nodes$datasets[grepl('egg', as.character(nodes$id))] <- 'Outcome' 

net <- graph_from_data_frame(d=links, vertices=nodes, directed=FALSE) 
E(net)$color <- links$Color
V(net)$color <- color.blocks[as.numeric(factor(nodes$datasets))]
E(net)$weight <- abs(links$cor)
weight <- E(net)$weight
weight[weight == 0.5] <- 0.3

pdf(paste0(name, '_network manual2.pdf'), width = 8, height = 8)
plot(net, edge.curved=.2, vertex.label.cex=1, vertex.size=6, 
      vertex.label.color="black")
dev.off()



corrvl <- Hmisc::rcorr(as.matrix(data_keep[,4:ncol(data_keep)]), type = 'pearson')
colnames(corrvl$r)
corr2 <- corrvl$r

stats <- corrvl$P
stats[corrvl$P > 0.05] <- ''
stats[corrvl$P < 0.05] <- '*'
stats[corrvl$P < 0.01] <- '**'
stats[corrvl$P < 0.001] <- '***'
diag(stats) <- ''


dev.off()
pdf('heatmap_correlations.pdf', width = 9, height = 9)
pheatmap(corr2, border_color = NA)
dev.off()


data_eo <- read.delim('../../All_ITTCHI_LumFC.txt')
data_eo2 <- data_eo[!is.na(data_eo$Protected),1:43]
data_eo3 <- data_eo2[data_eo2$VolID != '875',]
data_eo4 <- data_eo3[,colSums(apply(data_eo3, 2, is.na)) == 0]

colnames(data_eo4)
times <- gsub('Eosinophils_A', '' , colnames(data_eo4))
times <- gsub('Eosinophils_B', '' , times)
times <- as.numeric(times)
times[23:37] <- times[23:37] + 13

rowSums(data_eo4[,17:36]/29)

data_eo4$AUCall <- 0
for(i in 1:17){
  auc <- pkr::AUC(times[17:36], as.numeric(data_eo4[i,17:36]))
  data_eo4$AUCall[i] <- auc[20,1]
}

data_eo4$AUCB <- 0
for(i in 1:17){
  auc <- pkr::AUC(times[23:36], as.numeric(data_eo4[i,23:36]))
  data_eo4$AUCB[i] <- auc[14,1]
}

data_eo4$AUC1216 <- 0
for(i in 1:17){
  auc <- pkr::AUC(times[32:36], as.numeric(data_eo4[i,32:36]))
  data_eo4$AUC1216[i] <- auc[5,1]
}



ggplot(data_eo4, aes(x=AUCall, y = Mean_egg_load, colour = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1))
cor.test(data_eo4$AUCall, data_eo4$Mean_egg_load)

ggplot(data_eo4, aes(x=AUCB, y = Mean_egg_load, colour = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1))
cor.test(data_eo4$AUCB, data_eo4$Mean_egg_load)

ggplot(data_eo4, aes(x=AUC1216, y = Mean_egg_load, colour = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1))
cor.test(data_eo4$AUC1216, data_eo4$Mean_egg_load)



ggplot(data_eo4, aes(x=AUC1216, y = Mean_egg_load, colour = Skin_AEs, shape = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1), colour= 'black') + 
  theme_bw() + theme(aspect.ratio = 1) +
  scale_colour_manual(values = c('gold', 'red')) + 
  scale_shape_manual(values = c(15,16))
ggsave('AUC_eosino_skin_vs_load2.pdf', width = 8, height = 8, units= 'cm')




ggplot(data_eo4, aes(x=AUCall, y = Mean_egg_load, colour = Skin_AEs)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1)) +
  facet_wrap(.~Group)

ggplot(data_eo4, aes(x=AUCall, y = Mean_egg_load, colour = Skin_AEs, shape = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1)) 

ggplot(data_eo4, aes(x=Group, y = AUCall, colour = Group)) + 
  geom_point() + stat_smooth(method = 'lm', aes(group = 1))

model <- lm(Mean_egg_load ~ as.factor(Skin_AEs) + AUC1216 , data_eo4)
summary(model)

cor.test(data_eo4$Mean_egg_load, data_eo4$AUC1216)
save.image(paste0(name, 'spls2.RData'))

