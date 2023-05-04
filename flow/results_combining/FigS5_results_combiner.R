library(tidyverse)
library(limma)

figure = "FigS5"

setwd("/location_of_data/results_tables/")

rep1 = read.table(paste0(figure,"_results_Rep1.txt"),header=T)
rep2 = read.table(paste0(figure,"_results_Rep2.txt"),header=T)
rep3 = read.table(paste0(figure,"_results_Rep3.txt"),header=T)

rep1$Replicate = "Rep1"
rep2$Replicate = "Rep2"
rep3$Replicate = "Rep3"


all = rbind(rep1,rep2)
all = rbind(all,rep3)
all = separate(all,exp_name,into=c("Enzyme","Activator","Target"),remove=F)
all = mutate(all,"PercPos_PhiYFP" = round(100*num_phiyfp/num_singlets,1),
             "PercPos_mCerulean" = round(100*num_mcerulean/num_singlets,1),
             "PercPos_mCherry" = round(100*num_mcherry/num_singlets,1))
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$Enzyme = "None"
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$Activator = "None"
all[grep("UT|PhiYFP|mCerulean|mCherry",all$exp_name),]$Target = "None"
all[grep("P1-only",all$exp_name),]$Target = "P1"
all[grep("P3-only",all$exp_name),]$Target = "P3"
write.table(all,paste0(figure,"_combined_table.txt"),row.names=F,quote=F,sep="\t")

all$Channel = NA
all$Channel[grep("P1|PhiYFP",all$exp_name)] = "PhiYFP"
all$Channel[grep("P3|mCherry",all$exp_name)] = "mCherry"
uts = dplyr::filter(all,exp_name == "UT")
all = rbind(all,uts)
all[all$exp_name == "UT" & all$Replicate == "Rep1",]$Channel = c("PhiYFP","mCherry")
all[all$exp_name == "UT" & all$Replicate == "Rep2",]$Channel = c("PhiYFP","mCherry")
all[all$exp_name == "UT" & all$Replicate == "Rep3",]$Channel = c("PhiYFP","mCherry")


all$PercPos = all$PercPos_PhiYFP
all$PercPos[grep("mCerulean",all$Channel)] = all$PercPos_mCerulean[grep("mCerulean",all$Channel)]
all$PercPos[grep("mCherry",all$Channel)] = all$PercPos_mCherry[grep("mCherry",all$Channel)]

all$MFI = all$phiyfp_MFI
all$MFI[grep("mCerulean",all$Channel)] = all$mcerulean_MFI[grep("mCerulean",all$Channel)]
all$MFI[grep("mCherry",all$Channel)] = all$mcherry_MFI[grep("mCherry",all$Channel)]
all_onechannel = select(all,exp_name:Target,Replicate,Channel:MFI)
colnames(all_onechannel)[grep("exp_name",colnames(all_onechannel))] = "Name"

write.table(all_onechannel,paste0(figure,"_combined_onechannel_table.txt"),row.names=F,quote=F,sep="\t")
