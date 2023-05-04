library(tidyverse)
library(limma)

figure = "Fig5NLS"

setwd("/location_of_data/results_tables/")

rep1 = read.table(paste0(figure,"_results_Rep1.txt"),header=T)
rep2 = read.table(paste0(figure,"_results_Rep2.txt"),header=T)
rep3 = read.table(paste0(figure,"_results_Rep3.txt"),header=T)

rep1$Replicate = "Rep1"
rep2$Replicate = "Rep2"
rep3$Replicate = "Rep3"


all = rbind(rep1,rep2)
all = rbind(all,rep3)
all = separate(all,exp_name,into=c("Enzyme","xNLS","Activator","Amount","Target"),remove=F)
all = mutate(all,"PercPos_PhiYFP" = round(100*num_phiyfp/num_singlets,1),
             "PercPos_mCerulean" = round(100*num_mcerulean/num_singlets,1),
             "PercPos_mCherry" = round(100*num_mcherry/num_singlets,1),
             "PercPos_PacificOrange" = round(100*num_pacorange/num_singlets,1))
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$Enzyme = "None"
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$Activator = "None"
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$xNLS = "0xNLS"
all[grep("UT|only|PhiYFP|mCerulean|mCherry",all$exp_name),]$Amount = "0ng"
all[grep("UT|PhiYFP|mCerulean|mCherry",all$exp_name),]$Target = "None"
all[grep("P2-only",all$exp_name),]$Target = "P2"

write.table(all,paste0(figure,"_combined_table.txt"),row.names=F,quote=F,sep="\t")
all$Channel = "PacificOrange"

all$PercPos = all$PercPos_PacificOrange

all$MFI = all$pacorange_MFI
all_onechannel = select(all,exp_name:Amount,Replicate,Channel:MFI)
colnames(all_onechannel)[grep("exp_name",colnames(all_onechannel))] = "Name"

write.table(all_onechannel,paste0(figure,"_combined_onechannel_table.txt"),row.names=F,quote=F,sep="\t")
