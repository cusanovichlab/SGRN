library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
setwd("/location_of_data/")

#Figure 1
table_fig1 = read.table("./results_tables/Fig1_combined_onechannel_table.txt",header=T)

p1_filtered = filter(table_fig1,Target == "P1" & Enzyme != "None")
p1_filtered$Amount = as.numeric(gsub("ng","",p1_filtered$Amount))
p1_dCas9_mean = p1_filtered |>
  filter(Enzyme == "dCas9") |>
  summarize(mean(PercPos))
p1_fitter = filter(p1_filtered,Enzyme != "dCas9")
p1_fit = lm(PercPos ~ Amount,p1_fitter)
p1_xval = approx(x = p1_fit$fitted.values, y = p1_fitter$Amount,xout=as.numeric(p1_dCas9_mean[1,1]))$y
p1_xval
#[1] 697.2112

ggplot(p1_fitter,aes(x=Amount, y=PercPos)) +
  geom_quasirandom(alpha=1, color="dodgerblue2",size=2,width=20) +
  geom_smooth(method=lm,se=FALSE,color="black",size=0.5) +
  geom_hline(yintercept=as.numeric(p1_dCas9_mean[1,1]),linetype = "dashed") +
  theme_bw()

p2_filtered = filter(table_fig1,Target == "P2" & Enzyme != "None")
p2_filtered$Amount = as.numeric(gsub("ng","",p2_filtered$Amount))
p2_dCas9_mean = p2_filtered |>
  filter(Enzyme == "dCas9") |>
  summarize(mean(PercPos))
p2_fitter = filter(p2_filtered,Enzyme != "dCas9")
p2_fit = lm(PercPos ~ Amount,p2_fitter)
p2_xval = approx(x = p2_fit$fitted.values, y = p2_fitter$Amount,xout=as.numeric(p2_dCas9_mean[1,1]))$y
p2_xval
#[1] 498.7975

ggplot(p2_fitter,aes(x=Amount, y=PercPos)) +
  geom_quasirandom(alpha=1, color="dodgerblue2",size=2,width=20) +
  geom_smooth(method=lm,se=FALSE,color="black",size=0.5) +
  geom_hline(yintercept=as.numeric(p2_dCas9_mean[1,1]),linetype = "dashed") +
  theme_bw()

p3_filtered = filter(table_fig1,Target == "P3" & Enzyme != "None")
p3_filtered$Amount = as.numeric(gsub("ng","",p3_filtered$Amount))
p3_dCas9_mean = p3_filtered |>
  filter(Enzyme == "dCas9") |>
  summarize(mean(PercPos))
p3_fitter = filter(p3_filtered,Enzyme != "dCas9")
p3_fit = lm(PercPos ~ Amount,p3_fitter)
p3_xval = approx(x = p3_fit$fitted.values, y = p3_fitter$Amount,xout=as.numeric(p3_dCas9_mean[1,1]))$y
p3_xval
#[1] 695.5582

ggplot(p3_fitter,aes(x=Amount, y=PercPos)) +
  geom_quasirandom(alpha=1, color="dodgerblue2",size=2,width=20) +
  geom_smooth(method=lm,se=FALSE,color="black",size=0.5) +
  geom_hline(yintercept=as.numeric(p3_dCas9_mean[1,1]),linetype = "dashed") +
  theme_bw()

p1_order = c("UT","PhiYFP","P1-only","dCas9-P1-250ng","dCas12a-P1-250ng","dCas12a-P1-500ng","dCas12a-P1-750ng")
p1_pos_plot = table_fig1 |> filter(Channel == "PhiYFP") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p1_order) +
  coord_cartesian(ylim=c(0, 70))
p1_mfi_plot = table_fig1 |> filter(Channel == "PhiYFP") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p1_order)

p2_order = c("UT","mCerulean","P2-only","dCas9-P2-250ng","dCas12a-P2-250ng","dCas12a-P2-500ng","dCas12a-P2-750ng")
p2_pos_plot = table_fig1 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_order) +
  coord_cartesian(ylim=c(0, 70))
p2_mfi_plot = table_fig1 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_order)

p3_order = c("UT","mCherry","P3-only","dCas9-P3-250ng","dCas12a-P3-250ng","dCas12a-P3-500ng","dCas12a-P3-750ng")
p3_pos_plot = table_fig1 |> filter(Channel == "mCherry") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_order) +
  coord_cartesian(ylim=c(0, 70))
p3_mfi_plot = table_fig1 |> filter(Channel == "mCherry") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_order)

pdf("./Fig1_PerPos_plots.pdf")
p1_pos_plot
p2_pos_plot
p3_pos_plot
dev.off()
pdf("./Fig1_MFI_plots.pdf")
p1_mfi_plot
p2_mfi_plot
p3_mfi_plot
dev.off()

###############################################
#Figure S3
###############################################
table_figS3 = read.table("./results_tables/FigS3_combined_onechannel_table.txt",header=T)

table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA2" & Promoter == "P1") |>
  summarize(mean(PercPos))
#1     0.4333333
table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA3" & Promoter == "P1") |>
  summarize(mean(PercPos))
#1     0.6666667
table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA1" & Promoter == "P2") |>
  summarize(mean(PercPos))
#1           2.1
table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA3" & Promoter == "P2") |>
  summarize(mean(PercPos))
#1      2.033333
table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA1" & Promoter == "P3") |>
  summarize(mean(PercPos))
#1      5.733333
table_figS3 |> filter(Enzyme == "dCas12a" & Guide == "crRNA2" & Promoter == "P3") |>
  summarize(mean(PercPos))
#1           5.7

table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA2" & Promoter == "P1") |>
  summarize(mean(PercPos))
#1     0.3333333
table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA3" & Promoter == "P1") |>
  summarize(mean(PercPos))
#1     0.3666667
table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA1" & Promoter == "P2") |>
  summarize(mean(PercPos))
#1      2.366667
table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA3" & Promoter == "P2") |>
  summarize(mean(PercPos))
#1      2.233333
table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA1" & Promoter == "P3") |>
  summarize(mean(PercPos))
#      5.166667
table_figS3 |> filter(Enzyme == "dCas9" & Guide == "gRNA2" & Promoter == "P3") |>
  summarize(mean(PercPos))
#1      6.233333

p1_order = c("UT","PhiYFP","P1-only","dCas9-gRNA1-P1","dCas9-gRNA2-P1","dCas9-gRNA3-P1","dCas12a-crRNA1-P1","dCas12a-crRNA2-P1","dCas12a-crRNA3-P1")
p1_xr_pos_plot = table_figS3 |> filter(Channel == "PhiYFP") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p1_order) +
  coord_cartesian(ylim=c(0, 70))
p1_xr_mfi_plot = table_figS3 |> filter(Channel == "PhiYFP") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p1_order)

p2_order = c("UT","mCerulean","P2-only","dCas9-gRNA1-P2","dCas9-gRNA2-P2","dCas9-gRNA3-P2","dCas12a-crRNA1-P2","dCas12a-crRNA2-P2","dCas12a-crRNA3-P2")
p2_xr_pos_plot = table_figS3 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_order) +
  coord_cartesian(ylim=c(0, 70))
p2_xr_mfi_plot = table_figS3 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_order)

p3_order = c("UT","mCherry","P3-only","dCas9-gRNA1-P3","dCas9-gRNA2-P3","dCas9-gRNA3-P3","dCas12a-crRNA1-P3","dCas12a-crRNA2-P3","dCas12a-crRNA3-P3")
p3_xr_pos_plot = table_figS3 |> filter(Channel == "mCherry") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_order) +
  coord_cartesian(ylim=c(0, 70))
p3_xr_mfi_plot = table_figS3 |> filter(Channel == "mCherry") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_order)

pdf("./FigS3_PercPos_plots.pdf")
p1_xr_pos_plot
p2_xr_pos_plot
p3_xr_pos_plot
dev.off()
pdf("./FigS3_MFI_plots.pdf")
p1_xr_mfi_plot
p2_xr_mfi_plot
p3_xr_mfi_plot
dev.off()

#####################################
#Figure 2
#####################################
table_fig2 = read.table("./results_tables/Fig2_combined_onechannel_table.txt",header=T)

table_fig2_dCas12a = filter(table_fig2,Enzyme == "dCas12a")
table_fig2_dCas9 = filter(table_fig2,Enzyme == "dCas9")

table_fig2_dCas12a$Length = as.numeric(gsub("nt","",table_fig2_dCas12a$Length))
table_fig2_dCas9$Length = as.numeric(gsub("nt","",table_fig2_dCas9$Length))

table_fig2_pos_cas12_fit = lm(PercPos*100 ~ as.numeric(Length),table_fig2_dCas12a)
summary(table_fig2_pos_cas12_fit)$coefficient[2,4]
#[1] 0.0001133547

table_fig2_pos_cas9_fit = lm(PercPos*100 ~ as.numeric(Length),table_fig2_dCas9)
summary(table_fig2_pos_cas9_fit)$coefficient[2,4]
#[1] 0.3140895

table_fig2_mfi_cas12_fit = lm(log10(MFI) ~ as.numeric(Length),table_fig2_dCas12a)
summary(table_fig2_mfi_cas12_fit)$coefficient[2,4]
#[1] 0.0101847

table_fig2_mfi_cas9_fit = lm(log10(MFI) ~ as.numeric(Length),table_fig2_dCas9)
summary(table_fig2_mfi_cas9_fit)$coefficient[2,4]
#[1] 7.533063e-05

p2_length_order = c("UT","mCerulean","P2-only","dCas9-16nt-P2","dCas9-18nt-P2","dCas9-20nt-P2","dCas12a-16nt-P2","dCas12a-18nt-P2","dCas12a-20nt-P2")
p2_length_pos_plot = table_fig2 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_length_order) +
  coord_cartesian(ylim=c(0, 70))
p2_length_mfi_plot = table_fig2 |> filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p2_length_order)

pdf("./Fig2_PercPos_plots.pdf")
p2_length_pos_plot
dev.off()
pdf("./Fig2_MFI_plots.pdf")
p2_length_mfi_plot
dev.off()

#########################################
#Figure 3
#########################################
table_fig3 = read.table("./results_tables/Fig3_combined_onechannel_table.txt",header=T)
table_fig3$TargetN = as.numeric(gsub("xTS","",table_fig3$TargetN))

table_fig3_cas12 = table_fig3 |> filter(Enzyme == "dCas12a")
table_fig3_mfi_cas12_fit = lm(log10(MFI) ~ as.numeric(TargetN),table_fig3_cas12)
summary(table_fig3_mfi_cas12_fit)$coefficient[2,4]
#[1] 4.117159e-05

table_fig3_cas9 = table_fig3 |> filter(Enzyme == "dCas9")
table_fig3_mfi_cas9_fit = lm(log10(MFI) ~ as.numeric(TargetN),table_fig3_cas9)
summary(table_fig3_mfi_cas9_fit)$coefficient[2,4]
#[1] 6.793302e-09

table_fig3_1x_cas9 = filter(table_fig3,Enzyme != "dCas12a" & TargetN == 1)
table_fig3_pos_1x_cas9_fit = lm(PercPos ~ Enzyme,table_fig3_1x_cas9)
summary(table_fig3_pos_1x_cas9_fit)$coefficient[2,4]
#[1] 0.003495145

table_fig3_1x_cas12 = filter(table_fig3,Enzyme != "dCas9" & TargetN == 1)
table_fig3_pos_1x_cas12_fit = lm(PercPos ~ Enzyme,table_fig3_1x_cas12)
summary(table_fig3_pos_1x_cas12_fit)$coefficient[2,4]
#[1] 0.000436783

p3_ts_order = c("UT","mCerulean","P2-1xTS-only","P2-2xTS-only","P2-3xTS-only",
                "P2-4xTS-only","P2-5xTS-only","P2-6xTS-only","P2-7xTS-only",
                "P2-Original-only","dCas9-P2-1xTS","dCas9-P2-2xTS","dCas9-P2-3xTS",
                "dCas9-P2-4xTS","dCas9-P2-5xTS","dCas9-P2-6xTS","dCas9-P2-7xTS",
                "dCas9-P2-Original","dCas12a-P2-1xTS","dCas12a-P2-2xTS",
                "dCas12a-P2-3xTS","dCas12a-P2-4xTS","dCas12a-P2-5xTS",
                "dCas12a-P2-6xTS","dCas12a-P2-7xTS","dCas12a-P2-Original")
p3_ts_pos_plot = table_fig3 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_ts_order) +
  coord_cartesian(ylim=c(0, 70))
p3_ts_mfi_plot = table_fig3 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p3_ts_order)

pdf("./Fig3_PercPos_plots.pdf")
p3_ts_pos_plot
dev.off()
pdf("./Fig3_MFI_plots.pdf")
p3_ts_mfi_plot
dev.off()

##########################################
#Figure 4
##########################################
table_fig4 = read.table("./results_tables/Fig4_combined_onechannel_table.txt",header=T)
table_fig4_dCas12a = dplyr::filter(table_fig4,Enzyme == "dCas12a")

table_fig4_mfi_fit = lm(log10(MFI) ~ Promoter + Activator,table_fig4_dCas12a)
summary(table_fig4_mfi_fit)$coefficient[2,4]
#[1] 1.238138e-07
summary(table_fig4_mfi_fit)$coefficient[3,4]
#[1] 4.71082e-05

table_fig4_pos_fit = lm(PercPos ~ Promoter + Activator,table_fig4_dCas12a)
summary(table_fig4_pos_fit)$coefficient[2,4]
#[1] 2.21417e-05
summary(table_fig4_pos_fit)$coefficient[3,4]
#[1] 0.02687992

#Reduction in MFI for UbCp
10^(summary(table_fig4_mfi_fit)$coefficient[1,1] + summary(table_fig4_mfi_fit)$coefficient[2,1]) - 10^(summary(table_fig4_mfi_fit)$coefficient[1,1])
#[1] -1757.96

#Increase in MFI for VP64
10^(summary(table_fig4_mfi_fit)$coefficient[1,1] + summary(table_fig4_mfi_fit)$coefficient[3,1]) - 10^(summary(table_fig4_mfi_fit)$coefficient[1,1])
#[1] 1806.608

p4_domains_order = c("UT","mCerulean","P2-only","dCas9-EF1a-VP64-P2","dCas12a-EF1a-Act-P2","dCas12a-EF1a-VP64-P2","dCas12a-UbCp-Act-P2","dCas12a-UbCp-VP64-P2")
p4_domains_pos_plot = table_fig4 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p4_domains_order) +
  coord_cartesian(ylim=c(0, 75))
p4_domains_mfi_plot = table_fig4 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p4_domains_order)

pdf("./Fig4_PercPos_plots.pdf")
p4_domains_pos_plot
dev.off()
pdf("./Fig4_MFI_plots.pdf")
p4_domains_mfi_plot
dev.off()

#################################
#Figure S5
#################################
table_figS5 = read.table("./results_tables/FigS5_combined_onechannel_table.txt",header=T)
table_figS5_dCas12a_p1 = dplyr::filter(table_figS5,Enzyme == "dCas12a" & Target == "P1")
table_figS5_dCas12a_p3 = dplyr::filter(table_figS5,Enzyme == "dCas12a" & Target == "P3")

table_figS5_p1_mfi_fit = lm(log10(MFI) ~ Activator,table_figS5_dCas12a_p1)
summary(table_figS5_p1_mfi_fit)$coefficient[2,4]
#[1] 0.0004176716

table_figS5_p1_pos_fit = lm(PercPos ~ Activator,table_figS5_dCas12a_p1)
summary(table_figS5_p1_pos_fit)$coefficient[2,4]
#[1] 0.3404665

table_figS5_p3_mfi_fit = lm(log10(MFI) ~ Activator,table_figS5_dCas12a_p3)
summary(table_figS5_p3_mfi_fit)$coefficient[2,4]
#[1] 0.925289

table_figS5_p3_pos_fit = lm(PercPos ~ Activator,table_figS5_dCas12a_p3)
summary(table_figS5_p3_pos_fit)$coefficient[2,4]
#[1] 0.3374567

#Reduction in MFI for VP64 for P1
10^(summary(table_figS5_p1_mfi_fit)$coefficient[1,1] + summary(table_figS5_p1_mfi_fit)$coefficient[2,1]) - 10^(summary(table_figS5_p1_mfi_fit)$coefficient[1,1])
#[1] -552.9182

p5_domains_p1_order = c("UT","PhiYFP","P1-only","dCas9-VP64-P1","dCas12a-Act-P1","dCas12a-VP64-P1")
p5_domains_p1_pos_plot = table_figS5 |> dplyr::filter(Channel == "PhiYFP") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_domains_p1_order) +
  coord_cartesian(ylim=c(0, 70))
p5_domains_p1_mfi_plot = table_figS5 |> dplyr::filter(Channel == "PhiYFP") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_domains_p1_order)

p5_domains_p3_order = c("UT","mCherry","P3-only","dCas9-VP64-P3","dCas12a-Act-P3","dCas12a-VP64-P3")
p5_domains_p3_pos_plot = table_figS5 |> dplyr::filter(Channel == "mCherry") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_domains_p3_order) +
  coord_cartesian(ylim=c(0, 70))
p5_domains_p3_mfi_plot = table_figS5 |> dplyr::filter(Channel == "mCherry") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_domains_p3_order)


pdf("./FigS5_PercPos_plots.pdf")
p5_domains_p1_pos_plot
p5_domains_p3_pos_plot
dev.off()
pdf("./FigS5_MFI_plots.pdf")
p5_domains_p1_mfi_plot
p5_domains_p3_mfi_plot
dev.off()

##########################################
#Figure S6
##########################################
table_figS6 = read.table("./results_tables/FigS6_combined_onechannel_table.txt",header=T)
table_figS6_1s = filter(table_figS6,TargetN == "1xTS" & Enzyme == "dCas12a")
summary(lm(log10(MFI) ~ Activator, data = table_figS6_1s))$coefficient[2,4]
#[1] 0.01664105
summary(lm(PercPos ~ Activator, data = table_figS6_1s))$coefficient[2,4]
#[1] 0.5087808

table_figS6_7s = filter(table_figS6,TargetN == "7xTS" & Enzyme == "dCas12a")
summary(lm(log10(MFI) ~ Activator, data = table_figS6_7s))$coefficient[2,4]
#[1] 0.01323353
summary(lm(PercPos ~ Activator, data = table_figS6_7s))$coefficient[2,4]
#[1] 0.7700199

ps6_steric_order = c("UT","mCerulean","P2-1xTS-only","P2-7xTS-only",
                     "dCas12a-Act-1xTS-P2","dCas12a-Act-7xTS-P2",
                     "dCas12a-VP64-1xTS-P2","dCas12a-VP64-7xTS-P2")
ps6_steric_pos_plot = table_figS6 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=ps6_steric_order) +
  coord_cartesian(ylim=c(0, 75))
ps6_steric_mfi_plot = table_figS6 |> dplyr::filter(Channel == "mCerulean") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 11000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=ps6_steric_order)

pdf("./FigS6_PercPos_plots.pdf")
ps6_steric_pos_plot
dev.off()
pdf("./FigS6_MFI_plots.pdf")
ps6_steric_mfi_plot
dev.off()

####################################
#Figure 5 NLS
####################################
table_fig5NLS = read.table("./results_tables/Fig5NLS_combined_onechannel_table.txt",header=T)
table_fig5NLS_dCas12a = filter(table_fig5NLS,Enzyme == "dCas12a" & Amount == "500ng")
table_fig5NLS_mfi_fit = lm(log10(MFI) ~ Activator + xNLS,table_fig5NLS_dCas12a)
#Activator p-value
summary(table_fig5NLS_mfi_fit)$coefficient[2,4]
#[1] 0.004287748

#xNLS p-value
summary(table_fig5NLS_mfi_fit)$coefficient[3,4]
#[1] 2.244476e-10

table_fig5NLS_pos_fit = lm(PercPos ~ Activator + xNLS,table_fig5NLS_dCas12a)
#Activator p-value
summary(table_fig5NLS_pos_fit)$coefficient[2,4]
#[1] 0.8329247

#xNLS p-value
summary(table_fig5NLS_pos_fit)$coefficient[3,4]
#[1] 2.437625e-07

#Increase in MFI with 3xNLS
10^(summary(table_fig5NLS_mfi_fit)$coefficient[1,1]+summary(table_fig5NLS_mfi_fit)$coefficient[3,1]) - 10^summary(table_fig5NLS_mfi_fit)$coefficient[1,1]
#[1] 6304.532

#Decrease in MFI with VP64
10^(summary(table_fig5NLS_mfi_fit)$coefficient[1,1] + summary(table_fig5NLS_mfi_fit)$coefficient[2,1]) - 10^summary(table_fig5NLS_mfi_fit)$coefficient[1,1]
#[1] -215.9345

p5_nls_order = c("UT","mCerulean","P2-only","dCas9-3xNLS-VP64-250ng-P2",
                 "dCas12a-1xNLS-Activ-500ng-P2","dCas12a-1xNLS-VP64-500ng-P2",
                 "dCas12a-3xNLS-Activ-500ng-P2","dCas12a-3xNLS-VP64-500ng-P2")
p5_nls_pos_plot = table_fig5NLS |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_nls_order) +
  coord_cartesian(ylim=c(0, 70))
p5_nls_mfi_plot = table_fig5NLS |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_nls_order)

pdf("./Fig5_NLS_PercPos_plots.pdf")
p5_nls_pos_plot
dev.off()
pdf("./Fig5_NLS_MFI_plots.pdf")
p5_nls_mfi_plot
dev.off()

##############################################
#Figure S8
##############################################
#(uses same data data table as Fig5NLS)
table_figS8_3x_dCas12a = filter(table_fig5NLS,Enzyme == "dCas12a" & xNLS == "3xNLS")
table_figS8_3x_dCas12a$Amount = gsub("ng","",table_figS8_3x_dCas12a$Amount)

table_figS8_titer_mfi_fit = lm(log10(MFI) ~ Activator + as.numeric(Amount),table_figS8_3x_dCas12a)
summary(table_figS8_titer_mfi_fit)$coefficient[2,4]
#[1] 5.221905e-05

table_figS8_titer_pos_fit = lm(PercPos ~ Activator + as.numeric(Amount),table_figS8_3x_dCas12a)
summary(table_figS8_titer_pos_fit)$coefficient[3,4]
#[1] 1.135245e-05

p5_nls_titer_order = c("dCas9-3xNLS-VP64-250ng-P2","dCas12a-1xNLS-Activ-250ng-P2",
                 "dCas12a-1xNLS-Activ-500ng-P2","dCas12a-1xNLS-Activ-750ng-P2",
                 "dCas12a-3xNLS-Activ-250ng-P2","dCas12a-3xNLS-Activ-500ng-P2",
                 "dCas12a-3xNLS-Activ-750ng-P2","dCas12a-3xNLS-VP64-250ng-P2",
                 "dCas12a-3xNLS-VP64-500ng-P2","dCas12a-3xNLS-VP64-750ng-P2")
p5_nls_titer_pos_plot = table_fig5NLS |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_nls_titer_order) +
  coord_cartesian(ylim=c(0, 70))
p5_nls_titer_mfi_plot = table_fig5NLS |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_nls_titer_order)

pdf("./FigS8_titer_PercPos_plots.pdf")
p5_nls_titer_pos_plot
dev.off()
pdf("./FigS8_titer_MFI_plots.pdf")
p5_nls_titer_mfi_plot
dev.off()

####################################
#Figure 5 Cascade
####################################
table_fig5Cascade = read.table("./results_tables/Fig5Cascade_combined_onechannel_table.txt",header=T)
table_fig5Cascade_functional = filter(table_fig5Cascade,Guide != "None" & Target1 != "None" & Target2 != "None")

table_fig5Cascade_cas = filter(table_fig5Cascade_functional,Enzyme == "dCas9" | Activator == "Activ" & Enzyme == "dCas12a")
table_fig5Cascade_activator = filter(table_fig5Cascade_functional,Enzyme == "dCas12a")

table_fig5Cascade_mfi_cas_fit = lm(log10(MFI) ~ Enzyme,table_fig5Cascade_cas)
summary(table_fig5Cascade_mfi_cas_fit)$coefficient[2,4]
#[1] 0.0004780915

table_fig5Cascade_pos_cas_fit = lm(PercPos ~ Enzyme,table_fig5Cascade_cas)
summary(table_fig5Cascade_pos_cas_fit)$coefficient[2,4]
#[1] 0.01013185

table_fig5Cascade_mfi_activator_fit = lm(log10(MFI) ~ Activator,table_fig5Cascade_activator)
summary(table_fig5Cascade_mfi_activator_fit)$coefficient[2,4]
#[1] 0.0006085531

table_fig5Cascade_pos_activator_fit = lm(PercPos ~ Activator,table_fig5Cascade_activator)
summary(table_fig5Cascade_pos_activator_fit)$coefficient[2,4]
#[1] 0.02113967

p5_cascade_order = c("UT","mCerulean","P2-only","P1-28nt-P2-only",
                     "dCas9-VP64-gRNA1-None-P2","dCas9-VP64-gRNA1-P1-P2",
                     "P1-DR-P2-only","dCas12a-VP64-crRNA1-None-P2",
                     "dCas12a-VP64-crRNA1-P1-P2","dCas12a-Activ-crRNA1-None-P2",
                     "dCas12a-Activ-crRNA1-P1-P2")
p5_cascade_pos_plot = table_fig5Cascade |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x=Name,y=PercPos)) +
  geom_quasirandom(aes(colour=Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x =element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_cascade_order) +
  coord_cartesian(ylim=c(0, 70))
p5_cascade_mfi_plot = table_fig5Cascade |> dplyr::filter(Channel == "PacificOrange") |>
  ggplot(aes(x = Name,y = MFI)) +
  scale_y_log10() + coord_cartesian(ylim = c(100, 10000)) +
  geom_quasirandom(aes(colour = Name),show.legend = FALSE) +
  theme_bw() + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_x_discrete(limits=p5_cascade_order)

pdf("./Fig5Cascade_PercPos_plots.pdf")
p5_cascade_pos_plot
dev.off()
pdf("./Fig5Cascade_MFI_plots.pdf")
p5_cascade_mfi_plot
dev.off()