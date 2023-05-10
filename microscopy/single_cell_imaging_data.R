library(tidyverse)
setwd("/path_to_script/")

images = read.csv("./data/Kang_et_al_2023_Nuclear_GFP_localization_final_data.csv",header=T)
images = images |>
  rename(Well = Wells) |>
  mutate(NuclearCytoplasmicIntensityRatio =
           round((NucSegSum488Int/NucSegArea)/(CytoSegSum488Int/CytoSegArea),4)) |>
  select(-Transfection)
#write("./data/Table_S9.txt",row.names=F,col.names=T,quote=F,sep="\t")

images |>
  filter(NucSegSum488Int + CytoSegSum488Int > 500000) |>
  ggplot(aes(x=factor(Vector,levels=c("dCas9-VP64","dCas12a-[Activ]","dCas12a-3xNLS-[Activ]","dCas12a-3xNLS-VP64")),y=log2(NuclearCytoplasmicIntensityRatio),fill=Vector_Rep)) +
  geom_violin() + scale_fill_manual(values=rep(c("indianred","dodgerblue2","orange","mediumorchid"),each=3)) +
  theme_bw()
ggsave("Single_cell_image_nuclear_ratios.pdf",height=2.5,width=7.5)

ratios = images |>
  filter(NucSegSum488Int + CytoSegSum488Int > 500000) |>
  group_by(Vector_Rep) |>
  summarize(ave = mean(log2(NuclearCytoplasmicIntensityRatio)), med = median(log2(NuclearCytoplasmicIntensityRatio))) |>
  separate(col=Vector_Rep,into=c("vector","rep"),sep="_")

orig_vs_3xActiv = filter(ratios,vector=="dCas12a-[Activ]" | vector=="dCas12a-3xNLS-[Activ]")
orig_vs_3xActiv_fit = lm(ave ~ vector,orig_vs_3xActiv)
summary(orig_vs_3xActiv_fit)$coefficients[2,4]
#[1] 8.336357e-06

orig_vs_3xVP64 = filter(ratios,vector=="dCas12a-[Activ]" | vector=="dCas12a-3xNLS-VP64")
orig_vs_3xVP64_fit = lm(ave ~ vector,orig_vs_3xVP64)
summary(orig_vs_3xVP64_fit)$coefficients[2,4]
#[1] 9.74146e-05

orig_vs_dCas9 = filter(ratios,vector=="dCas12a-[Activ]" | vector=="dCas9-VP64")
orig_vs_dCas9_fit = lm(ave ~ vector,orig_vs_dCas9)
summary(orig_vs_dCas9_fit)$coefficients[2,4]
#[1] 7.46487e-07
