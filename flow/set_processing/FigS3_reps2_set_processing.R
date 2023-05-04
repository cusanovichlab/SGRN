library(flowCore)
library(flowViz)
library(flowUtils)
library(ggcyto)
library(ggplot2)
library(viridis)
library(limma)
library(scales)

figure = "FigS3"
replicate = "Rep2"
exp_dir = paste0("/location_of_data/")
setwd(exp_dir)
create.dir("./results_plots/")
create.dir("./results_tables/")
create.dir("./gating_files/")
fullset = read.flowSet(path="./fcs_files/",pattern=paste0(figure,"*",replicate,"*.fcs$"))

trLogicle=logicletGml2(parameters = "FITC-A", T = 262144, M = 4.5, 
                           W = 0.5, A = 0, transformationId="trLogicle")
trLogicle2=logicletGml2(parameters = "AmCyan-A", T = 262144, M = 4.5, 
                           W = 0.5, A = 0, transformationId="trLogicle2")
trLogicle3=logicletGml2(parameters = "PE-Cy5-A", T = 262144, M = 4.5, 
                           W = 0.5, A = 0, transformationId="trLogicle3")

logicler = function(x){
  logicle_FITC_A  = eval(trLogicle)(exprs(x)[,4])
  logicle_AmCyan_A  = eval(trLogicle2)(exprs(x)[,5])
  logicle_PE_Cy5_A  = eval(trLogicle3)(exprs(x)[,6])
  new_logicles = cbind(logicle_FITC_A,logicle_AmCyan_A)
  new_logicles = cbind(new_logicles,logicle_PE_Cy5_A)
  colnames(new_logicles) = c("logicle-FITC-A","logicle-AmCyan-A","logicle-PE-Cy5-A")
  x = fr_append_cols(x, new_logicles)
}
fullset = fsApply(fullset,logicler)
shortnames = strsplit2(pData(fullset)$name,"_")[,2]
fullset@phenoData@data$name = shortnames

cells_coords = matrix(c(50000,40000,70000,184000,246000,222000,182000,50000,10000,16000,85000,212000),nrow=6,ncol=2)
colnames(cells_coords) = c("FSC-A","SSC-A")
cells_gate = polygonGate(filterId = "cells_gate", .gate = cells_coords)
p = ggcyto(fullset, aes(x = `FSC-A`, y =  `SSC-A`))
p = p + geom_hex(bins = 256) + theme_minimal() + scale_fill_viridis(trans="log10")
p = p + geom_gate(cells_gate) + geom_stats(size = 3)

cells = Subset(fullset,cells_gate)
singlets_coords = matrix(c(45000,100000,220000,250000,225000,35000,30000,10000,50000,160000,210000,215000,45000,25000),nrow=7,ncol=2)
colnames(singlets_coords) = c("FSC-A","FSC-H")
singlets_gate = polygonGate(filterId = "singlets_gate", .gate = singlets_coords)
p1 = ggcyto(cells, aes(x = `FSC-A`, y =  `FSC-H`)) + scale_y_continuous(labels = scientific) + scale_x_continuous(labels = scientific)
p1 = p1 + geom_hex(bins = 256) + theme_minimal() + scale_fill_viridis(trans="log10")
p1 = p1 + geom_gate(singlets_gate) + geom_stats(size = 3)

singlets = Subset(cells,singlets_gate)
phiyfp_gate = rectangleGate(filterId = "phiyfp",list("logicle-FITC-A"=c(0.2,1),"SSC-A"=c(5000,215000)))
p2 = ggcyto(singlets, aes(x = `logicle-FITC-A`, y =  `SSC-A`))
p2 = p2 + geom_hex(bins = 256) + theme_minimal() + scale_fill_viridis(trans="log10")
p2 = p2 + geom_gate(phiyfp_gate) + geom_stats(size = 3) + scale_y_continuous(labels = scientific)
p5 = ggcyto(singlets, aes(x = `logicle-FITC-A`, fill = factor(name))) + geom_density(alpha = 0.2) + theme_minimal()

mcerulean_coords = matrix(c(0.2,0.26,0.28,0.38,1,1,5000,40000,62000,215000,215000,5000),nrow=6,ncol=2)
colnames(mcerulean_coords) = c("logicle-AmCyan-A","SSC-A")
mcerulean_gate = polygonGate(filterId = "mcerulean", .gate = mcerulean_coords)
p3 = ggcyto(singlets, aes(x = `logicle-AmCyan-A`, y =  `SSC-A`)) + scale_y_continuous(labels = scientific)
p3 = p3 + geom_hex(bins = 256) + theme_minimal() + scale_fill_viridis(trans="log10")
p3 = p3 + geom_gate(mcerulean_gate) + geom_stats(size = 3)
p6 = ggcyto(singlets, aes(x = `logicle-AmCyan-A`, fill = factor(name))) + geom_density(alpha = 0.2) + theme_minimal()

mcherry_gate = rectangleGate(filterId = "mcherry",list("logicle-PE-Cy5-A"=c(0.22,1),"SSC-A"=c(5000,215000)))
p4 = ggcyto(singlets, aes(x = `logicle-PE-Cy5-A`, y =  `SSC-A`)) + scale_y_continuous(labels = scientific)
p4 = p4 + geom_hex(bins = 256) + theme_minimal() + scale_fill_viridis(trans="log10")
p4 = p4 + geom_gate(mcherry_gate) + geom_stats(size = 3)
p7 = ggcyto(singlets, aes(x = `logicle-PE-Cy5-A`, fill = factor(name))) + geom_density(alpha = 0.2) + theme_minimal()

pdf(paste0("./results_plots/",figure,"_gating_results_",replicate,".pdf"),height = 10,width=17)
p
p1
p2
p3
p4
p5
p6
p7
dev.off()

phiyfp_positives = Subset(singlets,phiyfp_gate)
mcerulean_positives = Subset(singlets,mcerulean_gate)
mcherry_positives = Subset(singlets,mcherry_gate)

num_events = fsApply(fullset,function(x) nrow(exprs(x)))
num_cells = fsApply(cells,function(x) nrow(exprs(x)))
num_singlets = fsApply(singlets,function(x) nrow(exprs(x)))
num_events = fsApply(fullset,function(x) nrow(exprs(x)))
num_cells = fsApply(cells,function(x) nrow(exprs(x)))
num_singlets = fsApply(singlets,function(x) nrow(exprs(x)))
num_phiyfp = fsApply(phiyfp_positives,function(x) nrow(exprs(x)))
num_mcerulean = fsApply(mcerulean_positives,function(x) nrow(exprs(x)))
num_mcherry = fsApply(mcherry_positives,function(x) nrow(exprs(x)))
phiyfp_MFI = fsApply(phiyfp_positives,function(x) median(exprs(x)[,4]))
mcerulean_MFI = fsApply(mcerulean_positives,function(x) median(exprs(x)[,5]))
mcherry_MFI = fsApply(mcherry_positives,function(x) median(exprs(x)[,6]))

results = data.frame("exp_name" = shortnames,"exp_date" = fullset[[1]]@description$`$DATE`,"num_events" = num_events,
                     "num_cells" = num_cells,"num_singlets" = num_singlets,"num_phiyfp" = num_phiyfp,
                     "num_mcerulean" = num_mcerulean, "num_mcherry" = num_mcherry,"phiyfp_MFI" = round(phiyfp_MFI,1),
                     "mcerulean_MFI" = round(mcerulean_MFI,1), "mcherry_MFI" = round(mcherry_MFI,1))
write.table(results,paste0("./results_tables/",figure,"_results_",replicate,".txt"),quote=F,sep="\t",row.names=F)

flowEnv=new.env()
flowEnv[['trLogicle']]=trLogicle
flowEnv[['trLogicle2']]=trLogicle2
flowEnv[['trLogicle3']]=trLogicle3
flowEnv[['cells_gate']] = cells_gate
flowEnv[['singlets_gate']] = singlets_gate
good_cells_gate = new("intersectFilter", filterId="good_cells_gate", filters=list(cells_gate, singlets_gate))
flowEnv[['good_cells_gate']] = good_cells_gate
tr1Pars=list(transformReference("trLogicle",flowEnv))
tr2Pars=list(transformReference("trLogicle2",flowEnv))
tr3Pars=list(transformReference("trLogicle3",flowEnv))
phiyfp_gate@parameters=new("parameters",tr1Pars)
flowEnv[['phiyfp_gate']]=phiyfp_gate
mcerulean_gate@parameters=new("parameters",tr2Pars)
flowEnv[['mcerulean_gate']]=mcerulean_gate
mcherry_gate@parameters=new("parameters",tr3Pars)
flowEnv[['mcherry_gate']]=mcherry_gate
write.gatingML(flowEnv,paste0("./",subdir,"renamed/",figure,"_gates_",replicate,".gating-ml2.xml"))
