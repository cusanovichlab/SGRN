
########################################################################################
########## Community Code For Nuclear Localization of GFP-Tagged dCas Enzymes ########## 
################################## by Reeba Varghese ###################################
################################ University of Arizona #################################

### Load the required packages and library to read the metadata using R

shhh <-suppressWarnings(suppressMessages(suppressPackageStartupMessages))
shhh(require(data.table)) #This loads the fread function
shhh(library(dplyr)) #This uses the tbl_df function
shhh(library(tidyr))

### Specify the directory and folder that contains the raw data. 
#All subsequent files created from this analysis will be saved in the specified folder. 

InDir <- "/location_of_data/"

### Specify and read the raw data file that will be analyzed. 
#This code reads .csv file formats. 

DataFile <- "Kang_et_al_2023_Nuclear_GFP_localization_raw_data"
    DataFile2 <- fread(paste0(InDir,DataFile,".csv"))
    Data <-tbl_df(DataFile2)
    Data <- Data[,colSums(is.na(Data))<nrow(Data)]# Remove all NA columns
    head(Data)

### Read the method file that will add metadata information to the raw data. 

MethodFile <- paste0(InDir,"Kang_et_al_2023_Nuclear_GFP_localization_MethodFile.csv")

#Read_metadata function allows you to read the metadata from the directory.
Read_metadata2 <- function(methodFile){
  Methods=read.table(methodFile,stringsAsFactors = FALSE,header=F,sep=',')
  #This will specifically allow you to read it in 96 well layout format
  Condition <- Methods[2:9,2:13]
  row.names(Condition) <- Methods[2:9,1]
  row <- Methods[2:9,1]
  col <- c(1:12)
  
#This code will remove NA or empty wells
  logic=is.na(Condition) | Condition==""
  indeces=which(!logic,arr.ind=TRUE, useNames = TRUE)
  Wells=paste0(rownames(indeces),indeces[,2],sep='')
  
#This code will add different Condition[indeces] information
  for (i in 1:(nrow(Methods)/9)){
    if (i==1){
      mydata <- cbind(Wells,Condition[indeces])
      nm=c(Methods[1,1])
    }   else{
      temp=Methods[seq(9*(i-1)+2,9*i,1),2:13]
      row.names(temp) <- Methods[2:9,1]
      nm <- c(nm,Methods[9*(i-1)+1,1])
      mydata <- cbind(mydata,temp[indeces])   
    }
  }
    Info <- data.frame(matrix(mydata,length( Condition[indeces]),(nrow(Methods)/9)+1,dimnames = list(c(),c("Wells",nm))))
  Plates_info <-data.frame(Methods[-c(1),14:ncol(Methods)])
  Plates_info[,1] <- as.integer(Plates_info[,1])
  names(Plates_info) <- Methods[1,14:ncol(Methods)]
  indeces2=which(!(is.na(Plates_info )|Plates_info==""),arr.ind=TRUE, useNames = TRUE)
  
  Plane <- Plates_info [unique(indeces2[,1]),]
    
  return(list("Info"=Info,"Plane"=Plane))
}

### Read the meta data

Meta_data <- Read_metadata2(MethodFile)
WellInfo <- data.frame(Meta_data["Info"])
PlateInfo <- data.frame(Meta_data["Plane"])
names(WellInfo) = gsub(pattern = "Info.", replacement = "", x = names(WellInfo))
names(PlateInfo) = gsub(pattern = "Plane.", replacement = "", x = names(PlateInfo))
WellInfo

### Now combine the raw data with the meta data. 
#Any additional calculations that need to be done can be done here.

FinalData<- Data[complete.cases(Data),]%>%
  separate(XYSiteName,c("Wells","Field"),sep="_")%>%
  mutate(Row=substr(Wells,1,1),Column=as.numeric(substr(Wells,2,3)),Wells=paste0(Row,Column),Field=as.numeric(Field))%>% 

  left_join(WellInfo,by="Wells")%>% 

    mutate(Transfection=ifelse(CellSegSum488Int>500000,"GFP-positive","GFP-negative"),               
          CytoSegSum488Int= (`CellSegSum488Int`-`NucSegSum488Int`), # This will calculate the sum pixel intensity of the cytoplasmic annulus ring around the nucleus              
          CytoSegArea= (`CellSegArea`-`NucSegArea`) # This will calculate the area of the cytoplasmic annulus ring around the nucleus
#          CytoplasmicIntensitytoArea=(`CytoSegSum488Int`/`CytoSegArea`)
          )
head (FinalData)

### Format the data to prepare it for analysis. 

FinalData2 <-FinalData%>%
    select(-MultiPointIndex,-Entity,-Row,-Column,-CellSegArea,-CellSegSum488Int) #this will remove any columns not needed

###If you want to reorder the columns to a specific order: 

col_order <- c("ObjectId", "Field", "Wells","Vector","Vector_Rep","NucSegArea","NucSegSum488Int","CytoSegArea","CytoSegSum488Int","Transfection")
FinalData3 <- FinalData2[, col_order]
FinalData3

###Create a new .CSV file with the raw data and the meta data combined

write.csv(file=paste0(InDir,"./data/Kang_et_al_2023_Nuclear_GFP_localization_final_data.csv"),FinalData3,row.names=F)


