library("vegan")

library("dplyr")

library("doParallel")

library("foreach")

library("mgcv")

library("reshape2")

library("ggplot2")

library("cowplot")

library("Rcpp")

library("RcppArmadillo")


rm(list = ls())

gc()

metadata_file = "all_meta.csv" 

count_matrix = "all_otu.csv" 

EM_iterations = 1000 #default value=1000

different_sources_flag = 0

print("Change directory path")

dir_path = paste("C:/ ...your path/ FEAST/") 
setwd(paste0(dir_path, "FEAST_src"))

source("src.R")

setwd(paste0(dir_path, "Data_files"))

metadata <- read.csv(metadata_file, header=T, sep = ",", row.names = 1)

otus <- read.csv(count_matrix, header=T, sep = ",", row.names = 1)

otus <- t(as.matrix(otus))

common.sample.ids <- intersect(rownames(metadata), rownames(otus))

otus <- otus[common.sample.ids,]

metadata <- metadata[common.sample.ids,]

if(length(common.sample.ids) <= 1) {

 message <- paste(sprintf('Error: there are %d sample ids in common '),

                  'between the metadata file and data table')

 stop(message)

}

if(different_sources_flag == 0){

 metadata$id[metadata$SourceSink == 'Source'] = NA

 metadata$id[metadata$SourceSink == 'Sink'] = c(1:length(which(metadata$SourceSink == 'Sink')))

}

envs <- metadata$Env

Ids <- na.omit(unique(metadata$id))

Proportions_est <- list()

for(it in 1:length(Ids)){

   train.ix <- which(metadata$SourceSink=='Source' & metadata$id == Ids[it])

   test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])

 }

 else{

   train.ix <- which(metadata$SourceSink=='Source')

   test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])

 }
 
 num_sources <- length(train.ix)

 COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user

 str(COVERAGE)

 sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))

 sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))

 print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))

 print(paste("Seq depth in the sources and sink samples = ",COVERAGE))

 print(paste("The sink is:", envs[test.ix]))

 FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)

 Proportions_est[[it]] <- FEAST_output$data_prop[,1]

 names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")

 if(length(Proportions_est[[it]]) < num_sources +1){

   tmp = Proportions_est[[it]]

   Proportions_est[[it]][num_sources] = NA

   Proportions_est[[it]][num_sources+1] = tmp[num_sources]

 }

 print("Source mixing proportions")

 print(Proportions_est[[it]])

}

print(Proportions_est)


FEAST_output = as.data.frame(Proportions_est)

colnames(FEAST_output) = paste("repeat_",Ids,sep = "") 

head(FEAST_output)

filename = paste(dir_path,"Result/FEAST.csv",sep = "")

write.csv(FEAST_output,filename,quote = F)
