library(ggplot2)

setwd("/Users/natorro/Desktop/curso_r_ifisica/example_code/chara_yucateca-genomicPolimor/CY-Denovo-Ref-Cb-Zf")
archivos <- list.files(pattern = "formato_sumstats*")
for (file in archivos){
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_sumstats_(.*)\\..*", "\\1", file)
    temp_dataset$Experiment <- filename
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_sumstats_(.*)\\..*", "\\1", file)
    dataset$Experiment <- filename
  }
}

denovo_ref <- ggplot(dataset, aes(factor(Experiment), Fis))
denovo_ref + geom_boxplot() + stat_boxplot ()
denovo_ref <- ggplot(dataset, aes(factor(Experiment), Pi))
denovo_ref + geom_boxplot () + stat_boxplot ()

ggplot(dataset) + geom_point(aes(x=ExpHet, y=ObsHet, colour=Experiment), alpha=I(1/10)) + facet_grid(. ~ Experiment)

maf <- ggplot(dataset, aes(x=P)) + facet_grid(. ~ Experiment)
maf + geom_histogram(aes(y = ..density..))
rm(list=ls())
###############
library(ggplot2)
setwd("~/CY-Denovo-Ref-Cb-Zf")
archivos <- list.files(pattern = "formato_hapstats*")
for (file in archivos){
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_hapstats_(.*)\\..*", "\\1", file)
    temp_dataset$Experiment <- filename
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_hapstats_(.*)\\..*", "\\1", file)
    dataset$Experiment <- filename
  }
}

denovo_ref <- ggplot(dataset, aes(factor(Experiment), GeneDiversity))
denovo_ref + geom_boxplot () + stat_boxplot ()
denovo_ref <- ggplot(dataset, aes(factor(Experiment), HaplotypeDiversity))
denovo_ref + geom_boxplot () + stat_boxplot ()
rm(list=ls())
######################
library(ggplot2)
setwd("~/CY-Denovo-Ref-Cb-Zf/ref_compara")
archivos <- list.files(pattern = "formato_hapstats*")
for (file in archivos){
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_hapstats_(.*)\\..*", "\\1", file)
    temp_dataset$Experiment <- filename
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_hapstats_(.*)\\..*", "\\1", file)
    dataset$Experiment <- filename
  }
}

ggplot(dataset) + geom_point(aes(x=SmoothedHaplotypeDiversity, y=SmoothedHaplotypeDiversityP.value, colour=Experiment), alpha=I(1/20)) + facet_grid(. ~ Experiment)
ggplot(dataset) + geom_point(aes(x=SmoothedGeneDiversity, y=SmoothedGeneDiversityP.value, colour=Experiment), alpha=I(1/20)) + facet_grid(. ~ Experiment)
ggplot(dataset) + geom_point(aes(x=BP, y=SmoothedHaplotypeDiversity, colour=Experiment), alpha=I(1/20)) + facet_grid(. ~ Experiment)
ggplot(dataset) + geom_point(aes(x=BP, y=SmoothedGeneDiversity, colour=Experiment), alpha=I(1/20)) + facet_grid(. ~ Experiment)
ggplot(dataset) + geom_point(aes(x=BP, y=N, colour=Experiment), alpha=I(1/20)) + facet_grid(. ~ Experiment)
rm(list=ls())
####################
library(ggplot2)
setwd("~/CY-Denovo-Ref-Cb-Zf/ref_compara")
archivos <- list.files(pattern = "formato_sumstats*")
for (file in archivos){
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_sumstats_(.*)\\..*", "\\1", file)
    temp_dataset$Experiment <- filename
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
    filename <- gsub(".*formato_sumstats_(.*)\\..*", "\\1", file)
    dataset$Experiment <- filename
  }
}

ggplot(dataset) + geom_point(aes(x=BP, y=SmoothedPi, colour=Experiment), alpha=I(1/10)) + facet_grid(. ~ Experiment)
rm(list=ls())
