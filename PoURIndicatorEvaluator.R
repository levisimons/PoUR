rm(list=ls())
require(ggplot2)
require(ggVennDiagram)
wd <- "~/Desktop/eDNA-metadata/"
#wd <- "/home1/alsimons/PoUR"
setwd(wd)

#Input taxonomic level to aggregate on.  Rank 7 is the most resolved, and Rank 1 is the least.
rank=7

#Get random forest model evaluations for the detection of taxa against environmental variables.
SEDIFiles <- list.files(path=wd,pattern=paste('RFEvaluation(.*?)Round(.*?)Rank',rank,sep=""))

#Iterate over sample rounds and primer sets.
#Primers
Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
#Aggregate model evaluations.
RFEvaluation <- data.frame()
for(i in 1:3){
  for(Primer in Primers){
    filename <- paste("RFEvaluation",Primer,"Round",i,"Rank",rank,".txt",sep="")
    tmp <- read.table(filename, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    if(nrow(tmp)>0){
      tmp$Round <- i
      tmp$Primer <- Primer
      RFEvaluation <- rbind(RFEvaluation,tmp) 
    }
  }
}
#Plot histograms of the SEDI scores per taxonomic group, panels are split by primer and sampling round.
ggplot(RFEvaluation,aes(MeanSEDI)) + geom_histogram() + facet_grid(Primer~Round)

#Plot Venn diagram of the number of indicator taxa held in common between sampling rounds.
#Indicator status is defined a taxonomic group with a SEDI score for its random forest model exceeding a certain threshold.
Primer <- "FITS"
SEDIThreshold = 0.5
RFEvaluationFiltered <- RFEvaluation[RFEvaluation$MeanSEDI >= SEDIThreshold & RFEvaluation$Primer==Primer,]
tmp <- with(RFEvaluationFiltered[,c("Taxa","Round")],split(Taxa,Round))
ggVennDiagram(tmp)+labs(title = "Indicator species by sample time point",subtitle=paste("Primer:",Primer,", SEDI >",SEDIThreshold))
                        
