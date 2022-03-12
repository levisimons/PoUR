rm(list=ls())
require(ggplot2)
require(ggVennDiagram)
require(plyr)
require(dplyr)
require(tidyr)
require(reshape2)
require(viridis)
require(RColorBrewer)

wd <- "~/Desktop/eDNA-metadata/"
#wd <- "/home1/alsimons/PoUR"
setwd(wd)

#Input taxonomic level to aggregate on.  Rank 7 is the most resolved, and Rank 1 is the least.
rank=7
TaxonomicRanks <- c("kingdom","phylum","class","order","family","genus","species")

#Get random forest model evaluations for the detection of taxa against environmental variables.
#These were generated here: https://github.com/levisimons/PoUR/blob/main/PoURIndicators.R
EvaluationFiles <- list.files(path=wd,pattern=paste('RFEvaluation(.*?)Round(.*?)Rank',rank,sep=""))

#Iterate over sample rounds and primer sets.  Aggregate random forest SDM evaluations.
#This includes accuracy and relative importance of variables.
if(length(list.files(path=wd,pattern="RFEvaluationAllRoundsSpecies.txt"))<1){
  #Primers
  Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
  #Aggregate model evaluations.
  RFEvaluation <- data.frame()
  for(i in 1:3){
    for(Primer in Primers){
      filename <- paste("RFEvaluation",Primer,"Round",i,"Rank",rank,".txt",sep="")
      if(filename %in% EvaluationFiles){
        tmp <- read.table(filename, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
        if(nrow(tmp)>0){
          tmp$Round <- i
          tmp$Primer <- Primer
          RFEvaluation <- rbind(RFEvaluation,tmp) 
        }
      }
    }
  }
  write.table(RFEvaluation,"RFEvaluationAllRoundsSpecies.txt",quote=FALSE,sep="\t",row.names = FALSE)
}

#Read in summarized random forest models evaluations.
RFEvaluation <- read.table("RFEvaluationAllRoundsSpecies.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#
RankTest <- RFEvaluation[,!(colnames(RFEvaluation) %in% c("Taxa","MeanSEDI","sdSEDI","MeanTSS","sdTSS","Round","Primer"))]
for(i in 1:nrow(RankTest)){
  RankTest[i,] <- rank(RankTest[i,])
}
#Create data frame for ploting relative rank importances for random forest variables.
RankTest <- cbind(RFEvaluation[,c("Round","Primer")],RankTest)
RankTestTmp <- melt(RankTest,id.var=c("Round","Primer"))

##Determine relative rank importance of variables in each species random forest model.
#Test if the distribution of these relative rank importances, within each primer, are significantly different between sampling rounds.
#Read in summarized random forest models evaluations.
RFEvaluation <- read.table("RFEvaluationAllRoundsSpecies.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
RankTest <- RFEvaluation[,!(colnames(RFEvaluation) %in% c("Taxa","MeanSEDI","sdSEDI","MeanTSS","sdTSS","Round","Primer"))]
for(i in 1:nrow(RankTest)){
  RankTest[i,] <- rank(RankTest[i,])
}
RankTest <- cbind(RFEvaluation[,c("Round","Primer")],RankTest)
KWTestSummary <- data.frame()
Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
TestVars <- colnames(RankTest[,!(colnames(RankTest) %in% c('Primer','Round'))])
for(Primer in Primers){
  RankTestSubset <- RankTest[RankTest$Primer==Primer,]
  for(TestVar in TestVars){
    KWTest <- kruskal.test(RankTestSubset[,TestVar]~RankTestSubset[,"Round"])
    KWTestRow <- data.frame(matrix(nrow=1,ncol=4))
    colnames(KWTestRow) <- c("Primer","Variable","p","df")
    KWTestRow$Primer <- Primer
    KWTestRow$Variable <- TestVar
    KWTestRow$p <- KWTest$p.value
    KWTestRow$df <- KWTest$parameter
    KWTestSummary <- rbind(KWTestSummary,KWTestRow)
  }
}
#Write out Kruskal-Wallis test results.
write.table(KWTestSummary,"PoURKWTests.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Create boxplots of the relative rank importances of variables for each random forest variable.
#Split data by primer and then make facets based on sampling round.
#Color variable labels so that black is significant, and red is not.  Using a Kruskal-Wallis test of each variable against sampling round.
RankTest <- dplyr::left_join(RankTestTmp,KWTestSummary[,c("Primer","Variable","p")],by=c("Primer","variable"="Variable"))
RankTest <- RankTest[!duplicated(RankTest),]
Primer <- "CO1"
RankTestPrimerSubset <- RankTest[RankTest$Primer==Primer,]
tmp <- RankTestPrimerSubset[,c("variable","p")]
tmp <- tmp[!duplicated(tmp),]
StatusColor <- ifelse(tmp$p <= 0.05,"black","red")
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(RankTestPrimerSubset$variable)))
ggplot(data=RankTestPrimerSubset,aes(x=variable,y=value,fill=factor(variable)))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=T)+
  facet_grid(rows=vars(Round),labeller=label_both)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color=StatusColor), legend.position = "none")+
  labs(title=paste("Rank importance of variables in predicting the presence of species\nclassified using a",Primer,"primer"),y = "Rank importance")

#Run Dunn tests to see which sampling rounds are driving the significant differences between
#relative rank importance distributions which are significantly different between sampling rounds.
KWTestSignificant <- KWTestSummary[KWTestSummary$p <= 0.05,]
RFEvaluation <- read.table("RFEvaluationAllRoundsSpecies.txt",header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
RankTest <- RFEvaluation[,!(colnames(RFEvaluation) %in% c("Taxa","MeanSEDI","sdSEDI","MeanTSS","sdTSS","Round","Primer"))]
for(i in 1:nrow(RankTest)){
  RankTest[i,] <- rank(RankTest[i,])
}
RankTest <- cbind(RFEvaluation[,c("Round","Primer")],RankTest)
DunnTestSummary <- data.frame()
for(i in 1:nrow(KWTestSignificant)){
  RankTestSubset <- RankTest[RankTest$Primer==KWTestSignificant[i,"Primer"],]
  DunnTest <- dunn.test::dunn.test(RankTestSubset[,c(KWTestSignificant[i,"Variable"],"Round")],method="bonferroni")
  DunnTestRow <- data.frame(matrix(nrow=1,ncol=6))
  colnames(DunnTestRow) <- c("Primer","Variable","chi2","Z","P.adjusted","comparisons")
  DunnTestRow$Primer <- KWTestSignificant[i,"Primer"]
  DunnTestRow$Variable <- KWTestSignificant[i,"Variable"]
  DunnTestRow$chi2 <- DunnTest$chi2
  DunnTestRow$Z <- DunnTest$Z
  DunnTestRow$P.adjusted <- DunnTest$P.adjusted
  DunnTestRow$comparisons <- DunnTest$comparisons
  DunnTestSummary <- rbind(DunnTestSummary,DunnTestRow)
}
#Write out Dunn test results.
write.table(DunnTestSummary,"PoURDunnTests.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Plot histograms of the SEDI scores per taxonomic group, panels are split by primer and sampling round.
ggplot(RFEvaluation,aes(MeanSEDI)) + geom_histogram() + facet_grid(Primer~Round)
#Plot histograms of the TSS scores per taxonomic group, panels are split by primer and sampling round.
ggplot(RFEvaluation,aes(MeanTSS)) + geom_histogram() + facet_grid(Primer~Round)

#Plot Venn diagram of the number of indicator taxa held in common between sampling rounds.
#Indicator status is defined a taxonomic group with a SEDI score for its random forest model exceeding a certain threshold.
Primer <- "FITS"
SEDIThreshold = 0.5
RFEvaluationFilteredSEDI <- RFEvaluation[RFEvaluation$MeanSEDI >= SEDIThreshold & RFEvaluation$Primer==Primer,]
tmp <- with(RFEvaluationFilteredSEDI[,c("Taxa","Round")],split(Taxa,Round))
if(SEDIThreshold>=0){
  ggVennDiagram(tmp)+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position="bottom")+labs(title=paste(TaxonomicRanks[rank],"by sample time point"),subtitle=paste("Primer:",Primer,", SEDI >",SEDIThreshold)) + scale_fill_gradient(low = "yellow", high = "blue")
} else{
  ggVennDiagram(tmp)+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position="bottom")+labs(title=paste(TaxonomicRanks[rank],"by sample time point"),subtitle=paste("Primer:",Primer,", No SEDI cut-off")) + scale_fill_gradient(low = "yellow", high = "blue")
}

#Plot Venn diagram of the number of indicator taxa held in common between sampling rounds.
#Indicator status is defined a taxonomic group with a TSS score for its random forest model exceeding a certain threshold.
Primer <- "16S"
TSSThreshold = 0.5
RFEvaluationFilteredTSS <- RFEvaluation[RFEvaluation$MeanTSS >= TSSThreshold & RFEvaluation$Primer==Primer,]
tmp <- with(RFEvaluationFilteredTSS[,c("Taxa","Round")],split(Taxa,Round))
if(TSSThreshold>=0){
  ggVennDiagram(tmp)+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position="bottom")+labs(title=paste(TaxonomicRanks[rank],"by sample time point"),subtitle=paste("Primer:",Primer,", TSS >",TSSThreshold)) + scale_fill_gradient(low = "yellow", high = "blue")
} else{
  ggVennDiagram(tmp)+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position="bottom")+labs(title=paste(TaxonomicRanks[rank],"by sample time point"),subtitle=paste("Primer:",Primer,", No TSS cut-off")) + scale_fill_gradient(low = "yellow", high = "blue")
}

#Evaluate eDNA results by their traditional observation score (TOS), whether they were used in
#random forest modeling, and whether their random forest models were accurate.
#Determine if indicator species have been observed in GBIF on a local or regional basis.
#Do this for species selected as indicators using either their TSS or SEDI scores for their SDMs. 
#Read in local species list.
LocalSpecies <- read.table("LACounty.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
LocalSpecies <- LocalSpecies[!is.na(LocalSpecies$numberOfOccurrences),]
LocalSpecies[LocalSpecies==""] <- NA
#Assign weight value for species, genera, and families.
LocalSpecies <- LocalSpecies %>% dplyr:: mutate(TaxonWeight = case_when(taxonRank=="SPECIES" ~ 4, taxonRank=="GENUS" ~ 2, taxonRank=="FAMILY" ~ 1, !(taxonRank %in% c("SPECIES","GENUS","FAMILY"))~0))
#Read in regional species list.
RegionalSpecies <- read.table("WesternNorthAmerica.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
RegionalSpecies <- RegionalSpecies[!is.na(RegionalSpecies$numberOfOccurrences),]
RegionalSpecies[RegionalSpecies==""] <- NA
#Assign weight value for species, genera, and families.
RegionalSpecies <- RegionalSpecies %>% dplyr:: mutate(TaxonWeight = case_when(taxonRank=="SPECIES" ~ 4, taxonRank=="GENUS" ~ 2, taxonRank=="FAMILY" ~ 1, !(taxonRank %in% c("SPECIES","GENUS","FAMILY"))~0))
#Definite taxonomic levels.
TaxonomicLevels <- c("superkingdom","phylum","class","order","family","genus","species")
SEDIThreshold <- 0.5
TSSThreshold <- 0.5
SampleInputTotal <- data.frame()
for(Primer in Primers){
  print(Primer)
  #Get samples classified using a particular primer.
  SampleFiles <- list.files(path=wd,pattern=paste("_physeq_",Primer,"_dct_noblanks_min20p",sep=""))
  #Read in sample data classified using a selected primer. Select taxa from all sampling rounds.
  SampleInputSummary <- data.frame()
  for(SampleFile in SampleFiles){
    SampleInput <- read.table(SampleFile, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    SampleInput$File <- SampleFile
    SampleInput <- SampleInput[,c("File","sum.taxonomy")]
    SampleInputSummary <- rbind(SampleInputSummary,SampleInput)
  }
  #Separate taxonomies.
  SampleInputSummary <- tidyr::separate(data=SampleInputSummary,col=sum.taxonomy,sep=";",into=TaxonomicLevels,remove=F)
  #Clean up data.
  SampleInputSummary[SampleInputSummary=="NA" | SampleInputSummary==""] <- NA
  #Count eDNA taxonomic resolution and weigh them.
  #Species = 4, Genus = 2, Family = 1.  Everything else = 0.
  SampleInputSummary$eDNAResolutionWeight <- 1*as.numeric(!is.na(SampleInputSummary$family))+2*as.numeric(!is.na(SampleInputSummary$genus))+4*as.numeric(!is.na(SampleInputSummary$species))
  #Determine if each taxon was used for random forest modeling, and if it counts as accurate.
  RFEvaluationSubset <- RFEvaluation[RFEvaluation$Primer==Primer,]
  RFEvaluationSubset$PresentRandomForest <- 1
  RFEvaluationSubset <- RFEvaluationSubset[,c("Taxa","PresentRandomForest","MeanSEDI","MeanTSS")]
  RFEvaluationSubset <- RFEvaluationSubset[!duplicated(RFEvaluationSubset),]
  RFEvaluationSubset$Accurate_TSS <- ifelse(RFEvaluationSubset$MeanTSS >= TSSThreshold,1,0)
  RFEvaluationSubset$Accurate_SEDI <- ifelse(RFEvaluationSubset$MeanSEDI >= SEDIThreshold,1,0)
  SampleInputSummary <- dplyr::left_join(SampleInputSummary,RFEvaluationSubset,by=c("species"="Taxa"))
  SampleInputSummary$PresentRandomForest[is.na(SampleInputSummary$PresentRandomForest)] <- 0
  #Check if any of the eDNA reads show up in the local set of GBIF family observations.
  SampleInputSummary$LocalFamilyPresentGBIF <- as.numeric(lapply(SampleInputSummary$family,is.element,unique(na.omit(LocalSpecies$family))))
  #Check if any of the eDNA reads show up in the regional set of GBIF family observations.
  SampleInputSummary$RegionalFamilyPresentGBIF <- as.numeric(lapply(SampleInputSummary$family,is.element,unique(na.omit(RegionalSpecies$family))))
  #Check if any of the eDNA reads show up in the local set of GBIF genus observations.
  SampleInputSummary$LocalGenusPresentGBIF <- as.numeric(lapply(SampleInputSummary$genus,is.element,unique(na.omit(LocalSpecies$genus))))
  #Check if any of the eDNA reads show up in the regional set of GBIF family observations.
  SampleInputSummary$RegionalGenusPresentGBIF <- as.numeric(lapply(SampleInputSummary$genus,is.element,unique(na.omit(RegionalSpecies$genus))))
  #Check if any of the eDNA reads show up in the local set of GBIF species observations.
  SampleInputSummary$LocalSpeciesPresentGBIF <- as.numeric(lapply(SampleInputSummary$species,is.element,unique(na.omit(LocalSpecies$species))))
  #Check if any of the eDNA reads show up in the regional set of GBIF family observations.
  SampleInputSummary$RegionalSpeciesPresentGBIF <- as.numeric(lapply(SampleInputSummary$species,is.element,unique(na.omit(RegionalSpecies$species))))
  #Assign local and regional TOS scores for GBIF results.
  SampleInputSummary$LocalTOS <- (1*SampleInputSummary$LocalFamilyPresentGBIF+2*SampleInputSummary$LocalGenusPresentGBIF+4*SampleInputSummary$LocalSpeciesPresentGBIF)/SampleInputSummary$eDNAResolutionWeight
  SampleInputSummary$RegionalTOS <- (1*SampleInputSummary$RegionalFamilyPresentGBIF+2*SampleInputSummary$RegionalGenusPresentGBIF+4*SampleInputSummary$RegionalSpeciesPresentGBIF)/SampleInputSummary$eDNAResolutionWeight
  SampleInputTotal <- rbind(SampleInputTotal,SampleInputSummary)
}
SampleInputTotal$LocalTOS[is.na(SampleInputTotal$LocalTOS)] <- 0
SampleInputTotal$RegionalTOS[is.na(SampleInputTotal$RegionalTOS)] <- 0
#Output summary to a file.
write.table(SampleInputTotal,"PoUReDNAEvaluations.txt",quote=FALSE,sep="\t",row.names = FALSE)

##Generate fasta files for each primer.  These are use as input in BLAST-n in order to check
#for likeliest matches between our ASVs and taxonomically identified sequences.
require(seqRFLP)
Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
#Aggregate model evaluations.
RFEvaluation <- data.frame()
for(i in 1:3){
  for(Primer in Primers){
    filename <- paste("RFEvaluation",Primer,"Round",i,"Rank",rank,".txt",sep="")
    if(filename %in% EvaluationFiles){
      tmp <- read.table(filename, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
      if(nrow(tmp)>0){
        tmp$Round <- i
        tmp$Primer <- Primer
        RFEvaluation <- rbind(RFEvaluation,tmp) 
      }
    }
  }
}
TSSThreshold = 0.5
SEDIThreshold = 0.5
TaxonomicLevel <- "species"
for(Primer in Primers){
  DetailedTaxonomies <- list.files(path=wd,pattern=paste("R(.*?)",Primer,"_ASV_taxonomy_detailed.txt",sep=""))
  SequenceInput <- data.frame()
  for(DetailedTaxonomy in DetailedTaxonomies){
    SequenceInputTmp <- read.table(DetailedTaxonomy, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    #Extract sample round from filename.
    Round <- as.numeric(gsub("R([0-9]+).*$","\\1",DetailedTaxonomy))
    #Remove sequences with missing taxonomy
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$taxonomy!="not_found",]
    #Remove missing sequence rows.
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$sequence!="",]
    #Remove sequences with an assigned species confidence below 80%.
    SequenceInputTmp$speciesBayesianConfidence <- as.numeric(gsub(".*species:([0-9]+).*;$","\\1",SequenceInputTmp$taxonomy_confidence))
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$speciesBayesianConfidence >= 80 & !is.na(SequenceInputTmp$speciesBayesianConfidence),]
    SequenceInputTmp$Round <- Round#Assign sampling round number to sequence data.
    #Select sequences which show up as potential indicators using a TSS score cutoff.
    SequenceInputTmp <- SequenceInputTmp[,c(paste(Primer,"_seq_number",sep=""),"Round","input_sequence_length","sequence","taxonomy","accessions")]
    #Clean up sequence data.
    SequenceInputTmp <- SequenceInputTmp[!duplicated(SequenceInputTmp),]
    SequenceInputTmp <- tidyr::separate(data=SequenceInputTmp,col=taxonomy,sep=";",into=c("superkingdom","phylum","class","order","family","genus","species"),remove=F)
    SequenceInputTmp[,TaxonomicLevel] <- gsub(paste(TaxonomicLevel,":",sep=""),"",SequenceInputTmp[,TaxonomicLevel])
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$taxonomy!="Unclassified" & !is.na(SequenceInputTmp$sequence) & !is.na(SequenceInputTmp[,TaxonomicLevel]),]
    SequenceInputTmp[SequenceInputTmp=="" | SequenceInputTmp=="NA"] <- NA
    SequenceInputTmp$FASTA_ID <- SequenceInputTmp[,colnames(SequenceInputTmp[grepl("_seq_number",colnames(SequenceInputTmp))])]
    #Get taxa which have accurate random forest models.
    RFEvaluationFiltered <- RFEvaluation[RFEvaluation$MeanTSS >= TSSThreshold | RFEvaluation$MeanSEDI >= SEDIThreshold,]
    RFEvaluationFiltered <- RFEvaluationFiltered[RFEvaluationFiltered$Primer==Primer & RFEvaluationFiltered$Round==Round,]
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$species %in% unique(RFEvaluationFiltered$Taxa),]
    SequenceInputTmp <- SequenceInputTmp[!duplicated(SequenceInputTmp),]
    SequenceInput <- rbind(SequenceInput,SequenceInputTmp)
  }
  #Export ASVs which are associated with taxa with accurate random forest models to FASTA files for further processing.
  SequenceFASTA <- dataframe2fas(SequenceInput[,c("FASTA_ID","sequence")])
  write.fasta(SequenceFASTA,file=paste("AccurateASVs_",Primer,".fasta",sep=""))
}

##Generate files giving the name, and common name, for environmentally sensitive ASVs.
#Include the bitscore and evalue for the best matches.
require(data.table)
require(taxize)
options(ENTREZ_KEY="3f9e29c7c03842f72cf7523e34390d9f2208")
#Determine if any primers have already been analyzed.
Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
CompletedFiles <- list.files(path=wd,pattern='AccurateASVAssessments_(.*?).txt')
CompletedPrimers <- lapply(X=CompletedFiles, FUN = function(t) gsub(pattern = "AccurateASVAssessments_", replacement = "", x = t, fixed = TRUE))
CompletedPrimers <- lapply(X=CompletedPrimers,FUN = function(t) gsub(pattern = ".txt", replacement = "", x = t, fixed = TRUE))
Primers <- Primers[!(Primers %in% CompletedPrimers)]

for(Primer in Primers){
  #Read in BLAST output from environmental sensitive ASVs by primer.
  BLASTInput <- read.table(paste("AccurateASVs_",Primer,".txt",sep=""), header=F, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  colnames(BLASTInput) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  BLASTFiltered <- BLASTInput %>% dplyr::group_by(qseqid) %>% dplyr::filter(evalue==min(evalue),bitscore==max(bitscore))
  BLASTFiltered <- BLASTFiltered[,c("qseqid","sseqid","evalue","bitscore")]
  
  SequenceIDs <- unique(BLASTFiltered$sseqid)
  if(length(list.files(path=wd,pattern=paste("AccurateASVAssessments_",Primer,".txt",sep="")))>0){
    CheckInput <- read.table(paste("AccurateASVAssessments_",Primer,".txt",sep=""),header=T, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    SequenceIDs <- unique(BLASTFiltered$sseqid[!(BLASTFiltered$sseqid %in% CheckInput$sseqid)])
  }
  
  
  #Determine the name, and common name, for the BLAST results for the most environmentally sensitive ASVs.
  for(SequenceID in SequenceIDs){
    TaxTmp <- data.frame(matrix(nrow=1,ncol=5))
    colnames(TaxTmp) <- c("uid","sseqid","name","rank","CommonName")
    TaxTmp$sseqid <- SequenceID
    tmp <-  as.data.frame(ncbi_get_taxon_summary(id=as.character(genbank2uid(id=SequenceID))))
    Sys.sleep(0.2)
    if(nrow(tmp)>0){
      TaxTmp$uid <- tmp$uid
      TaxTmp$name <- tmp$name
      TaxTmp$rank <- tmp$rank
      TaxID <- sci2comm(get_uid(tmp$name),db="ncbi")
      if(length(TaxID)>0){
        if(length(TaxID[[1]])<1){
          TaxTmp$CommonName <- NA
        } else{
          TaxTmp$CommonName <- TaxID[[1]]
        }
      } else{
        TaxTmp$CommonName <- NA
      } 
    } else{
      TaxTmp$uid <- NA
      TaxTmp$name <- NA
      TaxTmp$rank <- NA
      TaxTmp$CommonName <- NA
    }
    print(paste(Primer,TaxTmp$uid,TaxTmp$sseqid,TaxTmp$name,TaxTmp$rank,TaxTmp$CommonName))
    #Merge environmentally sensitive ASVs with the bitscore and evalues of the most likely BLAST matches.
    SequencesWithBLAST <- dplyr::left_join(BLASTFiltered[BLASTFiltered$sseqid==SequenceID,],TaxTmp,by=c("sseqid"))
    #Output summary of the best sequence matches, for environmentally sensitive ASVs, to a file.
    if(length(list.files(path=wd,pattern=paste("AccurateASVAssessments_",Primer,".txt",sep="")))==0){
      write.table(SequencesWithBLAST,paste("AccurateASVAssessments_",Primer,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
    } else{
      write.table(SequencesWithBLAST,paste("AccurateASVAssessments_",Primer,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE,col.names=F,append=T)
    }
  }
}

##Generate text files containing the full taxonomies of the input ASVs,
#along with the BLAST results of the best ASV matches which also correspond to potential indicator species.
RFEvaluation <- read.table("RFEvaluationAllRoundsSpecies.txt",header=T, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
TSSThreshold = 0.5
SEDIThreshold = 0.5
TaxonomicLevel <- "species"
Primers <- c("16S", "18S","CO1","VERT12S","PITS","FITS")
for(Primer in Primers){
  DetailedTaxonomies <- list.files(path=wd,pattern=paste("R(.*?)",Primer,"_ASV_taxonomy_detailed.txt",sep=""))
  SequenceInput <- data.frame()
  for(DetailedTaxonomy in DetailedTaxonomies){
    SequenceInputTmp <- read.table(DetailedTaxonomy, header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
    #Extract sample round from filename.
    Round <- as.numeric(gsub("R([0-9]+).*$","\\1",DetailedTaxonomy))
    #Remove sequences with missing taxonomy
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$taxonomy!="not_found",]
    #Remove missing sequence rows.
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$sequence!="",]
    #Remove sequences with an assigned species confidence below 80%.
    SequenceInputTmp$speciesBayesianConfidence <- as.numeric(gsub(".*species:([0-9]+).*;$","\\1",SequenceInputTmp$taxonomy_confidence))
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$speciesBayesianConfidence >= 80 & !is.na(SequenceInputTmp$speciesBayesianConfidence),]
    SequenceInputTmp$Round <- Round#Assign sampling round number to sequence data.
    #Select sequences which show up as potential indicators using a TSS score cutoff.
    SequenceInputTmp <- SequenceInputTmp[,c(paste(Primer,"_seq_number",sep=""),"Round","input_sequence_length","sequence","taxonomy","accessions")]
    #Clean up sequence data.
    SequenceInputTmp <- SequenceInputTmp[!duplicated(SequenceInputTmp),]
    SequenceInputTmp <- tidyr::separate(data=SequenceInputTmp,col=taxonomy,sep=";",into=c("superkingdom","phylum","class","order","family","genus","species"),remove=F)
    SequenceInputTmp[,c("superkingdom","phylum","class","order","family","genus","species")] <- as.data.frame(apply(SequenceInputTmp[,c("superkingdom","phylum","class","order","family","genus","species")],2,function(y) gsub(".*:","",y)))
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$taxonomy!="Unclassified" & !is.na(SequenceInputTmp$sequence) & !is.na(SequenceInputTmp[,TaxonomicLevel]),]
    SequenceInputTmp[SequenceInputTmp=="" | SequenceInputTmp=="NA"] <- NA
    SequenceInputTmp$FASTA_ID <- SequenceInputTmp[,colnames(SequenceInputTmp[grepl("_seq_number",colnames(SequenceInputTmp))])]
    #Get taxa which have accurate random forest models.
    RFEvaluationFiltered <- RFEvaluation[RFEvaluation$MeanTSS >= TSSThreshold | RFEvaluation$MeanSEDI >= SEDIThreshold,]
    RFEvaluationFiltered <- RFEvaluationFiltered[RFEvaluationFiltered$Primer==Primer & RFEvaluationFiltered$Round==Round,]
    SequenceInputTmp <- SequenceInputTmp[SequenceInputTmp$species %in% unique(RFEvaluationFiltered$Taxa),]
    SequenceInputTmp <- SequenceInputTmp[!duplicated(SequenceInputTmp),]
    SequenceInput <- rbind(SequenceInput,SequenceInputTmp)
  }
  BLASTAssessment <- read.table(paste('AccurateASVAssessments_',Primer,'.txt',sep=""),header=T, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  SequencesWithBLAST <- dplyr::left_join(SequenceInput,BLASTAssessment,by=c("FASTA_ID"="qseqid"))
  write.table(SequencesWithBLAST,paste("SequencesWithBLAST_",Primer,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}

##Generate files adding the full taxonomy of each indicator ASV, along with the full taxonomies
#of the best BLAST-n match.
CompletedFiles <- list.files(path=wd,pattern='SequencesWithBLAST_(.*?).txt')
for(CompletedFile in CompletedFiles){
  SequencesWithBLAST <- read.table(CompletedFile,header=T, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
  SequencesWithBLAST <- SequencesWithBLAST[!is.na(SequencesWithBLAST$uid),]
  FullTaxonomy <- data.frame()
  for(taxon in unique(SequencesWithBLAST$uid)){
    tmp <- classification(taxon,db='ncbi')
    tmp <- as.data.frame(tmp[1])
    names(tmp) <-  gsub(paste("X",taxon,".",sep=""),"",names(tmp))
    tmp <- tmp[tmp$rank %in% c("species","genus","family","order","class","phylum","superkingdom"),]
    ranks <- as.data.frame(c("species","genus","family","order","class","phylum","superkingdom"))
    colnames(ranks) <- c("rank")
    tmp <- merge(tmp,ranks,all.y=T)
    tmp$id <- NULL
    tmp <- as.data.frame(t(tmp))
    colnames(tmp) <- tmp[rownames(tmp)=='rank',]
    tmp <- tmp[!(rownames(tmp) %in% c("rank")),]
    tmp$uid <- taxon
    tmp <- tmp[,c("uid","species","genus","family","order","class","phylum","superkingdom")]
    FullTaxonomy <- rbind(FullTaxonomy,tmp)
    print(paste("Primer: ",Primer," taxon ",nrow(FullTaxonomy)," of ",length(unique(SequencesWithBLAST$uid)),sep=""))
  }
  colnames(FullTaxonomy) <- paste("Matched_",colnames(FullTaxonomy),sep="")
  SequencesWithFullTaxonomy <- dplyr::left_join(SequencesWithBLAST,FullTaxonomy,by=c("uid"="Matched_uid"))
  SequencesWithFullTaxonomy <- SequencesWithFullTaxonomy[,c("FASTA_ID","sseqid","Round","input_sequence_length","evalue","bitscore","uid","rank","CommonName","species","Matched_species","genus","Matched_genus","family","Matched_family","order","Matched_order","class","Matched_class","phylum","Matched_phylum","superkingdom","Matched_superkingdom")]
  write.table(SequencesWithFullTaxonomy,paste("SequencesWithFullTaxonomy_",Primer,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}
