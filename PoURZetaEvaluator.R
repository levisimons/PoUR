rm(list=ls())
require(ggplot2)
require(vegan)
require(viridis)

wd <- "~/Desktop/eDNA-metadata/"
#wd <- "/home1/alsimons/PoUR"
setwd(wd)

#Read in results of the relative likelihood of community assembly models
#generated here: https://github.com/levisimons/PoUR/blob/main/PoURZeta.R
zetaAnalysis <- read.table("PoURZetaAssembly.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")

#Set factor variables.
zetaAnalysis$Primer <- as.factor(zetaAnalysis$Primer)
zetaAnalysis$SampleRound <- as.factor(zetaAnalysis$SampleRound)
zetaAnalysis$Site <- as.factor(zetaAnalysis$Site)

#Plot trend in AIC(Exponential)-AIC(Power-law) versus the global human modification.
ZetaPlot <- ggplot(zetaAnalysis,aes(x=gHM,y=deltaAIC))+geom_point()+theme(text = element_text(size=25))+geom_smooth(method=glm, aes(fill=deltaAIC))
ZetaPlot+xlab("Global Human Modification")+ylab(expression(Delta))

#Estimate contributions of variation to our model of AIC(Exponential)-AIC(Power-law) versus the global human modification.
adonis(deltaAIC ~gHM+Primer+SampleRound+Site,data=zetaAnalysis,permutations=1000,method="manhattan")

#Read in results related to factors influencing variation in zeta_4
#generated here: https://github.com/levisimons/PoUR/blob/main/PoURZeta.R
#Generate histogram, split by primer and sampling round, for what contributes to the
#observed variations in zeta_4
zetaVarAnalysis <- read.table("PoURZetaFactors.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
tmp1 <- zetaVarAnalysis[c("Primer","SampleRound","Variable","VariationFromDistance")]
colnames(tmp1) <- c("Primer","SampleRound","Variable","Variation")
tmp1$Variable_type <- "Distance"
tmp2 <- zetaVarAnalysis[c("Primer","SampleRound","Variable","VariationFromEnvironment")]
colnames(tmp2) <- c("Primer","SampleRound","Variable","Variation")
tmp2$Variable_type <- "Environment"
tmp3 <- zetaVarAnalysis[c("Primer","SampleRound","Variable","VariationUnknown")]
colnames(tmp3) <- c("Primer","SampleRound","Variable","Variation")
tmp3$Variable_type <- "Unknown"
tmp <- rbind(tmp1,tmp2,tmp3)
ggplot(tmp,aes(x=Variation,fill=Variable_type))+geom_histogram(alpha=0.8)+facet_grid(Primer~SampleRound)+scale_fill_manual(values=viridis(3))
