#Set directory and call libraries
library(ggplot2)
library(data.table)
library(scales)

#Import Files
IbdSeqCOdf= read.delim(file = "CO_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqCRdf= read.delim(file = "CR_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqCLMdf= read.delim(file = "CLM_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqCEUdf= read.delim(file = "CEU_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqYRIdf= read.delim(file = "YRI_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqPURdf= read.delim(file = "PUR_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqPELdf= read.delim(file = "PEL_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqMEXdf= read.delim(file = "MEX_allchroms_IBDSeq_FormattedRscript.out")
IbdSeqFINdf= read.delim(file = "FIN_allchroms_IBDSeq_FormattedRscript.out")

AllPopulationsDF=rbind.data.frame(IbdSeqCOdf,IbdSeqCRdf,IbdSeqCEUdf,IbdSeqCLMdf,IbdSeqYRIdf, IbdSeqFINdf, IbdSeqMEXdf, IbdSeqPELdf, IbdSeqPURdf)
names(AllPopulationsDF)[8] = "Population"
AllPopulationsDF$Population = gsub('MEX','MXL',AllPopulationsDF$Population)
AllPopulationsDF$Population = factor(AllPopulationsDF$Population, levels=c("YRI","CEU","FIN","PEL","CLM","CO", "CR","MXL","PUR"))

#Plot Distribution of Original IBD Segments
P1 = ggplot(AllPopulationsDF, aes(x=Len/10^6, fill=Population)) + geom_density(aes(fill=Population)) + scale_fill_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + xlab("Length (Mb)") + theme_bw()  + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) 

#print(P1)

#QC IBD Segments
#Remove IBDSegments less than 3 cM and greater than 20
#Remove IBDSegments if the proportion of IBD segment covered by snps is not within a SD of the mean
AllPopulationsDF_min3cM = subset.data.frame(AllPopulationsDF, AllPopulationsDF$cM > 3 & AllPopulationsDF$cM <= 20)
z = data.table(AllPopulationsDF_min3cM)
z[,ToKeep := abs(AllPopulationsDF_min3cM$Prop_Cov - mean(AllPopulationsDF_min3cM$Prop_Cov)) < sd(AllPopulationsDF_min3cM$Prop_Cov)][ToKeep  == TRUE] #create variable that identifies what to drop or keep 
table(z$ToKeep, z$Population)/length(z$ToKeep) #Proportion that fall within drop or keep for each Population
TrueIBDSeqIBDSegs = subset(z, z$ToKeep == "TRUE")
IBDSeq_Scores = aggregate(as.numeric(TrueIBDSeqIBDSegs$cM), by = list(Population = TrueIBDSeqIBDSegs$Population), FUN = sum)#Calculate length per Population
IBDSeq_Scores$CountsPerPopulation = table(TrueIBDSeqIBDSegs$Population)
sampSize = 30
normalizationConstant = (choose(2*sampSize, 2)) - sampSize #From Nakatsuka et al. 2017 Nature Genetics
IBDSeq_Scores$NormalizedScore = IBDSeq_Scores$x/normalizationConstant

#Make Scores Relative to Finnish
y = grep("FIN", IBDSeq_Scores$Population) #Find them in DF
Relative2FIN = IBDSeq_Scores[y,4] #Grab the Finnish normalized score
IBDSeq_Scores$R2FIN = IBDSeq_Scores$NormalizedScore/Relative2FIN

#Plot new data
TrueIBDSeqIBDSegs$Population = factor(TrueIBDSeqIBDSegs$Population, levels=c("YRI","CEU","FIN","PEL","CLM","CO","CR","MXL","PUR"))

P2 = ggplot(TrueIBDSeqIBDSegs, aes(x=Len/10^6, y=SNP_Count)) + geom_point(aes(colour=Population)) + scale_colour_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + theme_bw() + xlab("Length (Mb)") + ylab("Total SNPs in IBD Segment") + ggtitle("Physical Distance vs Total SNPs in IBD Segment")  + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) 

P3 = ggplot(TrueIBDSeqIBDSegs, aes(x=cM, y=SNP_Count)) + geom_point(aes(colour=Population)) + scale_colour_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + theme_bw() + xlab("Length (cM)") + ylab("Total SNPs in IBD Segment") + ggtitle("Genetic Distance vs Total SNPs in IBD Segment")  + guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20))

P4 = ggplot(TrueIBDSeqIBDSegs, aes(x=Len/10^6, y=cM)) + geom_point(aes(colour=Population)) + scale_colour_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + theme_bw() + xlab("Length (Mb)") + ylab("Length (cM)") + ggtitle("Physical Distance vs Genetic Distance")  + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) 

P5 = ggplot(TrueIBDSeqIBDSegs, aes(x=Population, y=cM, fill=Population)) + geom_violin(scale = "width") + geom_boxplot(width=0.1) + scale_fill_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + theme_bw() + xlab("Population") + ylab("Length IBD Segment (cM)") + ggtitle("IBD Segements per Population")  + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) +  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) 

pdf("Scenarios_CountNeutralSites_Aug17.pdf", width=15.25, height=8.0)

print(P1)
print(P2)
print(P3)
print(P4)
print(P5)

dev.off()


BeforeFilt = ggplot(AllPopulationsDF, aes(x=Len/10^6, y=SNP_Count)) + geom_point(aes(colour=Population)) + scale_colour_manual(values=c(CO="goldenrod2", CR="coral2", CEU="cornflowerblue", CLM="tan4", YRI="darkseagreen3" , CO_RmHighAfrAnc="darkorchid4", PEL="deeppink4", PUR="darkviolet", FIN="cyan4", MXL="chocolate1")) + theme_bw() + xlab("Length of IBD segment (Mb)") + ylab("Total SNPs in IBD Segment") + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold", hjust=0.5), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1)) 

AfterFilt = P2 + ggtitle(NULL)

pdf("FilteringIBDSegments.pdf", width=15.25, height=8.0)

grid.arrange(BeforeFilt + theme(legend.position = "none") + labs(x=NULL, y=NULL) + scale_y_continuous(labels = scales::scientific), AfterFilt + labs(x=NULL, y=NULL), left=textGrob("Count of SNPs", gp=gpar(fontface="bold",fontsize=20), rot=90) , bottom=textGrob("Length of IBD segment (Mb)", gp=gpar(fontface="bold",fontsize=20), vjust=0.01), nrow =1)

dev.off()

#Testing Correlation Between physical and genetic distance
PhysVsGen = lm(TrueIBDSeqIBDSegs$cM~TrueIBDSeqIBDSegs$Len)
#summary(PhysVsGen)

#Permutation test for difference in total length aka score of ibd tracts between populations
#pop1 is the test population 
#pop2 is the ref population
#popnames should be supplied in quotes
PermutationTest = function(dataFrame, pop1Name, pop2Name, numberPerms){
  subdata = subset.data.frame(dataFrame, dataFrame$Population == pop1Name | dataFrame$Population == pop2Name)
  y = subdata$cM
  group = as.character(subdata$Population)
  testStat <- function(w, g) sum(as.numeric(w[g == pop1Name])) - sum(as.numeric(w[g == pop2Name])) #test whether there is a difference in sum of Length of IBD tracts in CR vs FIN
  observedStat = testStat(y, group)
  permutations = sapply(1 : numberPerms, function(i) testStat(y, sample(group)))
  pvalue = mean(permutations > observedStat) #pvalue 
  pvalue[pvalue == 0] = 1/numberPerms #generate pvalue when there are no perms greater than observed
  bins = seq(1:numberPerms)
  dfPerms = cbind.data.frame(bins,permutations)
  return(list(dfPerms, observedStat, pvalue))
}

#Conduct permutation test
OutputCO = PermutationTest(TrueIBDSeqIBDSegs, "CO", "FIN", 10000)
OutputCR = PermutationTest(TrueIBDSeqIBDSegs, "CR", "FIN", 10000)
OutputCLM = PermutationTest(TrueIBDSeqIBDSegs, "CLM", "FIN", 10000)
OutputPUR = PermutationTest(TrueIBDSeqIBDSegs, "PUR", "FIN", 10000)
OutputPEL = PermutationTest(TrueIBDSeqIBDSegs, "PEL", "FIN", 10000)
OutputMXL = PermutationTest(TrueIBDSeqIBDSegs, "MXL", "FIN", 10000)
OutputYRI = PermutationTest(TrueIBDSeqIBDSegs, "YRI", "FIN", 10000)
OutputCEU = PermutationTest(TrueIBDSeqIBDSegs, "CEU", "FIN", 10000)

dataframePvaluesPerms = t(cbind.data.frame(OutputPUR[[3]], OutputCO[[3]], OutputCR[[3]], OutputCLM[[3]], OutputMXL[[3]], OutputPEL[[3]], OutputYRI[[3]],OutputCEU[[3]]))
 #alternate plotting for CR

COPerm = ggplot(OutputCO[[1]], aes(OutputCO[[1]]$permutations)) + geom_histogram(binwidth = 30, breaks=seq(500, 3000, by =20), col="goldenrod2", fill="white") + geom_vline(xintercept = OutputCO[[2]], col="purple") + theme_bw() + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold"), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + labs(x= "Permutation Score", y= "Count") + ggtitle("Permutations for CO Relative to FIN")
 
CRPerm = ggplot(OutputCR[[1]], aes(OutputCR[[1]]$permutations)) + geom_histogram(binwidth = 30, breaks=seq(500, 3000, by =20), col="coral2", fill="white") + geom_vline(xintercept = OutputCR[[2]], col="purple") + theme_bw() + theme(axis.text.x = element_text(size  = 20), axis.text.y = element_text(size  = 20), plot.title=element_text(size=26, face = "bold"), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=20), legend.text=element_text(size=20)) + labs(x= "Permutation Score", y= "Count") + ggtitle("Permutations for CR Relative to FIN")
                                                                                                     
