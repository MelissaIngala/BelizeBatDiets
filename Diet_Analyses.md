#Download and load packages
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
packageVersion("phyloseq")
library("ggplot2")
#packageVersion("ggplot2")
library(vegan)
library(microbiome)
library(bipartite)
setwd(dir = '/Users/mingala/Documents/Dissertation/Chapter2/Results 3.23.20/BatDiets3.24/')
#install.packages("spaa")
library(spaa)
library(microbiomeSeq)

#Import Feature table and metadata
motu_table <-read.delim(file = "FINAL_MOTUs.tsv", sep = "\t", header = T, row.names = 1)
taxonomy <- read.delim(file= "FINAL_TAXONOMY.tsv", sep = "\t", header = T, row.names = 1)
metadata<- read.delim(file ="Metadata_diet.txt", sep = "\t", header= T, row.names = 2)

#Create phyloseq object by converting each file into component phyloseq class objects
SAM<-sample_data(metadata)
TAX<-tax_table(as.matrix(taxonomy))
OTU<- otu_table(motu_table, taxa_are_rows=TRUE)
physeq<-merge_phyloseq(OTU, TAX, SAM)

summarize_phyloseq(physeq)
sample_sums(physeq)

############ DECONTAM; REMOVE CONTAMINANTS USING INITIAL DNA CONCENTRATIONS ############
require(ggplot2)
#install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("decontam", version = "3.8")
library(decontam)

#Inspect library sizes
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$DNA_CONC <- sample_sums(physeq)
df <- df[order(df$DNA_CONC),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=DNA_CONC, color= Niche)) + geom_point()

#Detect contaminants by frequency
pruned<- subset_samples(physeq = physeq, DNA_CONC != "NA")
contamdf.freq <- isContaminant(pruned, method="frequency", conc="DNA_CONC")
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

#Remove all contaminants from the dataset and use this new phyloseq obj (noncontam) for all downstream analyses
noncontam <- prune_taxa(!contamdf.freq$contaminant, pruned)

################## SUMMARY PLOTS ################
#First remove samples that have 0 identified reads
phyloseq_richness_filter <- function(physeq, mintaxa = 1){
  sp <- estimate_richness(physeq, measures = "Observed")
  samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
  
  if(length(samples_to_keep) == 0){
    stop("All samples will be removed.\n")  
  }
  
  if(length(samples_to_keep) == nsamples(physeq)){
    cat("All samples will be preserved\n")
    res <- physeq
  }
  
  if(length(samples_to_keep) < nsamples(physeq)){
    res <- prune_samples(samples = samples_to_keep, x = physeq)
  }
  
  return(res)
}  

filtered<-phyloseq_richness_filter(physeq = noncontam, mintaxa = 4)
sample.sums<-as.vector(sample_sums(filtered))
median(sample.sums) #2630
range(sample.sums)

#Transform to relative abundance
relative  = transform_sample_counts(filtered, function(OTU) OTU / sum(OTU))
plot_bar(physeq = relative, fill = "order_name", facet_grid =~kingdom_name, x = "Binomen")

#Top 25 taxa
#Barplot of top 25 OTUs
Top25OTUs = names(sort(taxa_sums(filtered), TRUE)[1:25])
comparetop25 = prune_taxa(Top25OTUs, filtered)
relative  = transform_sample_counts(comparetop25, function(OTU) OTU / sum(OTU))
plot_bar(physeq = relative, fill = "order_name", facet_grid =~kingdom_name, x = "Binomen")

batspecies<-filtered@sam_data$Binomen
dupes<-duplicated(batspecies)
filtered@sam_data$duplicated<-dupes
multiple.spp.obs <- prune_samples(filtered@sam_data$duplicated==TRUE, filtered)
multiple.spp.obs <- prune_samples(sample_sums(multiple.spp.obs)>=1, multiple.spp.obs)

#Evaluate within-species variation
library(microbiomeSeq)
#select bat spp with multiple observations
multispp<-table(filtered@sam_data$Binomen)
multiple.spp.obs<-as.vector(names(which(multispp<2)))
filtered2<-subset_samples(physeq = filtered, Binomen!="Artibeus_intermedius"& Binomen!="Bauerus_dubiaquercus"     
                          & Binomen!= "Chrotopterus_auritus" & Binomen!= "Dermanura_phaeotis"        
                          & Binomen!= "Gardnernycteris_crenulatum" & Binomen!="Mimon_cozumelae"           
                          & Binomen!= "Mormoops_megalophylla" & Binomen!="Myotis_keasyi"             
                          & Binomen!= "Trachops_cirrhosus")
filtered2 <- taxa_level(physeq = filtered2,which_level = "order_name")
relative<- transform_sample_counts(filtered2, function(OTU) OTU / sum(OTU))
p<- plot_taxa(relative, grouping_column = "Binomen", method = 'hellinger', number.taxa = 10, filename = "LCBD_multispeciesobs.csv")
pdf(file = "LCBD_multispp.pdf", width = 10, height = 5)
p
dev.off()

#Agglomerate by species
#Remove confidence from Tax table (messes up agglomeration)
tax_table(filtered) <- tax_table(filtered)[,1:5]
mergedspp<-merge_samples(x = filtered, group = "Binomen")
relative.merged  = transform_sample_counts(mergedspp, function(OTU) OTU / sum(OTU))

#p<-plot_bar(physeq = relative.merged, facet_grid=~kingdom_name)
#pdf(file = "Inv_Pln_Bars.pdf", width = 10, height = 7)
#p + geom_bar(aes(fill=order_name), stat="identity", position="stack")
#dev.off()

#Alpha diversity
Obs.MOTU<-plot_richness(physeq = filtered, color = "Niche", measures = "Observed", x = "Binomen") + geom_boxplot(outlier.color = "black") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
Obs.MOTU

#Check correlation between # of species replicate & # MOTUs, controlling for outliers
nobs<-as.data.frame(multispp)
Obs.MOTU<-estimate_richness(mergedspp, measures = c("Observed", "Shannon"))
Obs.MOTU$Var1<-rownames(Obs.MOTU)
samplesize_testdr<-merge.data.frame(nobs, Obs.MOTU, by = "Var1")

outliers <- boxplot(samplesize_testdr$Observed, plot=TRUE)$out
outliers

samplesize_testdr<-samplesize_testdr[-c(15, 17, 18), ]

#Observed MOTUs
cor.test(samplesize_testdr$Freq, samplesize_testdr$Observed, method = "kendall", exact = FALSE)
plot(samplesize_testdr$Freq, samplesize_testdr$Observed)
#tau=0.33, p = 0.03

#Shannon diversity
cor.test(samplesize_testdr$Freq, samplesize_testdr$Shannon, method = "kendall", exact = FALSE)
#tau = 0.07, p = 0.64

############## DATA TRANSFORMATION AND ORDINATION #################
#Center Log transform the data (don't rarefy)
require(microbiome)
physeq.trans<-microbiome::transform(filtered, "log10")
physeq.trans<-subset_samples(physeq = physeq.trans, Niche!="Control")
physeq.corr<-subset_samples(physeq = physeq.trans, Sample_ID!=52)

funx.ord <- ordinate(
  physeq = physeq.trans, 
  method = "NMDS", 
  distance = "bray"
)

funx.ord.jac <- ordinate(
  physeq = physeq.trans, 
  method = "NMDS", 
  distance = "jaccard"
)

# Calculate distance matrix
study.bray <- phyloseq::distance(physeq.trans, method = "bray")
study.jaccard<- phyloseq::distance(physeq.trans, method = "jaccard")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.trans))

## Homogeneity of dispersion test
beta <- betadisper(study.bray, sampledf$Niche, sqrt.dist = T)
permutest(beta) #p = 0.001
#May be due to differences in dispersion

beta.jac <- betadisper(study.jaccard, sampledf$Niche, sqrt.dist = T)
permutest(beta.jac) #p = 0.001

require(viridis)
pdf(file = "NMDS_Diet.pdf", width = 7, height = 5)
plot_ordination(
  physeq = physeq.trans,
  ordination = funx.ord,
  axes = c(1,2), 
  color = "Niche", 
  title = "NMDS of Dietary MOTUs") +
  scale_color_viridis_d()+
  geom_point(aes(color = Niche), size = 2.5) +
  scale_shape_manual(values = 0:7) +
  theme_bw() +
  stat_ellipse() + geom_jitter()
dev.off()

############# PERMANOVA OF DIET #############
sampledf$Niche<-as.factor(sampledf$Niche)
adonis(study.bray ~ Niche * Binomen + Family, data = sampledf)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Niche      3    2.6846 0.89487  2.3708 0.09014  0.001 ***
#Binomen   20   11.6227 0.58114  1.5396 0.39025  0.001 ***
#Residuals 41   15.4757 0.37746         0.51962           
#Total     64   29.7830                 1.00000    

adonis(study.jaccard ~ Niche * Binomen + Family, data = sampledf)

#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Niche        3     2.322 0.77398  1.8645 0.07585  0.001 ***
  #Binomen   20    11.271 0.56355  1.3576 0.36818  0.001 ***
  #Residuals 41    17.020 0.41512         0.55597           
#Total       64    30.613                 1.00000           
---
################## BIPARTITE NETWORK AND FOOD WEB #####################
######## use package bipartite
#install.packages('bipartite')
require(bipartite)

######## CONVERT PHYLOSEQ MOTU TABLE AND TAX TABLE INTO BIPARTITE INTERACTION MATRIX 
merged_bats_Order<-tax_glom(physeq = mergedspp, taxrank = "order_name")
merged.motus<-as.data.frame(merged_bats_Order@otu_table)
merged.tax<-as.matrix((tax_table(merged_bats_Order)))
merged.tax<-data.frame(merged.tax[,1:2])
colnames(merged.motus)<-merged.tax[,2]

interactions<-as.matrix(merged.motus)

#Visualize web
#pdf(file = 'FinalBipartite Network.pdf', width = 16, height = 13)
plotweb(sortweb(interactions, sort.order="inc"), text.rot = 90)
#dev.off()
visweb(sortweb(interactions, sort.order="inc"))

#Null model comparison
#To test network metrics (e.g. H2 specialization index) against the Vazquez null
Iobs <- networklevel(web = interactions, index = "H2") #0.7652146
nulls <- nullmodel(web=interactions, N=100, method='vaznull') # takes a while!
Inulls <- sapply(nulls, networklevel, index="H2")
plot(density(Inulls), xlim=c(0, 1), lwd=2, main="H'2")
abline(v=Iobs, col="red", lwd=2)

shapiro.test(Inulls) # p = 0.1527, not signif diff from normal distribution
res <- t.test(Inulls, mu = 0.7652146) #One sample t-test comparing nulls to observed mean
#t = -5520.4, df = 99, p-value < 2.2e-16

######Shannon diversity comparison
Shan.obs <- networklevel(web = interactions, index = "Shannon diversity") #3.054202
Shan.nulls <- sapply(nulls, networklevel, index="Shannon diversity")
plot(density(Inulls), xlim=c(0, 1), lwd=2, main="Shannon")
abline(v=Iobs, col="red", lwd=2)

shapiro.test(Shan.nulls) # p = 0.29
t.test(Shan.nulls, mu = 3.054202)
#t = 1439.8, df = 99, p-value < 2.2e-16

#compute modularity
mod <- computeModules(interactions)
#pdf(file = 'Adjacency_Matrix.pdf', width = 7, height = 7)
plotModuleWeb(mod, labsize = 0.5)
#dev.off()

#Compute Shannon Richness to summarize dietary niche width
Shannon<-estimate_richness(physeq = mergedspp, measures = c('Shannon', 'Observed'))
Shannon

#Add to metadata
mergedspp@sam_data$Shannon<-Shannon$Shannon

#PLot
#pdf('Shannon_Obs.pdf', width = 7, height = 5)
plot_richness(filtered, x="Binomen", measures=c("Observed", "Shannon")) + geom_boxplot()
#dev.off()

#Check for correlation with isotopic data from Oelbaum et al 2019
#write.csv(Shannon, file = "shannon_obs_alpha.csv")
Shan_SEAb<-read.csv(file = "shannon_obs_alpha.csv", header = T)
#Shan_SEAb$Niche<-as.factor(mergedspp@sam_data$Niche)

cor.test(Shan_SEAb$Shannon, Shan_SEAb$SEAb,  method = "pearson", use = "complete.obs")
#t = 0.65073, df = 17, p-value = 0.5239, r=0.15 
cor.test(Shan_SEAb$Observed, Shan_SEAb$SEAb,  method = "pearson", use = "complete.obs")
#t = 1.14, df = 17, p-value = 0.2701, r=0.27

# independent 2-group Mann-Whitney U Test
wilcox.test(Shan_SEAb$Shannon, Shan_SEAb$SEAb) # where y and x are numeric
#W = 81, p-value = 0.0001857

require(cowplot)
Sh<-ggplot(Shan_SEAb, aes(x=Shannon, y=SEAb)) +
  geom_point() + 
  geom_smooth(method=lm) +
  theme_bw()
Obs<-ggplot(Shan_SEAb, aes(x=Observed, y=SEAb)) +
  geom_point() + 
  geom_smooth(method=lm) +
  theme_bw()
#pdf(file = "Shan_SEAb_corr.pdf", width = 7, height = 5)
plot_grid(Sh, Obs, labels = "AUTO")
#dev.off()

#Some other basic data summaries
tax.df<-data.frame(filtered@tax_table)
length(which(tax.df$kingdom_name == "Plants"))
length(which(tax.df$kingdom_name == "Invertebrates"))

#############Split dataset into bins of top taxa to show shallower differentiation#############

######### COLEOPTERA #############
Beetles = subset_taxa(filtered, order_name=="Coleoptera")
Beetles = merge_samples(x = Beetles, group = "Binomen")

merged_beetles_family<-tax_glom(physeq = Beetles, taxrank = "family_name")
merged.beetle.motus<-as.data.frame(merged_beetles_family@otu_table)
merged.beetle.tax<-as.matrix((tax_table(merged_beetles_family)))
merged.beetle.tax<-data.frame(merged.beetle.tax[,1:3])
colnames(merged.beetle.motus)<-merged.beetle.tax[,3]

beetle_interactions<-as.matrix(merged.beetle.motus)

#pdf(file = 'Beetles_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(beetle_interactions, sort.order="inc"), text.rot = 90)
#dev.off()

######### HEMIPTERA #############
Bugs = subset_taxa(filtered, order_name=="Hemiptera")
Bugs = merge_samples(x = Bugs, group = "Binomen")

merged_bugs_family<-tax_glom(physeq = Bugs, taxrank = "family_name")
merged.bug.motus<-as.data.frame(merged_bugs_family@otu_table)
merged.bug.tax<-as.matrix((tax_table(merged_bugs_family)))
merged.bug.tax<-data.frame(merged.bug.tax[,1:3])
colnames(merged.bug.motus)<-merged.bug.tax[,3]

bug_interactions<-as.matrix(merged.bug.motus)

#pdf(file = 'Bugs_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(bug_interactions, sort.order="inc"), text.rot = 90)
#dev.off()

######### Lepidoptera #############
Moths = subset_taxa(filtered, order_name=="Lepidoptera")
Moths = merge_samples(x = Moths, group = "Binomen")

merged_moths_family<-tax_glom(physeq = Moths, taxrank = "family_name")
merged.moth.motus<-as.data.frame(merged_moths_family@otu_table)
merged.moth.tax<-as.matrix((tax_table(merged_moths_family)))
merged.moth.tax<-data.frame(merged.moth.tax[,1:3])
colnames(merged.moth.motus)<-merged.moth.tax[,3]

moth_interactions<-as.matrix(merged.moth.motus)

#pdf(file = 'Moths_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(moth_interactions, sort.order="inc"), text.rot = 90)
#dev.off()

######### Diptera #############
Flies = subset_taxa(filtered, order_name=="Diptera")
Flies = merge_samples(x = Flies, group = "Binomen")

merged_flies_family<-tax_glom(physeq = Flies, taxrank = "family_name")
merged.fly.motus<-as.data.frame(merged_flies_family@otu_table)
merged.fly.tax<-as.matrix((tax_table(merged_flies_family)))
merged.fly.tax<-data.frame(merged.fly.tax[,1:3])
colnames(merged.fly.motus)<-merged.fly.tax[,3]

fly_interactions<-as.matrix(merged.fly.motus)

#pdf(file = 'Flies_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(fly_interactions, sort.order="inc"), text.rot = 90)
#dev.off()

######### Plants #############
plant_eaters<-subset_samples(physeq = filtered, Niche=="Frugivore")
Plants = subset_taxa(plant_eaters, kingdom_name=="Plants")
Plants = merge_samples(x = Plants, group = "Binomen")

merged_plants_family<-tax_glom(physeq = Plants, taxrank = "family_name")
merged.pln.motus<-as.data.frame(merged_plants_family@otu_table)
merged.pln.tax<-as.matrix((tax_table(merged_plants_family)))
merged.pln.tax<-data.frame(merged.pln.tax[,1:3])
colnames(merged.pln.motus)<-merged.pln.tax[,3]

plants_interactions<-as.matrix(merged.pln.motus)

#pdf(file = 'Plants_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(plants_interactions, sort.order="inc"), text.rot = 90)
#dev.off()

###############Re-do Adjacency Plot for whole Community at Family Level############
######## CONVERT PHYLOSEQ MOTU TABLE AND TAX TABLE INTO BIPARTITE INTERACTION MATRIX 
merged_bats_Fam<-tax_glom(physeq = mergedspp, taxrank = "family_name")
merged.fam.motus<-as.data.frame(merged_bats_Fam@otu_table)
merged.fam.tax<-as.matrix((tax_table(merged_bats_Fam)))
merged.fam.tax<-data.frame(merged.fam.tax[,1:3])
colnames(merged.fam.motus)<-merged.fam.tax[,3]

family_interactions<-as.matrix(merged.fam.motus)

#Visualize web
#pdf(file = 'Family_Bipartite_Network.pdf', width = 16, height = 13)
plotweb(sortweb(family_interactions, sort.order="inc"), text.rot = 90)
dev.off()
visweb(sortweb(family_interactions, sort.order="inc"))

#compute modularity
fam.mod <- computeModules(family_interactions)
#pdf(file = 'Family_Adjacency_Matrix.pdf', width = 7, height = 7)
plotModuleWeb(fam.mod, labsize = 0.1)
#dev.off()

###############Cluster species by diet#############
#Create hierarchical clustering dendrogram 
require(pvclust)
#Create hierarchical clustering dendrogram
#Bootstrapped UPGMA using bray curtis distances
#Bray curtis dist
dist<-distance(merged_bats_Fam, "bray", type = "samples")
dist<-as.matrix(dist)
result <- pvclust(dist, method.dist="cor", method.hclust="average", nboot=1000)
plot(result, print.pv = T)

############### CORRELATION WITH ISOTOPIC SEAb #####################
#calculate Levin's niche breadth
#first flip columns and rows of interaction matrix
interactions.swap<-t(interactions)
levins<-niche.width(interactions.swap, method = "levins")
shannon<-niche.width(interactions.swap, method = "shannon")
levins<-t(levins)
levins<-as.data.frame(levins)
shannon<-t(shannon)
shannon<-as.data.frame(shannon)

Oelbaum19<- read.csv("~/Documents/Bat_Diets/Niche breadth-2.csv", row.names = 1)
SEAb<-as.matrix(as.array(Oelbaum19[,2]))
rownames(SEAb)<-rownames(Oelbaum19)
SEAb<-as.data.frame(SEAb)

merged_nichebreadth<-merge(SEAb,levins ["V1"],by="row.names",all.x=TRUE)

rownames(merged_nichebreadth)<-merged_nichebreadth$Row.names
merged_nichebreadth<-setNames(merged_nichebreadth, c("Species", "SEAb", "Levins"))
#remove species missing SEAb or Levins
library(tidyr)
merged_nichebreadth<-merged_nichebreadth %>% drop_na(Levins)

#Nonparametric correlation
cor.test(x = merged_nichebreadth$SEAb, y = merged_nichebreadth$Levins, method = "kendall")
#T = 75, p-value = 0.9405
