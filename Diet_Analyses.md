

```r
#Download and load packages
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
#packageVersion("phyloseq")
library("ggplot2")
#packageVersion("ggplot2")
library(vegan)

setwd(dir = '/Users/mingala/Documents/Dissertation/Chapter2/Results 3.23.20/BatDiets3.24/')

#Import Feature table and metadata
motu_table <-read.delim(file = "FINAL_MOTUs.tsv", sep = "\t", header = T, row.names = 1)
taxonomy <- read.delim(file= "FINAL_TAXONOMY.tsv", sep = "\t", header = T, row.names = 1)
metadata<- read.delim(file ="Metadata_diet.txt", sep = "\t", header= T, row.names = 2)

#Create phyloseq object by converting each file into component phyloseq class objects
SAM<-sample_data(metadata)
TAX<-tax_table(as.matrix(taxonomy))
OTU<- otu_table(motu_table, taxa_are_rows=TRUE)
physeq<-merge_phyloseq(OTU, TAX, SAM)

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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

```r
#Detect contaminants by frequency
pruned<- subset_samples(physeq = physeq, DNA_CONC != "NA")
contamdf.freq <- isContaminant(pruned, method="frequency", conc="DNA_CONC")
```

```
## Warning in isContaminant(pruned, method = "frequency", conc = "DNA_CONC"): Removed 4
## samples with zero total counts (or frequency).
```

```r
head(contamdf.freq)
```

```
##                                                          freq prev    p.freq p.prev
## M03636:43:000000000-BG4JM:1:1103:22475:7113_CONS  0.077653596   26 0.2145112     NA
## M03636:43:000000000-BG4JM:1:1108:9536:15974_CONS  0.014170453    8 0.4242533     NA
## M03636:43:000000000-BG4JM:1:1101:10306:5691_CONS  0.080642517   20 0.5890707     NA
## M03636:43:000000000-BG4JM:1:1101:22443:3884_CONS  0.014501094    4 0.4704623     NA
## M03636:43:000000000-BG4JM:1:1103:22156:18630_CONS 0.009634287    2 0.5164774     NA
## M03636:43:000000000-BG4JM:1:1101:16195:2779_CONS  0.036906917    9 0.5582210     NA
##                                                           p contaminant
## M03636:43:000000000-BG4JM:1:1103:22475:7113_CONS  0.2145112       FALSE
## M03636:43:000000000-BG4JM:1:1108:9536:15974_CONS  0.4242533       FALSE
## M03636:43:000000000-BG4JM:1:1101:10306:5691_CONS  0.5890707       FALSE
## M03636:43:000000000-BG4JM:1:1101:22443:3884_CONS  0.4704623       FALSE
## M03636:43:000000000-BG4JM:1:1103:22156:18630_CONS 0.5164774       FALSE
## M03636:43:000000000-BG4JM:1:1101:16195:2779_CONS  0.5582210       FALSE
```

```r
table(contamdf.freq$contaminant)
```

```
## 
## FALSE  TRUE 
##   811    15
```

```r
head(which(contamdf.freq$contaminant))
```

```
## [1]  63 155 163 170 172 218
```

```r
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

#Transform to relative abundance
relative  = transform_sample_counts(filtered, function(OTU) OTU / sum(OTU))
plot_bar(physeq = relative, fill = "order_name", facet_grid =~kingdom_name, x = "Binomen")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png)

```r
#Top 25 taxa
#Barplot of top 25 OTUs
Top25OTUs = names(sort(taxa_sums(filtered), TRUE)[1:25])
comparetop25 = prune_taxa(Top25OTUs, filtered)
relative  = transform_sample_counts(comparetop25, function(OTU) OTU / sum(OTU))
plot_bar(physeq = relative, fill = "order_name", facet_grid =~kingdom_name, x = "Binomen")
```

```
## Warning: Removed 100 rows containing missing values (position_stack).
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-3.png)

```r
#Agglomerate by species
#Remove confidence from Tax table (messes up agglomeration)
tax_table(filtered) <- tax_table(filtered)[,1:5]
mergedspp<-merge_samples(x = filtered, group = "Binomen")
relative.merged  = transform_sample_counts(mergedspp, function(OTU) OTU / sum(OTU))

#p<-plot_bar(physeq = relative.merged, facet_grid=~kingdom_name)
#pdf(file = "Inv_Pln_Bars.pdf", width = 10, height = 7)
#p + geom_bar(aes(fill=order_name), stat="identity", position="stack")
#dev.off()

############## DATA TRANSFORMATION AND ORDINATION #################
#Center Log transform the data (don't rarefy)
require(microbiome)
physeq.trans<-microbiome::transform(filtered, "log10")
```

```
## Warning in microbiome::transform(filtered, "log10"): OTU table contains zeroes. Using
## log10(1 + x) transform.
```

```r
physeq.trans<-subset_samples(physeq = physeq.trans, Niche!="Control")

funx.ord <- ordinate(
  physeq = physeq.trans, 
  method = "NMDS", 
  distance = "bray"
)
```

```
## Run 0 stress 0.1613339 
## Run 1 stress 0.1644276 
## Run 2 stress 0.1615984 
## ... Procrustes: rmse 0.05901576  max resid 0.2495413 
## Run 3 stress 0.160883 
## ... New best solution
## ... Procrustes: rmse 0.03962871  max resid 0.1468315 
## Run 4 stress 0.1642467 
## Run 5 stress 0.1608314 
## ... New best solution
## ... Procrustes: rmse 0.007578199  max resid 0.04628627 
## Run 6 stress 0.1607892 
## ... New best solution
## ... Procrustes: rmse 0.0296131  max resid 0.1333369 
## Run 7 stress 0.1605482 
## ... New best solution
## ... Procrustes: rmse 0.03196537  max resid 0.1298247 
## Run 8 stress 0.1653758 
## Run 9 stress 0.162219 
## Run 10 stress 0.1618395 
## Run 11 stress 0.1619474 
## Run 12 stress 0.1630624 
## Run 13 stress 0.1642081 
## Run 14 stress 0.1655518 
## Run 15 stress 0.1644304 
## Run 16 stress 0.1643658 
## Run 17 stress 0.1604105 
## ... New best solution
## ... Procrustes: rmse 0.02635508  max resid 0.1296647 
## Run 18 stress 0.1648805 
## Run 19 stress 0.162455 
## Run 20 stress 0.1617708 
## *** No convergence -- monoMDS stopping criteria:
##      6: no. of iterations >= maxit
##     14: stress ratio > sratmax
```

```r
# Calculate distance matrix
study.bray <- phyloseq::distance(physeq.trans, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.trans))

require(viridis)
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
  stat_ellipse() + geom_jitter() + 
  geom_text(label=physeq.trans@sam_data$Binomen,
    nudge_x = 0.25, nudge_y = 0.25,size = 2,
    check_overlap = T)
```

```
## Too few points to calculate an ellipse
## Too few points to calculate an ellipse
```

```
## Warning: Removed 2 rows containing missing values (geom_path).
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-4.png)

```r
############# PERMANOVA OF DIET #############
sampledf$Niche<-as.factor(sampledf$Niche)
adonis(study.bray ~ Niche * Binomen + Family, data = sampledf)
```

```
## 
## Call:
## adonis(formula = study.bray ~ Niche * Binomen + Family, data = sampledf) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## Niche      3    2.6846 0.89487  2.3708 0.09014  0.001 ***
## Binomen   20   11.6227 0.58114  1.5396 0.39025  0.001 ***
## Residuals 41   15.4757 0.37746         0.51962           
## Total     64   29.7830                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Niche      3    2.6846 0.89487  2.3708 0.09014  0.001 ***
#Binomen   20   11.6227 0.58114  1.5396 0.39025  0.001 ***
#Residuals 41   15.4757 0.37746         0.51962           
#Total     64   29.7830                 1.00000          

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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-5.png)

```r
#dev.off()
visweb(sortweb(interactions, sort.order="inc"))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-6.png)

```r
#compute modularity
mod <- computeModules(interactions)
#pdf(file = 'Adjacency_Matrix.pdf', width = 7, height = 7)
plotModuleWeb(mod, labsize = 0.5)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-7.png)

```r
#dev.off()

#Compute Shannon Richness to summarize dietary niche width
Shannon<-estimate_richness(physeq = mergedspp, measures = c('Shannon', 'Observed'))
Shannon
```

```
##                            Observed   Shannon
## Artibeus_intermedius              9 1.7764022
## Artibeus_jamaicensis             43 0.6984114
## Artibeus_lituratus               41 0.9807187
## Bauerus_dubiaquercus             39 2.6697901
## Carollia_perspicillata           30 0.8038310
## Carollia_sowelli                 68 1.1337826
## Chrotopterus_auritus             38 0.8218339
## Dermanura_phaeotis               45 0.3254159
## Dermanura_watsoni                37 0.8962527
## Eptesicus_furinalis              47 0.8652175
## Gardnernycteris_crenulatum        4 1.0626812
## Glossophaga_soricina             25 1.4087635
## Lophostoma_evotis                10 1.7231883
## Mimon_cozumelae                  32 0.4773561
## Molossus_nigricans              363 3.2279977
## Mormoops_megalophylla            29 1.7012887
## Myotis_elegans                  144 2.8399350
## Myotis_keasyi                   148 2.7840347
## Natalus_mexicanus                67 1.9685395
## Pteronotus_mesoamericanus        45 1.6472794
## Rhogeessa_aenaeus                45 0.3812481
## Rhynchonycteris_naso             69 1.1321234
## Sturnira_parvidens               38 0.3675559
## Trachops_cirrhosus               15 0.4210451
```

```r
#Add to metadata
mergedspp@sam_data$Shannon<-Shannon$Shannon

#PLot
#pdf('Shannon_Obs.pdf', width = 7, height = 5)
plot_richness(filtered, x="Binomen", measures=c("Observed", "Shannon")) + geom_boxplot()
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-8.png)

```r
#dev.off()

#Check for correlation with isotopic data from Oelbaum et al 2019
#write.csv(Shannon, file = "shannon_obs_alpha.csv")
Shan_SEAb<-read.csv(file = "shannon_obs_alpha.csv", header = T)
#Shan_SEAb$Niche<-as.factor(mergedspp@sam_data$Niche)

cor.test(Shan_SEAb$Shannon, Shan_SEAb$SEAb,  method = "pearson", use = "complete.obs")
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Shan_SEAb$Shannon and Shan_SEAb$SEAb
## t = 0.65073, df = 17, p-value = 0.5239
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.3210467  0.5697606
## sample estimates:
##       cor 
## 0.1558956
```

```r
#t = 0.65073, df = 17, p-value = 0.5239, r=0.15 
cor.test(Shan_SEAb$Observed, Shan_SEAb$SEAb,  method = "pearson", use = "complete.obs")
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Shan_SEAb$Observed and Shan_SEAb$SEAb
## t = 1.14, df = 17, p-value = 0.2701
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.2135771  0.6428791
## sample estimates:
##       cor 
## 0.2664834
```

```r
#t = 1.14, df = 17, p-value = 0.2701, r=0.27

# independent 2-group Mann-Whitney U Test
wilcox.test(Shan_SEAb$Shannon, Shan_SEAb$SEAb) # where y and x are numeric
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  Shan_SEAb$Shannon and Shan_SEAb$SEAb
## W = 81, p-value = 0.0001857
## alternative hypothesis: true location shift is not equal to 0
```

```r
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
```

```
## Warning: Removed 5 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

```
## Warning: Removed 5 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-9.png)

```r
#dev.off()

#Some other basic data summaries
tax.df<-data.frame(filtered@tax_table)
length(which(tax.df$kingdom_name == "Plants"))
```

```
## [1] 194
```

```r
length(which(tax.df$kingdom_name == "Invertebrates"))
```

```
## [1] 617
```

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-10.png)

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-11.png)

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-12.png)

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-13.png)

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-14.png)

```r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-15.png)

```r
dev.off()
```

```
## RStudioGD 
##         2
```

```r
visweb(sortweb(family_interactions, sort.order="inc"))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-16.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-17.png)

```r
#compute modularity
fam.mod <- computeModules(family_interactions)
#pdf(file = 'Family_Adjacency_Matrix.pdf', width = 7, height = 7)
plotModuleWeb(fam.mod, labsize = 0.1)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-18.png)

```r
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
```

```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.58)... Done.
## Bootstrap (r = 0.67)... Done.
## Bootstrap (r = 0.79)... Done.
## Bootstrap (r = 0.88)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.08)... Done.
## Bootstrap (r = 1.17)... Done.
## Bootstrap (r = 1.29)... Done.
## Bootstrap (r = 1.38)... Done.
```

```r
plot(result, print.pv = T)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-19.png)

