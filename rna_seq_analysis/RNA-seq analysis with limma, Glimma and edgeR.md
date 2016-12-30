

# RNA-seq analysis with limma, Glimma and edgeR

## Overview

This is a pipeline that is published Nov 30 2016. Here is the quote to describe this pipeline.
'''In this workflow article, we analyse RNA-sequencing data from the mouse mammary gland, demonstrating use of the popular edgeR package to import, organise, filter and normalise the data, followed by the limma package with its voom method, linear modelling and empirical Bayes moderation to assess differential expression and perform gene set testing. This pipeline is further enhanced by the Glimma package which enables interactive exploration of the results so that individual samples and genes can be examined by the user. The complete analysis offered by these three packages highlights the ease with which researchers can turn the raw counts from an RNA-sequencing experiment into biological insights using Bioconductor.'''

### Set-up
~~~
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
~~~

### Data packaing
~~~
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)
~~~
edgeR has readDGE to directly read all different files into one matrix, a DGEList object
~~~
x <- readDGE(files, columns=c(1,3))
class(x)
~~~

#### Organizing sample information

~~~R
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane

~~~

#### Organizeng gene annotations

~~~
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)
~~~

check duplicates

~~~
genes <- genes[!duplicated(genes$ENTREZID),]
~~~

#### Data pre-processing

Transformatons from the raw-scale

popular transformations include counts per million (CPM), log2-counts per million (log-CPM), etc.

~~~
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
~~~

removing genes that are lowly expressed

~~~
table(rowSums(x$counts==0)==9)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
~~~

visualizing the filtered data

~~~
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
 den <- density(lcpm[,i])
 lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
   den <- density(lcpm[,i])
   lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
~~~

#### Normalizing gene expression distributions

~~~R
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
~~~

visualizing the normalization

~~~
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

## [1] 0.0547 6.1306 1.2293 1.1705 1.2149 1.0562 1.1459 1.2613 1.1170

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
~~~

#### Unsupervised clustering of samples

~~~
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
~~~

## Differential expresson analysis

#### Creating a design matrix and contrasts

~~~
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP, 
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML, 
   levels = colnames(design))
contr.matrix
~~~

#### Removing heterosedascity from count data

~~~
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
~~~

#### Examining the number of DE genes

~~~
summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

## [1] 2409

head(tfit$genes$SYMBOL[de.common], n=20)

##  [1] "Xkr4"          "Rgs20"         "Cpa6"          "Sulf1"         "Eya1"         
##  [6] "Msc"           "Sbspon"        "Pi15"          "Crispld1"      "Kcnq5"        
## [11] "Ptpn18"        "Arhgef4"       "2010300C02Rik" "Aff3"          "Npas2"        
## [16] "Tbc1d8"        "Creg2"         "Il1r1"         "Il18r1"        "Il18rap"

vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
~~~

#### Examining individual DE genes from top to bottom

~~~
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)

head(basal.vs.ml)
~~~

#### Useful graphical representations of differential expression results

~~~
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
~~~

~~~
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         id.column="ENTREZID", counts=x$counts, groups=group, launch=FALSE)
~~~

~~~
# clustered heatmap
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
   labRow=v$genes$SYMBOL[i], labCol=group, 
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")
~~~

#### Gene set testing with camera

~~~
load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)

##                                             NGenes Direction   PValue      FDR
## LIM_MAMMARY_STEM_CELL_UP                       739        Up 1.13e-18 5.36e-15
## LIM_MAMMARY_STEM_CELL_DN                       630      Down 1.57e-15 3.71e-12
## ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER    163        Up 1.44e-13 2.26e-10
## SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP         183        Up 2.18e-13 2.58e-10
## LIM_MAMMARY_LUMINAL_PROGENITOR_UP               87      Down 6.73e-13 6.36e-10

cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)

##                                             NGenes Direction   PValue      FDR
## LIM_MAMMARY_STEM_CELL_UP                       739        Up 5.09e-23 2.40e-19
## LIM_MAMMARY_STEM_CELL_DN                       630      Down 5.13e-19 1.21e-15
## LIM_MAMMARY_LUMINAL_MATURE_DN                  166        Up 8.88e-16 1.40e-12
## LIM_MAMMARY_LUMINAL_MATURE_UP                  180      Down 6.29e-13 7.43e-10
## ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER    163        Up 1.68e-12 1.59e-09

cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)

##                                         NGenes Direction   PValue      FDR
## LIM_MAMMARY_LUMINAL_MATURE_UP              180      Down 8.50e-14 3.40e-10
## LIM_MAMMARY_LUMINAL_MATURE_DN              166        Up 1.44e-13 3.40e-10
## LIM_MAMMARY_LUMINAL_PROGENITOR_UP           87        Up 3.84e-11 6.05e-08
## REACTOME_RESPIRATORY_ELECTRON_TRANSPORT     91      Down 2.66e-08 3.14e-05
## NABA_CORE_MATRISOME                        222        Up 4.43e-08 4.19e-05
~~~

barcodeplot

~~~
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
~~~

## SessionInfo()

~~~
## R version 3.3.1 (2016-06-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
## [9] base     
## 
## other attached packages:
##  [1] gplots_3.0.1                             RColorBrewer_1.1-2                      
##  [3] Mus.musculus_1.3.1                       TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.0
##  [5] org.Mm.eg.db_3.4.0                       GO.db_3.4.0                             
##  [7] OrganismDbi_1.16.0                       GenomicFeatures_1.26.0                  
##  [9] GenomicRanges_1.26.1                     GenomeInfoDb_1.10.0                     
## [11] AnnotationDbi_1.36.0                     IRanges_2.8.0                           
## [13] S4Vectors_0.12.0                         Biobase_2.34.0                          
## [15] BiocGenerics_0.20.0                      edgeR_3.16.0                            
## [17] Glimma_1.2.0                             limma_3.30.0                            
## [19] BiocStyle_2.3.1                          rmarkdown_1.1                           
## [21] knitr_1.14                              
## 
## loaded via a namespace (and not attached):
##  [1] splines_3.3.1              R.utils_2.4.0              gtools_3.5.0              
##  [4] Formula_1.2-1              assertthat_0.1             latticeExtra_0.6-28       
##  [7] RBGL_1.50.0                Rsamtools_1.26.0           yaml_2.1.13               
## [10] RSQLite_1.0.0              lattice_0.20-34            chron_2.3-47              
## [13] digest_0.6.10              XVector_0.14.0             colorspace_1.2-7          
## [16] htmltools_0.3.5            Matrix_1.2-6               R.oo_1.20.0               
## [19] plyr_1.8.4                 DESeq2_1.14.0              XML_3.98-1.4              
## [22] biomaRt_2.30.0             genefilter_1.56.0          zlibbioc_1.20.0           
## [25] xtable_1.8-2               scales_0.4.0               gdata_2.17.0              
## [28] BiocParallel_1.8.0         tibble_1.2                 annotate_1.52.0           
## [31] ggplot2_2.1.0              SummarizedExperiment_1.4.0 nnet_7.3-12               
## [34] survival_2.39-5            magrittr_1.5               evaluate_0.10             
## [37] R.methodsS3_1.7.1          foreign_0.8-66             graph_1.52.0              
## [40] BiocInstaller_1.24.0       tools_3.3.1                data.table_1.9.6          
## [43] formatR_1.4                stringr_1.1.0              munsell_0.4.3             
## [46] locfit_1.5-9.1             cluster_2.0.4              Biostrings_2.42.0         
## [49] caTools_1.17.1             grid_3.3.1                 RCurl_1.95-4.8            
## [52] bitops_1.0-6               gtable_0.2.0               DBI_0.5-1                 
## [55] GenomicAlignments_1.10.0   gridExtra_2.2.1            rtracklayer_1.34.0        
## [58] Hmisc_3.17-4               KernSmooth_2.23-15         stringi_1.1.2             
## [61] Rcpp_0.12.7                geneplotter_1.52.0         rpart_4.1-10              
## [64] acepack_1.4.0
~~~



















































































## sessionInfo()

