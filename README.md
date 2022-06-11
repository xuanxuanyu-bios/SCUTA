# scMLLM
 linear mixed modelling for single cell RNAseq data with multilevel longitudinal design
 
### install scMLLM
```
library(devtools)
install_github("xuanxuanyu-bios/scMLLM")
```
### Analysis flowchart

![pipelineflowchart](https://user-images.githubusercontent.com/66747045/169596929-e20da322-8ac5-4cdc-8e12-12fc04be6f8e.png)


### Model fitting tutorial
`scMLLM` is a tools that designed for fitting linear mixed models for single cell RNAseq datasets, especially with longitudinal multi-level design. The algorithm is based on `Dream` in `VariancePartition` package.
```
Load library and data
library("scMLLM")
```
The imput matrix is the count matrix and meta data. 
counts is expression matrix whhere columns represent cells, rows represent genes.
coldata is the meta data including condition, individual and time information of each cell. 
```
data(example.data)
```
Load libraries
```
library("zinbwave")
library("DESeq2")
library("edgeR")
library("variancePartition")
```
prepare DESeqDataSet file and specify the variables that are taken into acount handling dropout events. 
```
counts  <- example.data$counts
coldata <- example.data$coldata

fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition + time)
```
The zinbwave function is used to compute observational weights which unlock bulk RNA-seq tools for single-cell applications, as illustrated in [Van den Berge et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1406-4).
```
zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000, observationalWeights = TRUE)
weights <- assay(fluidigm_zinb, "weights")
```
`voomWithDreamWeights` function is used to transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and compute observation-level weights. Then, the zinb weights and mean-variance weights are combined together as the overall weights. At last, a three-level linear mixed model is fitted for the data.
```
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
form <- ~ condition + time + (1 + time |individual)
vobjDream = voomWithDreamWeights( d0, form, coldata )
vobjDream.weight<-vobjDream
vobjDream.weight$weights <- weights*vobjDream.weight$weights
```
`getContrast` is used to specify the contrast matrix for linear mixed model. The three-level linear mixed model is  fitted by suing scMLLM function.
```
fit.scMLLM     <- scMLLM( vobjDream.weight, form, coldata)
```

The `scMLLM()` function is modified to replace the 'dream()' function in variancePartition, so that any  function in variance partition that used combined with `dream()` function can be used in conjuction with 'scMLLM()' function. For example, the top 6 differentiall expressed genes between two conditions are 

```
fit.scMLLM.res <- topTable(fit.scMLLM, coef="condition2", number=nrow(counts) )
head(fit.scMLLM.res)

            logFC  AveExpr         t      P.Value    adj.P.Val     z.std
Gene558 -1.451388 4.991292 -9.796536 1.481199e-09 1.478237e-06 -6.046398
Gene569 -1.479936 4.852454 -8.988242 6.702487e-08 3.195357e-05 -5.398969
Gene740 -1.337853 5.140980 -8.551072 9.605281e-08 3.195357e-05 -5.334037
Gene379 -1.452054 5.376776 -7.665623 4.786163e-07 9.380289e-05 -5.034693
Gene98  -1.358545 5.261751 -7.681800 5.068412e-07 9.380289e-05 -5.023705
Gene213 -1.286129 5.877459 -7.457352 5.639453e-07 9.380289e-05 -5.003172
```
Where logFC is the log fold change comparing condition 2 with condition 1, P.Value and adj.P.Val are the unadjusted and FDR adjusted P values, respectively.
