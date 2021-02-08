Libraries
=========

``` r
library(bigstatsr)
library(data.table)
library(bigsnpr)
library(tidyverse) 
```

    ## ── Attaching packages ────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.1     ✓ dplyr   1.0.2
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ───────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::between()   masks data.table::between()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::first()     masks data.table::first()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::last()      masks data.table::last()
    ## x purrr::transpose() masks data.table::transpose()

``` r
theme_set(theme_bw())
```

Load genetic data
=================

`.bed` files were created by using plink with the following parameters:
\* MAF: 0.03 \* geno: 0.1

Convert plink files to `.rds` object
------------------------------------

This code will create a `.rds` object from the plink files and only
needs to be run once. An error will be thrown if the given directory
containing (`sample_dir` in this example) does not contain other
associated plink files e.g. `.fam` as well.

``` r
dir<-"sample_dir/sample.bed"
bed<-snp_readBed(dir) #creates an .rds file- neceassry for the following steps. 
```

Attach and impute genetic data
------------------------------

`.rds` file from the previous step will be attached for further
analyses.

``` r
dir<-"sample_dir/sample.rds" #rds file generated from the previous code chunk
obj.bigSNPall <- snp_attach(dir) 
G<-obj.bigSNPall$genotypes
G.imp<-snp_fastImputeSimple(G)
y<-obj.bigSNPall$fam$affection-1
big_counts(G.imp,ind.col = 1:10)
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## 0    4455 5197 5414 5386 6921 5623 6411 4394 4324  5463
    ## 1    2805 2298 2291 2244 1133 2431 1523 3055 3140  2348
    ## 2     794  559  349  424    0    0  120  605  590   243
    ## <NA>    0    0    0    0    0    0    0    0    0     0

Load Models
===========

SCT-PS
------

``` r
# load sct-ps betas
all_snps<-fread("Mariam_et_al_ACCORD_C4_SCT-PS.csv",data.table = F,
                skip=1)
sumstats<-all_snps[,c('chr', 'rsid', 'pos', 'a0', 'a1', 'beta', 'p')]
sumstats<-sumstats[which(sumstats$beta!=0),]
map <- obj.bigSNPall$map[,-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map)
```

    ## 178,674 variants to be matched.

    ## 0 ambiguous SNPs have been removed.

    ## 178,136 variants have been matched; 0 were flipped and 438 were reversed.

``` r
x<-obj.bigSNPall$map$marker.ID[which(obj.bigSNPall$map$marker.ID %in% info_snp$rsid)]
info_sct_ps<-info_snp[order(match(info_snp$rsid,x)),]
# load sct-ps model
final_model<-read_rds("Mariam_et_al_ACCORD_C4_SCT-PS.rds")
summary(final_model$mod)
```

    ## # A tibble: 3 x 9
    ##    alpha power_adaptive power_scale validation_loss intercept beta  nb_var
    ##    <dbl>          <dbl>       <dbl>           <dbl>     <dbl> <lis>  <int>
    ## 1 0.0001              0           1           0.317     -115. <dbl…  26362
    ## 2 0.01                0           1           0.308     -130. <dbl…   6282
    ## 3 1                   0           1           0.309     -136. <dbl…   2241
    ## # … with 2 more variables: message <list>, all_conv <lgl>

CT-PS Model
-----------

``` r
# load ct-ps betas
all_snps<-fread("Mariam_et_al_ACCORD_C4_CT-PS.csv",data.table = F,
                skip = 1)
sumstats<-all_snps[,c('chr', 'rsid', 'pos', 'a0', 'a1', 'beta', 'p')]
sumstats$beta<-as.numeric(sumstats$beta)
map <- obj.bigSNPall$map[,-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map)
```

    ## 72,925 variants to be matched.

    ## 0 ambiguous SNPs have been removed.

    ## 72,632 variants have been matched; 0 were flipped and 115 were reversed.

``` r
x<-obj.bigSNPall$map$marker.ID[which(obj.bigSNPall$map$marker.ID %in% info_snp$rsid)]
info_ct_ps<-info_snp[order(match(info_snp$rsid,x)),]
```

Apply PS
========

SCT-PS
------

1.  SNP data is subset down to contain only SNPs involved in SCT-PS.
2.  SCT-PS scores are calculated.
3.  A pre-specified threshold is used to call individuals predicted to
    benefit from intensive treatment.

``` r
G.imp.sub<-big_copy(G.imp,
                    ind.row = 1:nrow(G.imp),
                    ind.col=which(obj.bigSNPall$map$marker.ID %in% info_sct_ps$rsid))
CHR.sub <- obj.bigSNPall$map$chromosome[which(obj.bigSNPall$map$marker.ID %in% info_sct_ps$rsid)]
POS.sub <- obj.bigSNPall$map$physical.pos[which(obj.bigSNPall$map$marker.ID %in% info_sct_ps$rsid)]
beta <- (info_sct_ps$beta)
lpval <- -log10(info_sct_ps$p)
threshold<- -0.6060066
pred<-final_model$intercept+
  big_prodVec(G.imp.sub,beta)
preds.sct_ps <- data.frame(cbind(obj.bigSNPall$fam$family.ID,pred,
                            ifelse(pred>threshold, 1, 0)))
colnames(preds.sct_ps)<-c("ID","sctPS","prediction")
table(preds.sct_ps[,3])
```

    ## 
    ##    0    1 
    ## 4220 3834

CT-PS
-----

1.  SNP data is subset down to contain only SNPs involved in CT-PS.
2.  CT-PS scores are calculated.
3.  A pre-specified threshold is used to call individuals predicted to
    benefit from intensive treatment.

``` r
G.imp.sub<-big_copy(G.imp,
                    ind.row = 1:nrow(G.imp),
                    ind.col=which(obj.bigSNPall$map$marker.ID %in% info_ct_ps$rsid))
CHR.sub <- obj.bigSNPall$map$chromosome[which(obj.bigSNPall$map$marker.ID %in% info_ct_ps$rsid)]
POS.sub <- obj.bigSNPall$map$physical.pos[which(obj.bigSNPall$map$marker.ID %in% info_ct_ps$rsid)]
beta <- (info_ct_ps$beta)
lpval <- -log10(info_ct_ps$p)
threshold<- 5463.685
ct_ps<-snp_PRS(G.imp.sub, beta,
        lpS.keep = lpval, thr.list = 1.073004)
preds.ct_ps <- data.frame(cbind(obj.bigSNPall$fam$family.ID,ct_ps,
                            ifelse(ct_ps>threshold, 1, 0)))
colnames(preds.ct_ps)<-c("ID","ctPS","prediction")
table(preds.ct_ps[,3])
```

    ## 
    ##    0    1 
    ## 3336 4718

Save Output
===========

``` r
write.csv(preds.sct_ps,"sct_ps_preds.csv",row.names = F)
write.csv(preds.ct_ps,"ct_ps_preds.csv",row.names = F)
```

SessionInfo
===========

``` r
sessionInfo()
```

    ## R version 4.0.0 (2020-04-24)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib64/libblas.so.3.4.2
    ## LAPACK: /opt/R-4.0.0/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] forcats_0.5.0     stringr_1.4.0     dplyr_1.0.2       purrr_0.3.4      
    ##  [5] readr_1.3.1       tidyr_1.1.0       tibble_3.0.1      ggplot2_3.3.0    
    ##  [9] tidyverse_1.3.0   bigsnpr_1.5.2     data.table_1.13.6 bigstatsr_1.3.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.5         lubridate_1.7.8    lattice_0.20-41    utf8_1.1.4        
    ##  [5] assertthat_0.2.1   digest_0.6.25      foreach_1.5.0      R6_2.4.1          
    ##  [9] bigsparser_0.4.0   cellranger_1.1.0   backports_1.1.6    reprex_0.3.0      
    ## [13] evaluate_0.14      httr_1.4.2         pillar_1.4.3       flock_0.7         
    ## [17] rlang_0.4.7        readxl_1.3.1       rstudioapi_0.11    Matrix_1.2-18     
    ## [21] rmarkdown_2.1      bigparallelr_0.3.0 munsell_0.5.0      broom_0.5.6       
    ## [25] compiler_4.0.0     modelr_0.1.6       xfun_0.13          pkgconfig_2.0.3   
    ## [29] htmltools_0.4.0    tidyselect_1.1.0   codetools_0.2-16   fansi_0.4.1       
    ## [33] crayon_1.3.4       dbplyr_1.4.3       withr_2.2.0        grid_4.0.0        
    ## [37] nlme_3.1-147       jsonlite_1.7.1     gtable_0.3.0       lifecycle_0.2.0   
    ## [41] DBI_1.1.0          magrittr_1.5       scales_1.1.0       cli_2.0.2         
    ## [45] stringi_1.4.6      fs_1.4.1           doParallel_1.0.15  xml2_1.3.2        
    ## [49] ellipsis_0.3.1     generics_0.1.0     vctrs_0.3.2        cowplot_1.1.0     
    ## [53] iterators_1.0.12   tools_4.0.0        glue_1.4.0         hms_0.5.3         
    ## [57] parallel_4.0.0     yaml_2.2.1         colorspace_1.4-1   bigassertr_0.1.3  
    ## [61] rvest_0.3.5        knitr_1.28         haven_2.2.0
