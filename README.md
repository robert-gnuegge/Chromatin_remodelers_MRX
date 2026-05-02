# Chromatin_remodelers_MRX [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19980318.svg)](https://doi.org/10.5281/zenodo.19980318)

Code to analyze the data of our project investigating the impact of chromatin remodelers on MRX nicking at DSBs.
See publication for software versions. Find R package versions below (output of `sessionInfo()`):

```
R version 4.5.3 (2026-03-11)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                  LC_TIME=en_US.UTF-8           LC_COLLATE=en_US.UTF-8        LC_MONETARY=en_US.UTF-8      
 [6] LC_MESSAGES=en_US.UTF-8       LC_PAPER=en_US.UTF-8          LC_NAME=en_US.UTF-8           LC_ADDRESS=en_US.UTF-8        LC_TELEPHONE=en_US.UTF-8     
[11] LC_MEASUREMENT=en_US.UTF-8    LC_IDENTIFICATION=en_US.UTF-8

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] ggseqlogo_0.2.2             ggplot2_4.0.3               rmelting_1.24.0             DNAshapeR_1.36.0            BSgenome_1.78.0             BiocIO_1.20.0              
 [7] plotrix_3.8-14              Gviz_1.54.0                 nucleR_2.42.0               GenomicAlignments_1.46.0    Rsamtools_2.26.0            SummarizedExperiment_1.40.0
[13] Biobase_2.70.0              MatrixGenerics_1.22.0       matrixStats_1.5.0           Biostrings_2.78.0           XVector_0.50.0              errors_0.4.4               
[19] propagate_1.1-0             rtracklayer_1.70.1          GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.1           
[25] BiocGenerics_0.56.0         generics_0.1.4             

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3       rstudioapi_0.18.0        jsonlite_2.0.0           magrittr_2.0.5           GenomicFeatures_1.62.0   farver_2.1.2            
  [7] rmarkdown_2.31           fields_17.1              vctrs_0.7.3              memoise_2.0.1            RCurl_1.98-1.18          base64enc_0.1-6         
 [13] htmltools_0.5.9          S4Arrays_1.10.1          progress_1.2.3           curl_7.1.0               SparseArray_1.10.10      Formula_1.2-5           
 [19] htmlwidgets_1.6.4        httr2_1.2.2              copula_1.1-7             cachem_1.1.0             lifecycle_1.0.5          minpack.lm_1.2-4        
 [25] pkgconfig_2.0.3          Matrix_1.7-5             R6_2.6.1                 fastmap_1.2.0            rbibutils_2.4.1          digest_0.6.39           
 [31] numDeriv_2016.8-1.1      colorspace_2.1-2         ShortRead_1.68.0         AnnotationDbi_1.72.0     Hmisc_5.2-5              RSQLite_2.4.6           
 [37] hwriter_1.3.2.1          filelock_1.0.3           pspline_1.0-21           httr_1.4.8               abind_1.4-8              compiler_4.5.3          
 [43] withr_3.0.2              bit64_4.8.0              gsl_2.1-9                htmlTable_2.5.0          S7_0.2.2                 backports_1.5.1         
 [49] BiocParallel_1.44.0      DBI_1.3.0                maps_3.4.3               biomaRt_2.66.2           rappdirs_0.3.4           DelayedArray_0.36.1     
 [55] rjson_0.2.23             tools_4.5.3              foreign_0.8-90           nnet_7.3-20              glue_1.8.1               bspm_0.5.7              
 [61] stabledist_0.7-2         restfulr_0.0.16          checkmate_2.3.4          cluster_2.1.8.2          hdf5r_1.3.12             gtable_0.3.6            
 [67] ensembldb_2.34.0         data.table_1.18.2.1      hms_1.1.4                pillar_1.11.1            stringr_1.6.0            spam_2.11-3             
 [73] rJava_1.0-18             dplyr_1.2.1              BiocFileCache_3.0.0      lattice_0.22-9           bit_4.6.0                deldir_2.0-4            
 [79] biovizBase_1.58.0        tidyselect_1.2.1         ADGofTest_0.3            knitr_1.51               gridExtra_2.3            ProtGenerics_1.42.0     
 [85] xfun_0.57                stringi_1.8.7            UCSC.utils_1.6.1         lazyeval_0.2.3           yaml_2.3.12              evaluate_1.0.5          
 [91] codetools_0.2-20         cigarillo_1.0.0          interp_1.1-6             tibble_3.3.1             cli_3.6.6                rpart_4.1.27            
 [97] Rdpack_2.6.6             dichromat_2.0-0.1        Rcpp_1.1.1-1.1           GenomeInfoDb_1.46.2      dbplyr_2.5.2             png_0.1-9               
[103] XML_3.99-0.23            blob_1.3.0               prettyunits_1.2.0        dotCall64_1.2            latticeExtra_0.6-31      jpeg_0.1-11             
[109] AnnotationFilter_1.34.0  bitops_1.0-9             pwalign_1.6.0            viridisLite_0.4.3        mvtnorm_1.3-7            VariantAnnotation_1.56.0
[115] scales_1.4.0             pcaPP_2.0-5              crayon_1.5.3             rlang_1.2.0              KEGGREST_1.50.0    
``` 

This project is licensed under the terms of the MIT license.
