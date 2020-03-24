# Prototype to process boxcar data



Reference: [BoxCar acquisition method enables single-shot proteomics
at a depth of 10,000 proteins i
100minutes](https://www.nature.com/articles/s41592-018-0003-5)




Load required packages and functions.


```r
library("MSnbase")
library("magrittr")
source("S01-bc-code.R")
```

Read the data in as an `MSnExp`


```r
f <- "./boxcar.mzML"
x <- readMSData(f, mode = "onDisk")
x
```

```
## MSn experiment data ("OnDiskMSnExp")
## Object size in memory: 6.04 Mb
## - - - Spectra data - - -
##  MS level(s): 1 
##  Number of spectra: 13572 
##  MSn retention times: 0:0 - 70:0 minutes
## - - - Processing information - - -
## Data loaded [Tue Mar 24 20:03:59 2020] 
##  MSnbase version: 2.13.3 
## - - - Meta data  - - -
## phenoData
##   rowNames: boxcar.mzML
##   varLabels: sampleNames
##   varMetadata: labelDescription
## Loaded from:
##   boxcar.mzML 
## protocolData: none
## featureData
##   featureNames: F1.S00001 F1.S00002 ... F1.S13572 (13572 total)
##   fvarLabels: fileIdx spIdx ... spectrum (35 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
```

Data processing (see `S01-bc-code.R` for details):

1. Define boxcar groups


```r
fData(x)$groups <- bc_groups(x)
```

2. Keep only boxcar spectra


```r
xbc <- x[bc_is_boxcar(x)] 
```

3. Combine boxcar spectra and coerce result back to an `MSnExp` object


```r
res <- xbc %>%
    bc_zero_out_box(offset = 0.5) %>%
    combineSpectra(fcol = "groups",
                   method = boxcarCombine) %>%
    .as_MSnExp(f)

res
```

```
## MSn experiment data ("MSnExp")
## Object size in memory: 457.2 Mb
## - - - Spectra data - - -
##  MS level(s): 1 
##  Number of spectra: 3393 
##  MSn retention times: 0:0 - 69:60 minutes
## - - - Processing information - - -
## Data converted from Spectra: Tue Mar 24 19:58:56 2020 
## boxcar processed [Tue Mar 24 19:59:00 2020] 
##  MSnbase version: 2.13.3 
## - - - Meta data  - - -
## phenoData
##   rowNames: 1
##   varLabels: sampleNames
##   varMetadata: labelDescription
## Loaded from:
##   boxcar.mzML 
## protocolData: none
## featureData
##   featureNames: F1.S00002 F1.S00006 ... F1.S13570 (3393 total)
##   fvarLabels: fileIdx spIdx ... groups (36 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
```

We have 3393 spectra after processing.

Write to a new mzML file


```r
writeMSData(res, "boxcar_processed.mzML")
```

See also `sp_plot`,  `bc_plot` and `bc_plotly` for visualisation. 

The figure below shows a slice of a combined spectrum with the
different boxcar spectra in different colours.

<iframe width="900" height="800" frameborder="0" scrolling="no" src="//plotly.com/~lgatto/4.embed"></iframe>



```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
## [1] magrittr_1.5        MSnbase_2.13.3      ProtGenerics_1.19.3
## [4] S4Vectors_0.24.3    mzR_2.20.0          Rcpp_1.0.4         
## [7] Biobase_2.46.0      BiocGenerics_0.32.0
## 
## loaded via a namespace (and not attached):
##  [1] BiocManager_1.30.10   plyr_1.8.6            compiler_3.6.2       
##  [4] pillar_1.4.3          iterators_1.0.12      zlibbioc_1.32.0      
##  [7] tools_3.6.2           digest_0.6.25         ncdf4_1.17           
## [10] MALDIquant_1.19.3     preprocessCore_1.48.0 evaluate_0.14        
## [13] lifecycle_0.2.0       tibble_2.1.3          gtable_0.3.0         
## [16] lattice_0.20-40       pkgconfig_2.0.3       rlang_0.4.5          
## [19] foreach_1.4.8         xfun_0.12             stringr_1.4.0        
## [22] dplyr_0.8.5           knitr_1.28            IRanges_2.20.2       
## [25] grid_3.6.2            tidyselect_1.0.0      glue_1.3.2           
## [28] impute_1.60.0         R6_2.4.1              XML_3.99-0.3         
## [31] BiocParallel_1.20.1   limma_3.42.2          ggplot2_3.3.0        
## [34] purrr_0.3.3           pcaMethods_1.78.0     scales_1.1.0         
## [37] codetools_0.2-16      MASS_7.3-51.5         mzID_1.24.0          
## [40] assertthat_0.2.1      colorspace_1.4-1      affy_1.64.0          
## [43] stringi_1.4.6         doParallel_1.0.15     munsell_0.5.0        
## [46] vsn_3.54.0            crayon_1.3.4          affyio_1.56.0
```
