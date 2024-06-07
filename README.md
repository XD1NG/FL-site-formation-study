# Deciphering site formation processes in lakeshore environments: a case study of Early Pleistocene site of Feiliang (Nihewan Basin, China): data and scripts

This repository contains the data and R scripts used to explore the site formation processes at FL in the following paper:

Ding, X., Pei, S., Ma, D. et al. Deciphering site formation processes in lakeshore environments: a case study of early Pleistocene site of Feiliang (Nihewan Basin, China). Archaeol Anthropol Sci 16, 92 (2024). https://doi.org/10.1007/s12520-024-02003-7.

Three parts of analysis are included in this study, including taphonomic analysis, spatial analysis and fabric analysis.

## Data Sets and Data Preparation

All raw data used in the paper, including four CSV files and four shapefiles, can be found under the FL_data and FL_outline folder in the repository.

The file FL_basic data.csv contains the basic information of archaeological finds at FL, including find IDs, trenches, material types, coordinates, altitudes, orientations, inclinations, dimensions, weights, technological category, rock types and degrees of weathering. It is originally published as supplementary materials in the following paper:

Pei S, Xie F, Deng C, Jia Z, Wang X, Guan Y, et al. (2017) Early Pleistocene archaeological occurrences at the Feiliang site, and the archaeology of human origins in the Nihewan Basin, North China. PLoS ONE 12(11): e0187251. DOP: <https://doi.org/10.1371/journal.pone.0187251>.

The file FL_fabric index.csv is the results of fabric strength, flatness index, cluster-girdle index, K index, isotropy index and elongation index. They were calculated with three normalised eigenvectors (S1, S2 and S3) which were produced with the long-axis orientation and degrees of objects in the software Stereonet. They will be used in the three-dimensional fabric analysis.

The folder FL_outline contains shapefiles of trench outlines of FL, digitalising the original excavation maps with ArcGIS software. They will be used in the spatial analysis as observation windows.

Two files, FL_Refit_T1_weight.csv and FL_Refit_TOK_weight.csv, contain information on refit sets found at FL. The IDs, coordinates, altitudes, and weights of refitted pieces, as well as distances, orientations, degrees, azimuth and dips of refit lines are presented. They will be used as part of taphonomic analysis.

## Taphonomic analysis

Edge rounding, dimensions and refit sets of lithic artefacts were analysed in this part. Bar figures were plotted with ggplot2 to show the composition of rounding groups at each trench and within subgroups (length groups, weight groups and technological groups). Function drawrefits3viewweight.fl() was created to plot the three-dimensional maps of refit lines. The azimuths of refit lines were also plotted as rose diagrams and tested with Rayleigh and Kuiper test.

## Spatial analysis

Kernel density estimate, local autocorrelation (Local G\*), chi-square test and correlation analysis were applied in this part. Several functions were created here to conduct calculations and plot figures:

Function kdensitymap.xy.75contour4fossils() could plot the kernel density surface of all materials with black circles showing the cluster centre (kernel density higher than 75th percentile of all) of fossils.

Function kdensitymap.xy() would plot the kernel density surfaces of data input, which in this research lithic artefacts.

localGplot() could calculate and plot the local autocorrelations of input data, thus the rounding, length and weight of lithic artefacts were considered here.

chisquaretest.2half.fl() aimed to automatically separate the trench vertically (upper and lower), transversally (west and east) and horizontally (north and south) from the middle into two halves, and the frequency of fossils and lithic artefacts was summed for each part and then compared by chi-square test.

## Fabric analysis

In this part, two-dimensional long-axis orientation and three-dimensional orientation and degrees were analysed. Function CurrayL() was created to calculate the Curray' L index. 2D long-axis orientations were plotted as rose diagrams and tested with Rayleigh and Kuiper test, while results of strength C, flatness index, cluster-girdle index and K index were plotted in bar figure with ggplot2.

## R Session Info

```         
R version 4.3.3
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.5.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/London

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
 [1] renv_1.0.5             ggridges_0.5.4        
 [3] ggbiplot_0.6.2         ggplot2_3.5.0         
 [5] spdep_1.2-8            sf_1.0-14             
 [7] spData_2.3.0           Hmisc_5.1-0           
 [9] corrplot_0.92          circular_0.4-95       
[11] raster_3.6-23          rgeos_0.6-4           
[13] rgdal_1.6-7            colorRamps_2.3.1      
[15] spatstat_3.0-6         spatstat.linnet_3.1-1 
[17] spatstat.model_3.2-4   rpart_4.1.23          
[19] spatstat.explore_3.2-1 nlme_3.1-164          
[21] spatstat.random_3.1-5  spatstat.geom_3.2-4   
[23] spatstat.data_3.0-1    maptools_1.1-8        
[25] sp_2.0-0              

loaded via a namespace (and not attached):
 [1] farver_2.1.1          fastmap_1.1.1         digest_0.6.33      
 [4] lifecycle_1.0.3       cluster_2.1.6         terra_1.7-39       
 [7] magrittr_2.0.3        compiler_4.3.3        rlang_1.1.1        
[10] tools_4.3.3           utf8_1.2.3            data.table_1.14.8  
[13] knitr_1.43            labeling_0.4.2        htmlwidgets_1.6.2  
[16] classInt_0.4-9        abind_1.4-5           KernSmooth_2.23-22 
[19] withr_2.5.0           foreign_0.8-86        nnet_7.3-19        
[22] grid_4.3.3            polyclip_1.10-4       fansi_1.0.4        
[25] e1071_1.7-13          colorspace_2.1-0      scales_1.3.0       
[28] spatstat.utils_3.0-3  cli_3.6.1             mvtnorm_1.2-2      
[31] crayon_1.5.2          rmarkdown_2.24        rstudioapi_0.15.0  
[34] DBI_1.1.3             proxy_0.4-27          stringr_1.5.0      
[37] splines_4.3.3         s2_1.1.4              base64enc_0.1-3    
[40] vctrs_0.6.3           boot_1.3-29           Matrix_1.6-5       
[43] Formula_1.2-5         htmlTable_2.4.1       tensor_1.5         
[46] units_0.8-3           goftest_1.2-3         glue_1.6.2         
[49] codetools_0.2-19      stringi_1.7.12        gtable_0.3.3       
[52] deldir_1.0-9          munsell_0.5.0         tibble_3.2.1       
[55] pillar_1.9.0          htmltools_0.5.6       R6_2.5.1           
[58] wk_0.7.3              evaluate_0.21         lattice_0.22-5     
[61] backports_1.4.1       class_7.3-22          Rcpp_1.0.11        
[64] gridExtra_2.3         checkmate_2.2.0       spatstat.sparse_3.0-2
[67] mgcv_1.9-1            xfun_0.40             pkgconfig_2.0.3      
```
## Licence
CC BY 4.0
