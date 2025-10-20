# SpaceTrooper

SpaceTrooper is a R/Bioconductor package for the preprocessing and quality 
control of imaging-based spatial transcriptomics and proteomics data.

<p align="center">

<img src="https://raw.githubusercontent.com/drighelli/SpaceTrooper/main/inst/SpaceTrooper_logo.png" alt="SpaceTrooper Logo" width="200"/>

</p>

The package is in Bioconductor since version 3.22, and it can be installed by 
aid of the `BiocManager` package, if you are using the proper Bioconductor 
version (3.22).

By executing `BiocManager::install("SpaceTrooper)` you will install the 
latest version of the package on Bioconductor.

# Installing via Github

## Requirements

On te other hand, if you want to install the package from Github, it could be 
unstable and R version must be at least 4.4.0.

We always suggest to use the latest version of the package present on the 
official Release of Bioconductor.

## Installation

First try to run the following either in R or in a virgin environment:

`BiocManager::install("drighelli/SpaceTrooper", ref="main")`

If you want a not stable version of the package you can use:

`BiocManager::install("drighelli/SpaceTrooper", ref="devel")`

If you were able to install the package correctly, CONGRATULATIONS!!! 🎉🎉🎉

# Introduction to the package 

We provide two additional vignettes for the package usage, first one is focused
on imaging-based Spatial Transcriptomics technologies.

The second one is focused on CosMx protein technologies.

See the compiled version of the vignettes in the 
[official Bioconductor landing page](https://bioconductor.org/packages/devel/bioc/html/SpaceTrooper.html)
or visit the [Vignette folder](https://github.com/drighelli/SpaceTrooper/tree/main/vignettes).

# Issues

Please report any issue on GH [here](https://github.com/drighelli/SpaceTrooper/issues)!

Thanks for checking out!🌸
