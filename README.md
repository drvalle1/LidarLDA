
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LidarLDA

<!-- badges: start -->
<!-- badges: end -->

The goal of LidarLDA is to fit a modified version of the Latent
Dirichlet Allocation (LDA) model to LIDAR data. This model will a)
characterize each pixel in relation to the relative abundance of
clusters; and b) characterize each cluster in relation to its
absorptance probabilities.

## Installation

You can install this package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("drvalle1/LidarLDA")
```

## Example

This is a basic example which shows how to run LidarLDA:

``` r
library(LidarLDA)

## basic example code
P=1000 #number of pixels
H=10   #number of height bins

#data
y=matrix(1,P,H) #number of pulses returned in each pixel and height bin
n=matrix(2,P,H) #number of incoming pulses in each pixel and height bin
nclust=10       #maximum number of clusters

Model.Results=LidarLDA(y=y,
                       n=n,
                       nclust=nclust,
                       a.phi=1,b.phi=1,
                       gamma=0.1,ngibbs=100,
                       nburn=50,theta.post=F,phi.post=F)
```

A more elaborate example can be found in the vignette of this package.
