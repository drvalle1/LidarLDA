
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Authors

Denis Valle, Carlos Silva, Marcos Longo, Paulo Brando

# LidarLDA

<!-- badges: start -->
<!-- badges: end -->

The goal of LidarLDA is to fit a modified version of the Latent
Dirichlet Allocation (LDA) model to LiDAR data. The main benefit of
using this mixed-membership model is that it allows for some grid cells
to have varying proportions of each cluster whereas more standard
hard-clustering methods force grid cells to belong to a single cluster.
This model estimates two sets of parameters: one set characterizes each
grid cell in relation to the relative abundance of clusters while the
other set characterizes each cluster in relation to its absorptance
probabilities.

## Installation

You can install this package from [GitHub](https://github.com/) with:

``` r
library("devtools")
devtools::install_github("drvalle1/LidarLDA",build_vignettes=T)
```

## Fitting LidarLDA to simulated data

We start by showing how to fit the model based on simulated data with 5
clusters. The simulated datasets contain information for 2,601 pixels
(rows) and 30 height bins (columns), labeled z1, z2, …, z30.

The data in `sim_y5` consist of the number of returned light pulses
whereas the data in `sim_n5` consist of the number of incoming light
pulses, in each pixel and height bin. As a result, `sim_y5` is always
smaller or equal to `sim_n5`. These data were simulated with 5 clusters.

``` r
#library('devtools')
#devtools::install_github("drvalle1/LidarLDA",build_vignettes=T)
library(LidarLDA)

#basic characteristics of simulated data
dim(sim_y5) #sim_y5 is provided in the LidarLDA package
#> [1] 2000   50
dim(sim_n5) #sim_n5 is provided in the LidarLDA package
#> [1] 2000   50

colnames(sim_y5)
#>  [1] "z1"  "z2"  "z3"  "z4"  "z5"  "z6"  "z7"  "z8"  "z9"  "z10" "z11" "z12"
#> [13] "z13" "z14" "z15" "z16" "z17" "z18" "z19" "z20" "z21" "z22" "z23" "z24"
#> [25] "z25" "z26" "z27" "z28" "z29" "z30" "z31" "z32" "z33" "z34" "z35" "z36"
#> [37] "z37" "z38" "z39" "z40" "z41" "z42" "z43" "z44" "z45" "z46" "z47" "z48"
#> [49] "z49" "z50"
colnames(sim_n5)
#>  [1] "z1"  "z2"  "z3"  "z4"  "z5"  "z6"  "z7"  "z8"  "z9"  "z10" "z11" "z12"
#> [13] "z13" "z14" "z15" "z16" "z17" "z18" "z19" "z20" "z21" "z22" "z23" "z24"
#> [25] "z25" "z26" "z27" "z28" "z29" "z30" "z31" "z32" "z33" "z34" "z35" "z36"
#> [37] "z37" "z38" "z39" "z40" "z41" "z42" "z43" "z44" "z45" "z46" "z47" "z48"
#> [49] "z49" "z50"

#remove columns with coordinates and topography information
ind=which(colnames(sim_y5)%in%c('x','y','topo'))
sim_y5a=sim_y5[,-ind]
coord=sim_y5[,ind]

ind=which(colnames(sim_n5)%in%c('x','y','topo'))
sim_n5a=sim_n5[,-ind]
mean(sim_y5a<=sim_n5a)
#> [1] NaN
```

We fit these simulated data using the code below. In this code, we
assume a maximum of 10 clusters and we rely on 10000 iterations of the
gibbs sampler with a burn-in of 9000 iterations. Finally, we just return
the posterior mean parameter estimates instead of all the posterior
samples by specifying `theta.post=F` and `phi.post=F`.

``` r
Model.Results=LidarLDA(y=data.matrix(sim_y5a),
                       n=data.matrix(sim_n5a),
                       nclust=10,
                       a.phi=1,b.phi=1,
                       gamma=0.1,ngibbs=10000,
                       nburn=9000,theta.post=F,phi.post=F)
```

We can assess convergence by examining the trace-plot of the
log-likelihood. This plot suggests that the algorithm has converged.

``` r
plot(Model.Results$llk,type='l',xlab='Iterations',
     ylab='Log-likelihood')
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

According to the `theta` matrix (i.e., the matrix that shows the
relative abundance of each cluster for each pixel), our model has
identified 5 (out of a maximum of 10) main clusters. These 5 first
clusters, on average, represent 99.8% of all observations in each pixel.

``` r
boxplot(Model.Results$theta,xlab='Cluster id',ylab='theta')
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
sum1=apply(Model.Results$theta[,1:5],1,sum)
mean(sum1)
#> [1] 0.9977585
```

Because this is based on simulated data, there is a nice pattern
regarding how the relative abundance of each cluster changes as a
function of topography. This is shown below.
