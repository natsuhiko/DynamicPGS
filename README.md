# DynamicPGS
Computing dynamic polygenic scores across time using Gaussian process regression models

## 0. Constructing population average model (optional)

We provide a population average model of BMI along child age. If you think your target population is different to the Japanese population, you can create your own model to provide an average of the target trait in your population over time. 

$y=(y_1^\top,\ldots,y_N^\top)^\top$ where $y_{i}=(y_{i1},\ldots,y_{i,n_i})^\top$

$x=(x_1^\top,\ldots,x_N^\top)^\top$ where $x_{i}=(x_{i1},\ldots,x_{i,n_i})^\top$

```
avg = getAvg(y=y, x=x, id=id, X=X)
```

## 1. Computing dynamic PGS using our dynamic PGS model

This is how we estimate PGS at child age $x^*$ in month using their genotype data (given by the plink2 PGEN format). To compute PGS at age of 0 to 54 months, you can simply run:
```
# you need to stay at the package home dir
source("Data/getDynamicPGS.R")
pgs = getDynamicPGS(xstar=0:54, pgen_dir="/path/to/your/pgen/dir/")
```
Here PGEN files must be named as
```
chr1.pgen chr2.pgen ... chr22.pgen chrX.pgen
```
It is also the case that no genotype data on X chromosome is available, the R code simply skips the missing chromosome(s) to compute PGS. 

The output object `pgs` contains 4 branches, `xstar`, `avg`, `E` and `SE`. First, you can check the population average estimated from our data by
```
plot(pgs$xstar, pgs$avg)
```  
then you might want to check PGS for the first 10 individual by
```
matplot(pgs$xstar, pgs$E[,1:10], type="l")
```
Note that it accepts any number, including decimal value (`xstar=11:19/10` for 1.1 to 1.9 months old), and can predict PGS outside of the constructed GP model (although not recommended).



