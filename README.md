# DynamicPGS
Computing dynamic polygenic scores across time using Gaussian process regression models

## 0. Constructing population average model (optional)

We provide a population average model of BMI along child age. If you think your target population is different to the Japanese population, you can create your own model to provide an average of the target trait in your population over time. 

$y=(y_1^\top,\ldots,y_N^\top)^\top$ where $y_{i}=(y_{i1},\ldots,y_{i,n_i})^\top$

$x=(x_1^\top,\ldots,x_N^\top)^\top$ where $x_{i}=(x_{i1},\ldots,x_{i,n_i})^\top$

## 1. Computing dynamic PGS using our dynamic PGS model

```
getDynamicPGS(x=x, pgen_dir="/path/to/your/pgen/dir/", avg_model=avg_jecs)
```




