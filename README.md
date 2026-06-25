# DynamicPGS

DynamicPGS is an R package for modelling longitudinal phenotypes with Gaussian-process regression, performing dynamic genetic association testing, and computing dynamic polygenic scores (dynamic PGS) over a continuous index such as age or time.

The typical workflow is:

```r
getData() -> gpreg1() -> gpreg2() -> getP() -> getDynamicPGS()
```

## Installation

Install the package from a local source directory:

```bash
R CMD INSTALL DynamicPGS
```

During development, use `devtools::load_all()` from the package root:

```r
devtools::load_all()
```

Do not source individual files such as `source("R/getData.R")` during package development, because this may create conflicts between objects in the global environment and functions loaded from the package.

## Package structure

Recommended source-tree structure:

```text
DynamicPGS/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── class.R
│   ├── data.R
│   ├── gpreg1.R
│   ├── gpreg2.R
│   ├── prep_assoc.R
│   ├── getP.R
│   ├── getDynamicPGS.R
│   ├── getDoseFromVCF.R
│   └── util.R
├── man/
├── README.md
└── tests/
    └── testthat/
```

Help files under `man/` should be generated from roxygen2 comments:

```r
devtools::document()
```

Then functions can be documented and called as:

```r
?getData
?gpreg1
?gpreg2
?getP
?getDynamicPGS
?getDoseFromVCF
```

## Input data

### Longitudinal phenotype data

The phenotype table must contain at least the following columns:

| Column | Description |
|---|---|
| `IID` | Individual ID |
| `x` | Continuous index, for example age in months |
| `y` | Phenotype value |

Example:

```text
IID     x      y
id001   6      16.2
id001   12     18.4
id002   6      15.8
id002   18     20.1
```

### Covariates

Covariates can be supplied as a separate data.frame or file. Numeric covariates are standardised internally. Character or factor covariates are expanded into dummy variables.

The number of rows in the covariate table should match the number of rows in the phenotype table.

### Relatedness information

If a KING relatedness file is supplied, DynamicPGS uses the relatedness information to construct family blocks. The KING file should contain at least:

| Column | Description |
|---|---|
| `ID1` | First individual ID |
| `ID2` | Second individual ID |
| `Kinship` | KING kinship estimate |

If no KING file is supplied, all individuals are treated as unrelated.

## Basic workflow

### 1. Create a DynamicPGS object

```r
library(DynamicPGS)

adata <- getData(
  Data = "phenotype.tsv",
  Covariates = "covariates.tsv",
  king_file = "king.kin0",
  inducing_points = seq(0, 60, by = 6)
)
```

The returned object has class `DynamicPGS`.

```r
adata
```

This prints a short summary of the number of observations, number of individuals, maximum family size, support of the continuous index, and covariate structure.

### 2. Fit the population-level GP model

```r
adata <- gpreg1(
  adata,
  rho_a = 0.01,
  rho_b = 5,
  Verbose = TRUE,
  Plot = FALSE,
  MAXITR = 100
)
```

This step estimates the population-level smooth trajectory, covariate variance components, GP length-scale parameter, and residual variance.

### 3. Fit individual-level dynamic deviation components

```r
adata <- gpreg2(
  adata,
  delta2d0 = c(0.94, 1.7),
  ncore = 4,
  Verbose = TRUE,
  Plot = FALSE,
  MAXITR = 100
)
```

This step estimates individual-level dynamic deviation components. It also calls `prep_assoc()` internally to prepare matrices for association mapping.

If needed, association-mapping quantities can be recomputed explicitly:

```r
adata <- prep_assoc(
  adata,
  r_rho = 1,
  r_delta2d = 1,
  ncore = 4
)
```

## Genotype dosage input

### From an indexed VCF

DynamicPGS can read genotype dosages from a bgzip-compressed and tabix-indexed VCF file.

```r
Gall <- getDoseFromVCF(
  vcf = "imputed.vcf.gz",
  region = "chr1:100000-200000"
)
```

The returned object is a dosage matrix:

```r
dim(Gall)
head(rownames(Gall))
head(colnames(Gall))
```

Rows are variants and columns are samples. Variant row names are constructed as:

```text
CHROM:POS:REF:ALT
```

Variant metadata are stored as an attribute:

```r
attr(Gall, "variant_info")
```

If the VCF contains a `DS` field, dosages are read from `DS`. If `DS` is absent, `GT` is converted to alternate-allele dosage.

### Expected dosage matrix format

Association testing expects a numeric matrix with variants in rows and individuals in columns:

```text
          id001  id002  id003
variant1  0.02   1.01   2.00
variant2  1.00   0.00   0.98
```

Column names should correspond to individuals in the fitted `DynamicPGS` object.

## Dynamic association testing

Run variant-level dynamic association testing:

```r
adata <- getP(
  adata,
  Gall = Gall,
  delta2g = 0.01,
  Beta = TRUE,
  Sinv = TRUE,
  ncore = 4
)
```

The association p-values are stored in:

```r
head(adata$pval)
```

If `Beta = TRUE`, posterior dynamic effect estimates are stored in:

```r
adata$Beta
```

If `Sinv = TRUE`, inverse covariance matrices are stored in:

```r
adata$Sinv
```

`Beta = TRUE` and `Sinv = TRUE` are recommended if the next step is to compute dynamic PGS with uncertainty estimates.

## Computing dynamic PGS

After association testing, dynamic PGS can be evaluated at arbitrary values of the continuous index.

For example, to compute dynamic PGS from 0 to 60 months:

```r
dpgs <- getDynamicPGS(
  xstar = 0:60,
  Gall = Gall,
  adata = adata
)
```

The output is a list with the following elements:

| Element | Description |
|---|---|
| `xstar` | Target index values |
| `avg` | Estimated population-average trajectory |
| `E` | Dynamic PGS matrix, with rows corresponding to `xstar` and columns to samples |
| `SE` | Approximate standard errors of dynamic PGS |
| `sigma2` | Residual variance estimate |

Plot dynamic PGS for the first 10 individuals:

```r
matplot(dpgs$xstar, dpgs$E[, 1:10], type = "l",
        xlab = "x", ylab = "Dynamic PGS")
```

Plot a dynamic PGS with an approximate 95% interval:

```r
i <- 1
upper <- dpgs$E[, i] + 1.96 * dpgs$SE[, i]
lower <- dpgs$E[, i] - 1.96 * dpgs$SE[, i]

plot(dpgs$xstar, dpgs$E[, i], type = "l",
     xlab = "x", ylab = "Dynamic PGS")
polygon(
  c(dpgs$xstar, rev(dpgs$xstar)),
  c(upper, rev(lower)),
  col = rgb(1, 0, 0, 0.1),
  border = NA
)
lines(dpgs$xstar, dpgs$E[, i])
```

Plot predicted phenotype values with an approximate prediction interval:

```r
i <- 1
pred <- as.numeric(dpgs$avg + dpgs$E[, i])
upper <- pred + 1.96 * sqrt(dpgs$sigma2 + dpgs$SE[, i]^2)
lower <- pred - 1.96 * sqrt(dpgs$sigma2 + dpgs$SE[, i]^2)

plot(dpgs$xstar, pred, type = "l",
     xlab = "x", ylab = "Predicted phenotype")
polygon(
  c(dpgs$xstar, rev(dpgs$xstar)),
  c(upper, rev(lower)),
  col = rgb(1, 0, 0, 0.1),
  border = NA
)
lines(dpgs$xstar, pred)
```

## Development notes

### Updating documentation

Add roxygen2 comments immediately above each exported function, then run:

```r
devtools::document()
```

This updates:

```text
NAMESPACE
man/*.Rd
```

### Loading during development

Use:

```r
devtools::load_all()
```

rather than sourcing files manually.

If you accidentally sourced functions and see messages such as:

```text
getData masks DynamicPGS::getData()
gpreg1 masks DynamicPGS::gpreg1()
```

remove the global objects:

```r
rm(list = c("getData", "gpreg1"))
devtools::load_all()
```

or restart the R session and run `devtools::load_all()` again.

### Checking the package

Run:

```r
devtools::check()
```

Before submission or release, also check that all exported functions have examples and help pages.

## Main exported functions

| Function | Purpose |
|---|---|
| `getData()` | Create a `DynamicPGS` object from longitudinal phenotype, covariate, and relatedness data |
| `gpreg1()` | Fit the population-level Gaussian-process regression model |
| `gpreg2()` | Fit individual-level dynamic deviation components |
| `prep_assoc()` | Prepare matrices for dynamic association mapping |
| `getP()` | Perform dynamic association testing for genotype dosages |
| `getDynamicPGS()` | Compute dynamic PGS over a continuous index |
| `getDoseFromVCF()` | Extract genotype dosages from a tabix-indexed VCF |
| `print.DynamicPGS()` | Print a summary of a `DynamicPGS` object |

## Dependencies

DynamicPGS currently uses the following R packages:

```r
Matrix
CompQuadForm
parallel
```

`tabix` is required when using `getDoseFromVCF()` with VCF files.

## Citation

Citation information will be added after manuscript or package release.
