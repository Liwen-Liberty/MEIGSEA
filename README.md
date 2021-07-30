
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MEIGSEA

<!-- badges: start -->
<!-- badges: end -->

The goal of MEIGSEA is to …

## Installation

You can install the released version of MEIGSEA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Liwen-Liberty/MEIGSEA")
```

## Usage

    MEIGSEA(raw.count.profile, is.rawcount = TRUE, gene.length = NULL, tpm.exp = NULL, mutationCalling.map, 

        signature.df, case = "1", control = "0", 
        
        cor.cutoff = NULL, cor.method = "spearman",
        
        min.sz = 15, perm.times = 100, pAdjustMethod = "BH", 
        
        cond.corr.anal = TRUE, sig.cutoff = 0.01, 
        
        output = TRUE, output.dir = getwd(), interested.mut = NULL)

## Arguments

**raw.count.profile** RNA-seq data in raw count or microarray expression
profile with row as genes, column as samples.

**is.rawcount** Whether or not the ‘raw.count.profile’ is RNA-seq raw
count. For microarray expression data as input, is.rawcount = FALSE.

**gene.length** Gene length annotation file, the first column is the
gene name, the second column is the gene length, if ‘gene.length’ is not
provided, the ‘tpm.exp’ needs to be given.

**tpm.exp** Given tpm expression profile (optional). If the
‘gene.lengths’ given, tpm.exp is not needed.

**mutationCalling.map** A matrix of gene mutation status (0 or 1, in
character), row as genes, column as samples.

**signature.df** signature A data frame containing immunophenotype
signature: the first column (named ‘Gene’) is the gene ID, the second
column (named ‘setAnno’) is the immunophenotype name.

**case/control** In default, case = “1”, control = “0”, or the specific
character representing mutant (case) and wild type (control) according
to ‘mutationCalling.map’.

**cor.cutoff** Setting the cutoff of correlation coefficient when
refining the immunophenotype, cor.cutoff = NULL by default, means all
genes with statistical significance will be retained

**cor.method** Method used to compute the correlation, either
“spearman”, “pearson” or “kendall”, “spearman” by default

**min.sz** At least how many immunophenotype genes are matched in the
expression profile, min.sz = 15 by default.

**perm.times** Perumutation times when the normalization of ssGSEA score
is needed, 100 by default.

**pAdjustMethod** Method used to adjust p values, “BH” by default.

**cond.corr.anal** Whether to perform conditional association analysis,
cond.corr.anal=TRUE by default.

**sig.cutoff** Cutoff for significant associations, sig.cutoff = 0.01 by
default.

**output** Whether to store the results locally, output = TRUE by
default.

**output.dir** Character representing the storage path of results.

**interested.mut** Vector of interested gene mutations, for which the
detailed correlation results will be given as independent files

## Value

A list including detailed corrlation results, refined immunophenotype
signature and all correlated mutation-signature pairs

## Author(s)

Liwen Xu, Shiwei Zhu

## Example

``` r
library(MEIGSEA)
data(THCA_MutCountTPM)
data(CIBERSORT_cell_types_markers)

## Mutation status data
THCA_MutCountTPM$mutationCalling.map[,1:5]
#>      TCGA-DJ-A3UX TCGA-EL-A3D0 TCGA-EM-A22Q TCGA-E3-A3E2 TCGA-ET-A3BW
#> BRAF "1"          "1"          "0"          "1"          "1"         
#> HRAS "0"          "0"          "0"          "0"          "0"         
#> NRAS "0"          "0"          "1"          "0"          "0"

## Predetermined immunophenotype signatures
head(CIBERSORT_cell_types_markers)
#>       Gene        setAnno    Source
#> 553 ADAM28 B.cells.memory CIBERSORT
#> 558   AIM2 B.cells.memory CIBERSORT
#> 560  ALOX5 B.cells.memory CIBERSORT
#> 577  BACH2 B.cells.memory CIBERSORT
#> 578  BANK1 B.cells.memory CIBERSORT
#> 587    BLK B.cells.memory CIBERSORT

## Running the MEIGSEA 
#MEIGSEA.res <- MEIGSEA(raw.count.profile =THCA_MutCountTPM$rawCount, 
#                       is.rawcount = TRUE, 
#                       gene.length = NULL, 
#                       tpm.exp = THCA_MutCountTPM$TPM, 
#                       mutationCalling.map = THCA_MutCountTPM$mutationCalling.map,
#                       signature.df=CIBERSORT_cell_types_markers)
```
