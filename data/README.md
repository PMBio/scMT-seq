---
title: "Preprocessed data of scM&T"
author: "Christof Angermueller"
date: "2016-07-08"
output:
  html_document:
    toc: no
---

<style>
img {
    max-width: none;
}
</style>






```r
scMT <- readRDS('data.rds')
cols <- c('name', 'id_.x', 'id_.y', 'sample', 'met', 'weight', 'expr',
  'chromo.x', 'start.x', 'end.x', 'chromo.y', 'start.y', 'end.y', 'strand',
  'ens_id', 'gene_id')
scMT <- scMT[, c(cols)]
glimpse(scMT)
```

```
## Observations: 7,709,302
## Variables: 16
## $ name     (fctr) H3K27ac, H3K27ac, H3K27ac, H3K27ac, H3K27ac, H3K27ac...
## $ id_.x    (int) 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1...
## $ id_.y    (int) 9848, 9848, 9848, 9848, 9848, 9848, 9848, 9848, 9848,...
## $ sample   (fctr) CSCP3_SERUM_A02, CSCP3_SERUM_A03, CSCP3_SERUM_A04, C...
## $ met      (dbl) 70.588234, 100.000000, NaN, 8.333333, 16.666666, 0.00...
## $ weight   (dbl) 17, 12, 0, 12, 6, 6, 0, 5, 0, 40, 17, 0, 5, 15, 13, 3...
## $ expr     (dbl) 0.4506199, 1.2592708, 0.4329859, 0.4880455, 0.3933367...
## $ chromo.x (fctr) 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, ...
## $ start.x  (int) 108920349, 108920349, 108920349, 108920349, 108920349...
## $ end.x    (int) 108950783, 108950783, 108950783, 108950783, 108950783...
## $ chromo.y (fctr) 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, ...
## $ start.y  (int) 108909110, 108909110, 108909110, 108909110, 108909110...
## $ end.y    (int) 108912460, 108912460, 108912460, 108912460, 108912460...
## $ strand   (fctr) +, +, +, +, +, +, +, +, +, +, +, +, +, +, +, +, +, +...
## $ ens_id   (fctr) ENSMUSG00000000142, ENSMUSG00000000142, ENSMUSG00000...
## $ gene_id  (fctr) Axin2, Axin2, Axin2, Axin2, Axin2, Axin2, Axin2, Axi...
```

* `name`: Name of annotation.
* `id_.x`: Identifier of methylated region.
* `id_.y`: Identifier of gene matched to methylated region.
* `sample`: Cell identifier.
* `met`: Methylation rate of methylated region.
* `weight`: Weight of methylated region proportional to the number of covered
  CpG sites in that region.
* `expr: Expression rate of gene.
* `ens_id`: ENSEMBL gene identifier.
* `gene_id`: Gene name.
* `chromo.x`: Chromosome of methylated region.
* `start.x`: Start position of methylated region.
* `end.x`: End position of methylated region.
* `chromo.y`: Chromosome of gene.
* `start.y`: Start position of gene.
* `end.y`: End position of gene.
