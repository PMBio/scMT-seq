scM&T-Seq
=========

Source code of the manuscript ***Parallel single-cell sequencing links
transcriptional and epigenetic heterogeneity*** ([Nature Methods](http://www.nature.com/nmeth/journal/v13/n3/full/nmeth.3728.html)).

Abstract
--------
We report scM&T–seq, a method for parallel single–cell genome–wide methylome and
transcriptome sequencing, allowing discovery of associations between
transcriptional and epigenetic variation. Profiling of 61 mouse embryonic stem
cells confirmed known links between DNA methylation and transcription. Notably,
the method reveals novel associations between heterogeneous methylation of
distal regulatory elements and transcriptional heterogeneity of key pluripotency
genes.

Content
-------
* `/cca/`: Canonical Correlation Analysis
* `/clust/`: Clustering scM&T-Seq and scBS-Seq cells
* `/cov/`: Coverage analysis
* `/cor_bulk/`: Correlation scM&T-Seq methylation rates with bulk methylation rate
* `/data/`: Data directory
* `index.Rmd`: Table of content
* `/gene/`: Gene-specific correlation analysis
* `/gene_mean/`: Correlating mean methylation with gene expression
* `/gene_robust/`: Robustness analysis gene-specific correlation
* `/heat/`: Visualizing methylation and expression heatmap
* `/lib/`: Library functions
* `/sample/`: Sample-specific correlation analysis
* `/qc/`: Quality control DNA methylation
* `/var/`: Comparison methylation variability in context
* `/zoom/`: Visualizing Esrrb gene

`data/join/data.rds` contains the pre-processed and joined methylation and
expression data, which were used for the correlation analysis reported in the
manuscript. The raw data and intermediate output files can be downloaded
from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74535).

Contact
-------
* Christof Angermueller
* cangermueller@ebi.ac.uk
* https://cangermueller.com
