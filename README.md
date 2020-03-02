## Synchrony and idiosyncrasy in the gut microbiome of wild primates

#### Load data
The analyzed data (i.e. the ASV table, metadata and taxonomy) is stored in a phyloseq object.

You can load this data in `R` by typing

```
library(phyloseq)
ps <- loadRDS("ps.RDS")
```

If you want to side step phyloseq specific `R` code and work directly on the data, you can extract each component by typing 

```
ftbl <- as(otu_table(ps),"matrix") 
metadata <- as(sample_data(ps), "data.frame")
tax <- as(tax_table(ps2),"matrix")
```

#### Dealing with compositionality
To account for the compositional nature of 16S amplicon sequence data and differences in sequencing depth, we followed the recommendations of Gloor et al. (2017) and Morton et al. (2019) and used the centered log-ratio (clr) transforation (see also Pawlowsky-Glahn & Egozcue, 2006). An important property of the clr-transformation is scale invariance; this means that the same ratio is expected in two identical samples--one with a few reads, and one with many reads (Gloor et al. 2017). An important consequence of this property is that count normalization is unnecessary as it only leads to a loss of information and precision (Gloor et al. 2017; McMurdie & Holmes, 2014). Lastly, the Euclidean distance between clr-transformed samples--the Aitchison distance--has been shown to be superior to other commonly used distance metrics such as Bray-Curtis and Jensen-Shannon divergence when applied on compositional data (Martino et al. 2019; Aitchison et al. 2000). 

The first three microbiome PCs (`PC1`, `PC2`, `PC3`) you can find in the data.frame `metadata` are computed from the below `R` code.     

```
pseudocount <- 0.65
clr_pca <- function(ftbl) {
  dclr <- t(apply(t(ftbl)+pseudocount, 2, compositions::clr))
  pcx <- prcomp(dclr)
  out <- list(clr=dclr, pca=pcx)
  return(out)
}
```

### Microbiome features
The `metadata` contain 6 different microbiome features. 
* `PC1-PC3:` the first three PCs from clr-transformed realtive abundances (see above).
* `richness; shannon; simpson:` alpha diversity indicies corresponding to the first three Hill numbers, i.e. Species richness, Shannon entropy and the inverse Simpson index.

[NOTE LOG-RATIOS ARE IN DIFFERENT DATASET?]

#### Covariates
At the **population level**, `metadata` contain the columns `rain_monthly`, `rain_annual`, `tempmax_monthly` and `tempmax_annual` for total monthly and annual rainfall and maximum monthly and annual temperature. 

At the **social group level**, `metadata` contain 18 columns: 
* `goup_size:` each group's total number of individuals.
* `lon:` each group's monthly average longitude. 
* `lat:` each group's monthly average latitude.
* `area_tot_sqkm:` each group's monthly total home range area (km2).
* `frac_uniq:` each group's non-overlapping home range area (0-1).
* `diet_PC1-diet_PC13:` diet composition was computed for each social group and fecal collection date combination. More specifically, we computed the average relative abundance of each diet category from 1,000 random subsamples of all foraging observations recorded one month prior to the collection date. We then performed a PCA on clr-transformed diet compositions. 

At the **host level**, `metadata` contain 3 columns: 
* `age:` an individual host's age (in years).
* `sex:` an individual host's sex.
* `proportional_rank:` an individual host's social rank as determined by assigning wins and losses in dyadic agonistic interactions between animals that live in the same social group.


#### Analyses

Below you can find a brief description of the analyses in each **R Notebook:**

* `hierarchical_smooths.Rmd:` code to reproduce plots showing dynamics (in terms of smooth functions) over time at the population, group and host level.

* `taxa_balances.Rmd:` code to reproduce plots showing log-ratios (also called balances) of the most abundant phyla and families plotted over time. **[STILL HAVE TO UPLOAD!]**

* `time_decay.Rmd:` code to reproduce plots showing time-decay in community similarity for samples from (1) the same host in the same group; (2) different hosts in the same group; and (3) different hosts in different groups.    

* `gams.Rmd:` code to fit our hierarchical Generalized Additive Models (GAMs).

* `gams_sanity_checks.Rmd:` code for a series of GAM sanity checks. Especially to show that the increase in deviance explained is not an artifact from data aggregation or increasing model complexity.

 

