# HLAfreq

The aim of this tool is to be able to estimate the HLA allele
frequencies at large multi population scale, e.g. combine data from
several studies within a country or combine countries.
Allele frequency data from
[allele frequencies.net](http://www.allelefrequencies.net/).
This gives frequencies by population, which can be combined across
studies. Estimates are combined by modelling allele frequency as a 
Dirichlet distribution which defines the probability of drawing each
allele. When combining studies are weighted by 2 x sample size by
default. Sample size is doubled as each person in the study
contributes two alleles. Alternative weightings can be used
for example population size when averaging across countries.

```
import sys
sys.path.append("/home/dwells/work/tools/scrapeAF/code")
import scrapeAF
```

When selecting alleles to represent a population, all supertypes
should be represented.

When combining estimates of allele frequency larger samples should
weigh more heavily. However, when combining estiamtes that disagree
the uncertainty on the estimate should be larger than when combining
equivalent samples that agree.

When combining studies we must weight averages by sample size but
when combining countries are we not better to weight by national population?
Should sample size also have an effect?

IEDB provide reference sets of HLA alleles designed to cover the
majority of the "general population". However, all populations are
not likely to be equally well covered by them. Using population
specific allele frequencies we should be able to select more
appropriate HLA sets. However, allele frequency is not the only
determinent when selecting a pannel of alleles. The pannel must
also contain representatives of commonly shared binding specificities
(i.e., supertypes).

# Dependencies
scrapeAF
umap
maps
