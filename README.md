# HLAfreq

`HLAfreq` allows you to download and combine HLA allele
frequencies from multiple datasets, e.g. combine data from
several studies within a country or combine countries.

Automated download of allele frequency data download from 
[allele frequencies.net](http://www.allelefrequencies.net/).

## Details
Estimates are combined by modelling allele frequency as a 
Dirichlet distribution which defines the probability of drawing each
allele. When combining studies their estimates are weighted as 2x sample size by
default. Sample size is doubled as each person in the study
contributes two alleles. Alternative weightings can be used
for example population size when averaging across countries.

A note on credible intervals... **(to be written)**

When selecting a panel of HLA alleles to represent a population,
allele frequency is not the only thing to consider. Depending on
the purpose of the panel, you should include a range of loci and
supertypes (groups alleles sharing binding specificies).

## Install
```
pip install HLAfreq
```

## Minimal example
Download HLA data using `makeURL()` and `getAFdata()`.
All arguments that can be specified in the webpage form are available,
see `help(HLAfreq.makeURL)` for details (press `q` to exit).
```
import HLAfreq
base_url = HLAfreq.makeURL("Uganda", locus="A")
aftab = HLAfreq.getAFdata(base_url)
```

After downloading the data, it must be filtered so that all studies
sum to allele frequency 1 (within tolerence). Then we must ensure
that all studies report alleles at the same resolution.
Finaly we can combine frequency estimates.
```
aftab = HLAfreq.only_complete(aftab)
aftab = HLAfreq.decrease_resolution(aftab, 2)
caf = HLAfreq.combineAF(aftab)
```

## Detailed examples
For more detailed walkthroughs see [HLAfreq/examples](https://github.com/Vaccitech/HLAfreq/tree/main/examples).

- [Single country](https://github.com/Vaccitech/HLAfreq/blob/main/examples/single_country.ipynb) download and combine
- [Multi-country](https://github.com/Vaccitech/HLAfreq/blob/main/examples/multi_country.ipynb) download and combine, weight by population coverage
- [Using priors](https://github.com/Vaccitech/HLAfreq/blob/main/examples/working_with_priors.ipynb)

## Docs
For help on specific functions view the docstring, `help(function_name)`.
Full documentation API at [HLAfreq/docs](https://github.com/Vaccitech/HLAfreq/blob/main/docs/HLAfreq.md)
created with pdoc3 in pdf mode.

<!-- ## Developer notes
Update documentation with `pdoc --pdf -o docs/ src/HLAfreq/ > docs/HLAfreq.md`.
Build with `pythom -m build`.
Run tests `pytest` -->

## Citation
*In prep.*
