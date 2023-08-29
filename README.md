# HLAfreq

`HLAfreq` allows you to download and combine HLA allele
frequencies from multiple datasets, e.g. combine data from
several studies within a country or combine countries.
Useful for studying regional diversity in immune genes
and, when paired with epitope prediction, estimating a population's
ability to mount an immune response to specific epitopes.

Automated download of allele frequency data download from 
[allele frequencies.net](http://www.allelefrequencies.net/).

## Details
Estimates are combined by modelling allele frequency as a 
Dirichlet distribution which defines the probability of drawing each
allele. When combining studies their estimates are weighted as 2x sample size by
default. Sample size is doubled as each person in the study
contributes two alleles. Alternative weightings can be used,
for example population size when averaging across countries.

When selecting a panel of HLA alleles to represent a population,
allele frequency is not the only thing to consider. Depending on
the purpose of the panel, you should include a range of loci and
supertypes (groups alleles sharing binding specificies).

## Install
`HLAfreq` is a `python` package available on windows, mac, and linux. We recommend installing
with `conda`.
```
conda create -n hlafreq -c bioconda -c conda-forge hlafreq
conda activate hlafreq
```
If you're new to conda see the miniconda [installation guide](https://conda.io/projects/conda/en/stable/user-guide/install/index.html) and [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/index.html)
to get started with `conda`.
Enter the above command into your conda prompt to create and
activate a conda environment with `HLAfreq` installed.
Typing `python` into this activated environment will start
a python session where you can enter your python code such as
the HLAfreq [minimal example](#minimal-example) below.

If you prefer to write your python code as scripts using an IDE such as
PyCharm or VScode, you'll need to look up how to configure a conda
virtual environment with those tools.

### Troubleshooting
`HLAfreq` uses `pymc` to estimate credible intervals,
which is the source of most installation difficulty, see
[pymc installation guide](https://www.pymc.io/projects/docs/en/stable/installation.html).

At time of writing `pymc` doesn't play nice with python 3.11, so
you can try installing a specific `python` version
and then add `HLAfreq` with pip or conda.
For example
```
conda create -n hlafreq
conda activate hlafreq
conda install -c conda-forge python=3.10
conda install -c conda-forge pymc
conda install -c bioconda hlafreq
```

`HLAfreq` requires `python>=3.8` and `pymc>=3`.
If you do run into trouble please open an [issue](https://github.com/Vaccitech/HLAfreq/issues).

If you don't intend to use credible intervals you can install
with pip: `pip install HLAfreq`.
However, if you do import `HLAfreq_pymc` you may get warnings
about degraded performance.

See the [pip documentation](https://pip.pypa.io/en/stable/)
to get started with pip. If you do have issues with pip,
try installing with conda as described above.

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
- [Credible intervals](https://github.com/Vaccitech/HLAfreq/blob/main/examples/credible_intervals.ipynb)

## Docs
For help on specific functions view the docstring, `help(function_name)`.
Full documentation API at [HLAfreq/docs](https://github.com/Vaccitech/HLAfreq/blob/main/docs/HLAfreq.md)
created with pdoc3 in pdf mode.

<!-- ## Developer notes
# Install in dev mode
pip install -e HLAfreq

Update version in setup.py

Update documentation with `pdoc --pdf -o docs/ src/HLAfreq/ > docs/HLAfreq.md`.

Run tests `pytest` 

# Clear old build info
rm -rf build dist src/*.egg-info 

Build with `python -m build`.

twine check dist/*

# Upload to test pypi
twine upload --repository testpypi dist/*

# Install from test pypi
python3 -m pip install --extra-index-url https://test.pypi.org/simple/ HLAfreq

# Upload to pypi
twine upload dist/*
-->

## Citation
*In prep.*
