# HLAfreq

`HLAfreq` allows you to download and combine HLA allele
frequencies from multiple datasets, e.g. combine data from
several studies within a country or combine countries.
Useful for studying regional diversity in immune genes
and, when paired with epitope prediction, estimating a population's
ability to mount an immune response to specific epitopes.

Automated download of allele frequency data download from 
[allelefrequencies.net](http://www.allelefrequencies.net/).

Full documentation at [HLAfreq/docs](https://BarinthusBio.github.io/HLAfreq/HLAfreq.html). Source code is available at [BarinthusBio/HLAfreq](https://github.com/BarinthusBio/HLAfreq).

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
supertypes (grouped alleles sharing binding specificies).

## Install
`HLAfreq` is a `python` package available on windows, mac, and linux. We recommend installing
with `conda`.
```
conda create -n hlafreq -c conda-forge -c bioconda hlafreq
conda activate hlafreq
```

### Troubleshooting
`HLAfreq` uses `pymc` to estimate credible intervals,
which is the source of most installation difficulty, see
[pymc installation guide](https://www.pymc.io/projects/docs/en/stable/installation.html) and [tips and tricks](https://conda-forge.org/docs/user/tipsandtricks/#using-multiple-channels).

You may see an error about g++ and degraded performance:
```
WARNING (pytensor.configdefaults): g++ not detected!  PyTensor will be unable to compile C-implementations and will default to Python. Performance may be severely degraded. To remove this warning, set PyTensor flags cxx to an empty string.
```

This means that one of the pymc backends is missing and estimating confidence
intervals will be very slow. But don't worry, try one of these fixes below:

- Set the channel priority to strict, then install as above (using conda-forge then bioconda channels).
```
conda config --set channel_priority strict
```

- Install a conda compiler to handle g++ based on your os.
```
conda create -n hlafreq -c conda-forge -c bioconda hlafreq cxx-compiler
```

When running entire scripts on windows, you may see an error about
"Safe importing of main module", multiprocessing, and starting
new processes. To fix this, main guard your code with
`if __name__ == "__main__":` after the `import`s as demonstrated in
[`examples/quickstart.py`](https://github.com/BarinthusBio/HLAfreq/blob/main/examples/quickstart.py).

If you do run into trouble please open an [issue](https://github.com/BarinthusBio/HLAfreq/issues).

### conda
If you're new to conda see the miniconda [installation guide](https://conda.io/projects/conda/en/stable/user-guide/install/index.html) and [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/index.html)
to get started with `conda`.

Enter the install command from above into your conda prompt to create and
activate a conda environment with `HLAfreq` installed.
Typing `python` into this activated environment will start
a python session where you can enter your python code such as
the HLAfreq [minimal example](#minimal-example) below.

If you prefer to write your python code as scripts using an IDE such as
PyCharm or VScode, you'll need to look up how to configure a conda
virtual environment with those tools.

### pip
If you don't intend to use credible intervals you can install
with pip: `pip install HLAfreq`.
However, if you do import `HLAfreq_pymc` you may get warnings
about degraded performance.

See the [pip documentation](https://pip.pypa.io/en/stable/)
to get started with pip. If you do have issues with pip,
try installing with conda as described above.

## Minimal example
Download HLA data using `HLAfreq.HLAfreq.makeURL()` and `HLAfreq.HLAfreq.getAFdata()`.
All arguments that can be specified in the webpage form are available,
see the [`makeURL()` docs](https://barinthusbio.github.io/HLAfreq/HLAfreq/HLAfreq.html#makeURL) for details.
```
import HLAfreq
base_url = HLAfreq.makeURL("Uganda", locus="A")
aftab = HLAfreq.getAFdata(base_url)
```

After downloading the data, it must be filtered so that all studies
sum to allele frequency 1 (within tolerence). Then we must ensure
that all studies report alleles at the same resolution.
Finaly we can combine frequency estimates, for more details see
the [`combineAF()` api documentation](https://barinthusbio.github.io/HLAfreq/HLAfreq/HLAfreq.html#combineAF).
```
aftab = HLAfreq.only_complete(aftab)
aftab = HLAfreq.decrease_resolution(aftab, 2)
caf = HLAfreq.combineAF(aftab)
```

To add confidence intervals to estimates see
[`examples/quickstart.py`](https://github.com/BarinthusBio/HLAfreq/blob/main/examples/quickstart.py).

## Detailed examples
For more detailed walkthroughs see [HLAfreq/examples](https://barinthusbio.github.io/HLAfreq/HLAfreq/examples.html).

- [Single country](https://BarinthusBio.github.io/HLAfreq/HLAfreq/examples/single_country.html) download and combine
- [Multi-country](https://BarinthusBio.github.io/HLAfreq/HLAfreq/examples/multi_country.html) download and combine, weight by population coverage
- [Using priors](https://BarinthusBio.github.io/HLAfreq/HLAfreq/examples/working_with_priors.html)
- [Credible intervals](https://BarinthusBio.github.io/HLAfreq/HLAfreq/examples/credible_intervals.html)

## Docs
Full documentation at [HLAfreq/docs](https://BarinthusBio.github.io/HLAfreq/HLAfreq.html).
API documentation for functions are under the submodules on the left.
- `HLAfreq.HLAfreq` documents most functions, specifically download and combine
allele data.
- `HLAfreq.HLAfreq_pymc` is functions using pymc to acurately estimate credible intervals on allele frequency estimates.

For help on specific functions view the docstring, `help(function_name)`.

Run `pdoc -d google -o docs/ HLAfreq` to generate the
documentation in `./docs`.
<!-- Documentation generated by pdoc should not be commited
as it is auto generated by a github action. -->


<!-- ## Developer notes
Install in dev mode
pip install -e HLAfreq
pip install -e .

Update version in setup.py

Update documentation with: `pdoc -d google -o docs/ HLAfreq`.
Note that github actions will automatically run this when pushed
to `main` branch.

Run tests `pytest`
Or allow nox to do it `nox`. Nox will also run linting.
On push github actions will run linting and pytest

Clear old build info
rm -rf build dist src/*.egg-info 

Build with `python -m build`.

twine check dist/*

Upload to test pypi
twine upload --repository testpypi dist/*

Install from test pypi
python3 -m pip install --extra-index-url https://test.pypi.org/simple/ HLAfreq

Upload to pypi
twine upload dist/*
-->

## Citation
Wells, D. A., & McAuley, M. (2023). HLAfreq: Download and combine HLA allele frequency data. bioRxiv, 2023-09. https://doi.org/10.1101/2023.09.15.557761 
