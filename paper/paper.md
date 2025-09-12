---
title: 'HLAfreq: Download and combine HLA allele frequency data'
tags:
  - Python
  - HLA
  - Human leukocyte antigen
  - allele frequency
  - population genetics
  - hierarchical Bayesian model
authors:
    - name: David A. Wells
      orcid: 0000-0002-4531-5968
      affiliation: '1'
    - name: Michael McAuley
      orcid: 0000-0003-2486-931X
      affiliation: '2'
affiliations:
    - name: Barinthus Biotherapeutics, United Kingdom
      index: 1
    - name: School of Mathematics and Statistics, Technological University Dublin, Dublin, Ireland
      index: 2
date: 06 August 2025
bibliography: paper.bib
---

# Summary
Human leukocyte antigen (HLA) genes encode cell-surface proteins which play an important role in immunity. Since different HLA alleles enable different immune responses, the population frequency of HLA alleles is often considered when designing vaccines [@gulukota1996hla]. Specific HLA alleles have been linked to autoimmune disease [@simmonds2007hla] and associated with adverse drug reactions [@fan2017hla]. Further, the success of solid organ and stem cell transplants is related to HLA matching between donor and recipient [@morishima2002clinical; furst2019hla].

The [Allele Frequency Net Database](www.allelefrequencies.net) is a publicly available repository for human immune gene frequency data from across the world [@Gonzalez-Galarza2020]. However, difficulties downloading and combining data from multiple studies make it hard for researchers to study larger regions or even single countries where the data is split across many sources. To address this gap, we present `HLAfreq`: a Python package which can be used to download, combine and analyse datasets from the Allele Frequency Net Database.

# Statement of need
The Allele Frequency Net Database is an excellent resource; however, downloading data from a large number of studies is currently manual and slow. After downloading multiple studies, combining them is hindered by different allele resolutions, missing alleles, and incomplete studies. `HLAfreq` provides functions to identify incomplete studies, handle missing alleles, harmonise allele resolution, calculate population coverage, and estimate allele frequencies and uncertainty using a Bayesian framework. Allele frequency plots can be generated to identify anomalous datasets and interesting diversity in a set of populations. To get started, see the guide and examples at [github.com/BarinthusBio/HLAfreq](https://github.com/BarinthusBio/HLAfreq).

The function `combineAF()` is used to combine multiple datasets and estimate the allele frequencies according to the 'default model' described below. So that larger studies contribute more to the combined allele frequency estimate, each dataset is weighted by twice the sample size by default (because each individual is diploid). Alternatively, any supplied weighting can be used; for example, population size can be used for a multi-country estimate, see the [multi-country example](https://barinthusbio.github.io/HLAfreq/HLAfreq/examples/multi_country.html). `combineAF()` performs several automated checks on the dataset. Studies are flagged as incomplete if the total allele frequency is outside a specified range (0.95-1.1 by default). Unmeasured alleles are added with a frequency of zero to ensure all populations report the same set of alleles before combining. Allele resolution is automatically checked because only alleles of the same resolution can be combined. The allele frequencies and credible intervals can be estimated according to the 'compound model' described below with `AFhdi()` imported from the submodule `HLAfreq.HLAfreq_pymc`.

# Methods

## Statistical methods
`HLAfreq` uses a Bayesian framework to estimate allele frequency statistics from combined datasets for a specific population. The user can select from two possible statistical models. The simpler 'default model' gives point estimates for allele frequencies. The more sophisticated 'compound model' gives both point estimates and credible intervals.

### Default model
Let $p_k$ be the frequency of the $k$-th allele of a particular gene in a given population (e.g. a country). The default model assumes that the observations from all datasets for the population are drawn independently and that the probability of being the $k$-th allele is $p_k$. In other words, each observation is drawn from a categorical distribution with parameters $(p_1,\dots,p_K)$ where $K$ is the total number of alleles. The prior for $(p_1,\dots,p_K)$ is taken to be a Dirichlet distribution with parameters $\alpha_1,\dots,\alpha_K$ which are chosen by the user. The Dirichlet distribution is a generalisation of the Beta distribution to higher dimensions; its properties are described in Section 4.6.3 of [@murphy2022probabilistic].

The Dirichlet distribution is conjugate to the categorical distribution, meaning that the posterior distribution for the default model is also Dirichlet. More precisely, if the combined datasets contain $x_k$ observations of the $k$-th allele (for $k=1,\dots,K$) then the posterior distribution is Dirichlet with parameters $\alpha_1+x_1,\dots,\alpha_K+x_K$. The posterior mean for the frequency of allele $j$ is then given by
$$
    \frac{\alpha_j+x_j}{\sum_{k=1}^K(\alpha_k+x_k)}.
$$

By default, `HLAfreq` takes the prior parameters to be $\alpha_1=\dots=\alpha_K=1$. This results in a uniform prior on $(p_1,\dots,p_K)$ subject to the constraints that $p_1,\dots,p_K\geq 0$ and $p_1+\dots+p_K=1$. The user can specify alternative values for $\alpha_1,\dots,\alpha_K$. These parameters may be interpreted as a `pseudocount' in the sense that choosing the prior $\alpha_1,\dots,\alpha_K$ is equivalent to taking a uniform prior and then observing a dataset with $\alpha_k-1$ observations of the $k$-th allele. (Intuitively the uniform prior corresponds to one observation of each allele). This can be used as a heuristic for choosing prior parameters based on external information.

The default model could in principle be used to estimate credible intervals for the frequencies $p_1,\dots,p_K$ however this option is not provided by `HLAfreq` because in practice such intervals are frequently unrealistically narrow. This is because the default model does not account for variance between studies. In the next subsection we outline how we can account for this variation and obtain accurate credible intervals. The current model is chosen as the default because it is simpler and we expect its point estimates to be sufficient for the majority of use cases.

### Compound model
The default model assumes that all observations are sampled from a homogeneous population; however, observations within a single study are more likely to be similar e.g. they may be sampled at the same time or place. To account for this, `HLAfreq` provides a 'compound model' which accounts for the grouping of observations within studies and allows the allele frequencies of study populations to differ from each other. The additional uncertainty of the model results in wider credible intervals which better reflect the true state of knowledge regarding population frequencies. This falls within the general class of hierarchical Bayesian models: see Chapter 5 [@gelman2014] for further details and background.

The compound model makes the following assumptions. As before, $p_k$ denotes the frequency of the $k$-th allele in the population and the prior distribution for $p_1,\dots,p_K$ is Dirichlet with parameters $\alpha_1,\dots,\alpha_K$. A concentration parameter $\gamma\geq 0$ is given with a standard log-normal prior distribution. For the $j$-th data source, a vector $\beta^{(j)}=(\beta^{(j)}_1,\dots,\beta^{(j)}_K)$ is sampled independently from a Dirichlet distribution with parameters $\gamma p_1,\dots,\gamma p_K$. Observations from the $j$-th data source are then sampled from a categorical distribution with parameters $\beta^{(j)}_1,\dots,\beta^{(j)}_K$. (Equivalently, the $j$-th data source as a whole is sampled from a multinomial distribution.)

Idiosyncratic sampling biases are captured by the different values of $\beta^{(j)}$, which result in different probabilities of sampling particular alleles for each data source. If $\gamma$ is large, then $\beta^{(j)}$ is likely to concentrate around $(p_1,\dots,p_K)$ which means that different studies tend to have similar allele frequencies.

The posterior distributions of $p_1,\dots,p_K$ and $\gamma$ do not have a closed form and so are estimated numerically using `PyMC` [@salvatier2016probabilistic]. The `HLAfreq` function `AFhdi` outputs posterior means and credible intervals for allele frequencies.

# Acknowledgements
MM was supported by the European Research Council (ERC) Advanced Grant QFPROBA (grant number 741487). DW is employed by Barinthus Biotherapeutics (UK) Ltd.

# References