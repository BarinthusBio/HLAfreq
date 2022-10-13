"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA, UMAP, t-SNE). Then plot
those points coloured by geographic area.

For this each country must have an estimate for the same set
of alleles.
"""

import os
import code.HLAfreq as HLAfreq
import pandas as pd
import numpy as np

# Countries in regions as defined on 
# http://www.allelefrequencies.net/datasets.asp#tag_4
regions = pd.read_csv("data/example/countries.csv")

countries = regions.Country.tolist()

# Download HLA allele frequencies
# Not all countries have data
# those without will print "Failed to get data for..."
for country in countries:
    print()
    print(country)
    if not os.path.exists("data/example/globalPCA/%s_raw.csv" %country):
        base_url = HLAfreq.makeURL(country, standard="g", locus="A")
        try:
            aftab = HLAfreq.getAFdata(base_url)
            aftab.to_csv("data/example/globalPCA/%s_raw.csv" %country, index=False)
        except:
            print(f"Failed to get data for {country}: {base_url}")

cafs = []
for country in countries:
    try:
        df = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
        df = HLAfreq.only_complete(df)
        df = HLAfreq.decrease_resolution(df, 2)
        caf = HLAfreq.combineAF(df)
        caf['country'] = country
        cafs.append(caf)
    except:
        pass

cafs = pd.concat(cafs, axis=0).reset_index(drop=True)

# Give record for all alleles to all countries
cafs = HLAfreq.unmeasured_alleles(cafs, 'country')

# Check all countries have the same number of alleles
cafs.groupby('country').allele.unique().apply(len).unique()

# Sort by allele then select a single country,
# check that its alleles match a specified list
# then use the frequencies as features
cafs = cafs.sort_values('allele')

sorted_alleles = sorted(list(cafs.allele.unique()))

AFeatures = []
for country in countries:
    mask = cafs['country'] == country
    # Only generate features for countries with records
    # e.g. Thailand has no loci A records
    if any(mask):
        assert all(cafs[mask].allele.unique() == sorted_alleles)
        features = [country] + cafs[mask].allele_freq.tolist()
        AFeatures.append(features)

# Dataframe of allele frequencies read for dimension reduction
AFeatures = pd.DataFrame(AFeatures, columns = ['country'] + sorted_alleles)
AFeatures.to_csv("data/example/globalPCA/AF_features.csv", index=False)
