"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA, UMAP, t-SNE). Then plot
those points coloured by geographic area.

For this each country must have an estimate for the same set
of alleles.
"""

# Use umap conda env

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import umap

AFeatures = pd.read_csv("data/example/globalPCA/AF_features.csv")
regions = pd.read_csv("data/example/countries.csv")

pca = PCA(n_components=2)
pca_embedding = pca.fit_transform(AFeatures.loc[:,'A*01:01':].values)
AFeatures['pca0'] = pca_embedding[:,0]
AFeatures['pca1'] = pca_embedding[:,1]

umapper = umap.UMAP(transform_seed=1234)
umap_embedding = umapper.fit_transform(AFeatures.loc[:,'A*01:01':].values)
AFeatures['umap0'] = umap_embedding[:,0]
AFeatures['umap1'] = umap_embedding[:,1]

# Add region to countries
AFeatures = pd.merge(AFeatures, regions, how="left", left_on="country", right_on="Country")

# Select Venezuela With Prior
vwp = AFeatures[AFeatures.country == "Venezuela with prior"]
# Drop countries without largeRegion, this is Venezuela with a prior
AFeatures = AFeatures.dropna(subset=['largeRegion'])

# plt.style.use('fivethirtyeight')

regions = AFeatures.largeRegion.unique()
regions.sort()
for region in regions:
    mask = AFeatures.largeRegion == region
    plt.scatter(AFeatures.pca0[mask], AFeatures.pca1[mask], label=region)
    for i in AFeatures[mask].index:
        plt.annotate(AFeatures.loc[i].country.replace("+", " "), AFeatures.loc[i][['pca0','pca1']])
plt.legend()
# Arrow shows change of Venezuela with prior
plt.annotate(
    "",
    xy=(vwp.pca0, vwp.pca1),
    xytext=(AFeatures[AFeatures.country == "Venezuela"].pca0, AFeatures[AFeatures.country == "Venezuela"].pca1),
    arrowprops=dict(arrowstyle="->")
    )
plt.xlabel('1st principal component'); plt.ylabel('2nd principal component')
plt.show()

for region in AFeatures.largeRegion.unique():
    mask = AFeatures.largeRegion == region
    plt.scatter(AFeatures.umap0[mask], AFeatures.umap1[mask], label=region)
    for i in AFeatures[mask].index:
        plt.annotate(AFeatures.loc[i].country.replace("+", " "), AFeatures.loc[i][['umap0','umap1']])
plt.legend()
plt.show()
