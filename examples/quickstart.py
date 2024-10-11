import HLAfreq
from HLAfreq import HLAfreq_pymc as HLAhdi
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Windows can have issues with multiprocessing in python
# main guard and spawn method help to execute everything in order
import multiprocessing
if __name__ == "__main__":
    multiprocessing.set_start_method('spawn')

    # Get HLA frequency data
    base_url = HLAfreq.makeURL("Uganda", locus="A")
    aftab = HLAfreq.getAFdata(base_url)

    # Preprocess data, checks etc
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)

    # Combine frequency estimates
    caf = HLAfreq.combineAF(aftab)
    # Calculate credible intervals
    hdi = HLAhdi.AFhdi(aftab, credible_interval=0.95)
    # Add credible intervals and posterior mean to caf
    caf = pd.merge(caf, hdi, how="left", on="allele")

    HLAfreq.plotAF(caf, aftab, hdi=hdi, compound_mean=hdi)
