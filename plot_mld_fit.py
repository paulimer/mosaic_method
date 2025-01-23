#!/usr/bin/env python3
"Plots the fit of the mosaic model to the MLDs. Translation of plot_mld_fit.R."

import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from fit.fit import *

def plot_mld_fit(binned_mld, muc, mus, delta, fitted_params, level, L0, outfile=None):
    """Plots the fit of the mosaic model to the MLDs."""
    max_x = binned_mld["match_length"].max()
    min_y = binned_mld[binned_mld["freq"] > 0]["freq"].min()
    print(f"min_y : {min_y}")
    if binned_mld["freq"].sum() == 0:
        # mld is empty, do not plot
        return
    r = np.logspace(0, np.log10(max_x), 1000)
    mh, mc = theoretical_mld(fitted_params, 0.1, r, mus, muc, delta, L0, False)

    fig, ax = plt.subplots()
    ax.plot(binned_mld["match_length"], binned_mld["freq"], 'o', label="Observed", color="black")
    ax.plot(r, mh, label="mh", color="red")
    ax.plot(r, mc, label="mc", color="blue")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(min_y/10, None)
    ax.legend()
    ax.set_title(f"MLD fit for {level[0]} vs {level[1]}")
    ax.set_xlabel("Match length")
    ax.set_ylabel("Frequency")
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()
