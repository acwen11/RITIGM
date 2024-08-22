#!/usr/bin/env python
#
# Plots the average neutrino energy times the lapse function.
# This should be constant in the Kerr-Schild metric if a = 0.

import matplotlib.pyplot as plt
import numpy as np
import scidata.carpet.hdf5 as h5
import os

OUTDIR = "alp_eave"
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)

fnames = ["alp.xz.h5", "rE[0].xz.h5", "rN[0].xz.h5"]

dset = h5.dataset(fnames)
grid = dset.get_reflevel()
x, z = grid.mesh()

for it in dset.iterations:
    print("Processing iteration {}...".format(it)),

    alp = dset.get_reflevel_data(grid, variable="ADMBASE::alp", iteration=it)
    rE = dset.get_reflevel_data(grid, variable="THC_M1::rE[0]", iteration=it)
    rN = dset.get_reflevel_data(grid, variable="THC_M1::rN[0]", iteration=it)

    alp_eave = alp*rE/rN
    alp_eave /= np.max(alp_eave[0,:])

    plt.contourf(x, z, alp_eave, levels=np.linspace(0.6, 1.4, 10), cmap='jet')
    plt.colorbar()
    plt.savefig(OUTDIR + "/" + str(it).zfill(4) + ".png")
    plt.close()

    print("done!")
