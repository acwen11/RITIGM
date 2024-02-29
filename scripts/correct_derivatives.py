import os
import sys
import bisect
import h5py as h5
import numpy as np

def main(extended_tablepath, ref_table):

    print("Input  table file: ",extended_tablepath)
    eos_file = h5.File(extended_tablepath,"r+")
    ref_file = h5.File(ref_table,"r")

    # Open EOS file, read in values
    tab_Ye   = np.array(eos_file['ye'       ])[:]
    tab_logtemp = np.array(eos_file['logtemp'  ])[:]
    tab_logrho  = np.array(eos_file['logrho'   ])[:]
    dlogrho = tab_logrho[1] - tab_logrho[0]
    dlogT = tab_logtemp[1] - tab_logtemp[0]

    tab_P = 10**np.array(eos_file['logpress'])[:]
    tab_shift = eos_file['energy_shift'][0]
    tab_eps  = 10**np.array(eos_file['logenergy'])[:] - tab_shift

    # Calculate Derivatives. Currently only setting what is used by WVU_EOS
    print("Calculating Derivatives")
    dT = 10**tab_logtemp
    drho = 10**tab_logrho

    dedT = np.gradient(tab_eps, dT, axis=1)
    dPdT = np.gradient(tab_P, dT, axis=1)
    dPde = dPdT / dedT
    dPdrho = np.gradient(tab_P, drho, axis=2)
    dedrho = np.gradient(tab_eps, drho, axis=2)
    dPdrhoe = dPdrho - dPdT * dedrho / dedT

    # Default to original values when possible. May enhance C2P
    nrho_orig = np.array(ref_file["pointsrho"])[0] 
    ntemp_orig = np.array(ref_file["pointstemp"])[0] 
    nrho_extend = np.array(eos_file["pointsrho"])[0] 
    ntemp_extend = np.array(eos_file["pointstemp"])[0] 
    num_new_rho = nrho_extend - nrho_orig
    num_new_temp = ntemp_extend - ntemp_orig

    ## Load in original data
    orig_dedt = np.array(ref_file["dedt"])[:]
    orig_dPde = np.array(ref_file["dpderho"])[:]
    orig_dPdrhoe = np.array(ref_file["dpdrhoe"])[:]

    dedT[:,num_new_temp:,num_new_rho:] = orig_dedt
    dPde[:,num_new_temp:,num_new_rho:] = orig_dPde
    dPdrhoe[:,num_new_temp:,num_new_rho:] = orig_dPdrhoe

    # Write to file
    eos_file["dedt"][...] = dedT
    eos_file["dpdrhoe"][...] = dPdrhoe
    eos_file["dpderho"][...] = dPde
    
    eos_file.close()

    return

if __name__ == '__main__':
    from sys import argv

    if len(argv) != 3:
        print(f"Correct usage: python3 {argv[0]} <Table Path> <Reference Table>")
        exit(1)
    
    main(argv[1], argv[2])
    print("Done :)")

