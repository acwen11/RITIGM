import os
import sys
import bisect
import h5py as h5
import numpy as np

# UNIT CONVERSIONS
## APR stores entropy in units of k_B / baryon, while helmholtz uses erg / K / g. Convert between the two using m_n / k_B.
ENTFAC = 1.2131450484808232e-08
K_TO_MEV = 8.61733326214518e-11 

def stitch(chi, tab, helm):
    return chi * tab + (1 - chi) * helm

def extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, orig_tab):
    # Assume composition does not change in low density, temperature region (see Hayashi et al 2023)
    new_tab = np.full(new_shape, np.nan)
    new_tab[:, num_new_pts_T:, num_new_pts_rho:] = orig_tab
    for iYe in range(new_shape[0]):
        for iT in range(new_shape[1]):
            for irho in range(new_shape[2]):
                if np.isnan(new_tab[iYe, iT, irho]):
                    irho_star = max(0, irho - num_new_pts_rho)
                    iT_star = max(0, iT - num_new_pts_T)
                    new_tab[iYe, iT, irho] = orig_tab[iYe, iT_star, irho_star]
    return new_tab

def extend_std_table_at_Ye(iYe, nT_new, nrho_new, irp, irm, iTp, iTm, chi_rho, chi_temp, extend_tab, orig_tab, h_tab):
    # In dens merge region
    extend_tab[iYe, nT_new:, irm:irp+1] = chi_rho * orig_tab[iYe, :, irm-nrho_new:irp-nrho_new+1] + (1 - chi_rho) * h_tab[nT_new:, irm:irp+1]

    # In low dens region
    extend_tab[iYe,nT_new:, :irm+1] = h_tab[nT_new:, :irm+1]

    # In temp merge region
    extend_tab[iYe, iTm:iTp+1, :] = chi_temp * extend_tab[iYe, iTm:iTp+1, :] + (1 - chi_temp) * h_tab[iTm:iTp+1, :]

    # low temp region
    extend_tab[iYe,:iTm+1,:] = h_tab[:iTm+1, :]

    return extend_tab

def main(original_tablepath):
    try:
        import helmholtz
    except:
        print("Did you install python-helmholtz?")
        return
 
    # Set New Table Parameters

    # Final Params
    logrho_min = -4
    logT_min = -6

    # TODO: Debug Params
    # logrho_min = 2.5
    # logT_min = -2.5

    # Create New Table
    output_tablepath   = original_tablepath.split(".h5")[0]+"_extended.h5"

    os.system("cp %s %s"%(original_tablepath,output_tablepath))

    print("Input  table file: ",original_tablepath)
    print("Output table file: ",output_tablepath)
    os.remove(output_tablepath)
    eos_file = h5.File(original_tablepath,"r")
    out_file = h5.File(output_tablepath, "a")

    # Open EOS file, read in values
    tab_Ye   = np.array(eos_file['ye'       ])[:]
    tab_logtemp = np.array(eos_file['logtemp'  ])[:]
    tab_logrho  = np.array(eos_file['logrho'   ])[:]

    # TODO: reset debug
    # tab_Ye   = np.array(eos_file['ye'       ])[:2]
    # tab_logtemp = np.array(eos_file['logtemp'  ])[:20]
    # tab_logrho  = np.array(eos_file['logrho'   ])[:20]

    tab_logrho_min = tab_logrho[0]
    tab_logT_min = tab_logtemp[0]

    dYe = np.abs(tab_Ye[1] - tab_Ye[0])
    dlogT = np.abs(tab_logtemp[1] - tab_logtemp[0])
    dlogrho = np.abs(tab_logrho[1] - tab_logrho[0])
    
    tab_P    = 10**np.array(eos_file['logpress' ])[:]
    tab_shift = eos_file['energy_shift'][0]
    tab_eps  = 10**np.array(eos_file['logenergy'])[:] - tab_shift 
    tab_S = np.array(eos_file['entropy'])[:]
    tab_cs2 = np.array(eos_file['cs2'])[:]
    tab_muhat = np.array(eos_file['muhat'])[:]
    tab_mu_n = np.array(eos_file['mu_n'])[:]
    tab_mu_p = np.array(eos_file['mu_p'])[:]
    tab_mu_e = np.array(eos_file['mu_e'])[:]
    tab_munu = np.array(eos_file['munu'])[:]
    tab_Xn = np.array(eos_file['Xn'])[:]
    tab_Xp = np.array(eos_file['Xp'])[:]
    tab_Xa = np.array(eos_file['Xa'])[:]
    tab_Xh = np.array(eos_file['Xh'])[:]
    tab_Abar = np.array(eos_file['Abar'])[:]
    tab_Zbar = np.array(eos_file['Zbar'])[:]
    tab_gamma = np.array(eos_file['gamma'])[:]

    # We want the new table to line up with the old one where they overlap
    num_new_pts_rho = int((tab_logrho[0] - logrho_min) / dlogrho)
    num_new_pts_T = int((tab_logtemp[0] - logT_min) / dlogT)
    rho_min_adj = tab_logrho[0] - num_new_pts_rho * dlogrho
    T_min_adj = tab_logtemp[0] - num_new_pts_T * dlogT
    total_pts_rho = num_new_pts_rho + len(tab_logrho)
    total_pts_T = num_new_pts_T + len(tab_logtemp)
    rho_space_final = np.linspace(rho_min_adj, tab_logrho[-1], total_pts_rho)
    T_space_final = np.linspace(T_min_adj, tab_logtemp[-1], total_pts_T)

    new_shape = (len(tab_Ye), len(T_space_final), len(rho_space_final))

    # STITCHING PARAMETERS
    rho_stitch = 2.5 # log10(rho in g cm^-3)
    T_stitch = -3.4  # log10(T in MeV)
    width = 0.25
    tol = 5e-2

    # TODO: reset
    # rho_stitch = 3.5 # log10(rho in g cm^-3)
    # T_stitch = -1.75  # log10(T in MeV)
    # width = 0.15
    # tol = 1e-1

    ## Transition region bounds. See MERGE in SROEOS User Guide
    rho_t_plus = rho_stitch - width * np.arctanh(2 * tol - 1)
    rho_t_minus =  rho_stitch + width * np.arctanh(2 * tol - 1)
    irho_plus = bisect.bisect_left(rho_space_final, rho_t_plus)
    irho_minus = bisect.bisect_left(rho_space_final, rho_t_minus)

    T_t_plus = T_stitch - width * np.arctanh(2 * tol - 1)
    T_t_minus =  T_stitch + width * np.arctanh(2 * tol - 1)
    iT_plus = bisect.bisect_left(T_space_final, T_t_plus)
    iT_minus = bisect.bisect_left(T_space_final, T_t_minus)

    rho_space_orig = np.arange(tab_logrho[0], tab_logrho[-1] + dlogrho, dlogrho)
    T_space_orig = np.arange(tab_logtemp[0], tab_logtemp[-1] + dlogT, dlogT)

    ## Stitch param check
    print("npts rho = {:d}; npts T = {:d}".format(len(rho_space_final), len(T_space_final)))
    print("rho_t_+ = {:.6f}; rho_t_- = {:.6f}; T_t_+ = {:.6f}; T_t_- = {:.6f}".format(rho_t_plus, rho_t_minus, T_t_plus, T_t_minus))
    print("irho_+ = {:d}; irho_- = {:d}; iT_+ = {:d}; iT_- = {:d}; ".format(irho_plus, irho_minus, iT_plus, iT_minus))
    print("new_pts_rho = {:d}; new_pts_T = {:d}".format(num_new_pts_rho, num_new_pts_T))
    print("Consistency check: rho[irho_+] = {:.4e}; rho[irho_-] = {:.4e}; T[iT_+] = {:.4e}; T[iT_-] = {:.4e}; ".format(rho_space_final[irho_plus], rho_space_final[irho_minus], T_space_final[iT_plus], T_space_final[iT_minus]))
    if (iT_minus < num_new_pts_T) or (irho_minus < num_new_pts_rho) or (iT_plus > len(T_space_final)) or (irho_plus > len(rho_space_final)):
        print("Error: transition region out of bounds.")
        return 
    

    # CREATE NEW COMPOSITION TABLES
    muhat = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_muhat)
    mu_n = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_n)
    mu_p = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_p)
    mu_e = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_mu_e)
    munu = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_munu)
    Xp = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xp)
    Xn = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xn)
    Xa = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xa)
    Xh = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Xh)
    Abar = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Abar)
    Zbar = extend_composition(new_shape, num_new_pts_T, num_new_pts_rho, tab_Zbar)

    # INITIALIZE EXTENDED TABLES
    P = np.empty(new_shape)
    S = np.empty(new_shape)
    cs2 = np.empty(new_shape)
    gamma = np.empty(new_shape)
    # eps needs special treatment due to eps_nuc (see Hayashi et al 2023)
    eps = np.full(new_shape, np.nan)

    # FILL ORIGINAL REGION
    P[:, num_new_pts_T:, num_new_pts_rho:] = tab_P
    eps[:, num_new_pts_T:, num_new_pts_rho:] = tab_eps
    S[:, num_new_pts_T:, num_new_pts_rho:] = tab_S
    cs2[:, num_new_pts_T:, num_new_pts_rho:] = tab_cs2
    gamma[:, num_new_pts_T:, num_new_pts_rho:] = tab_gamma

    # EXTEND REMAINING TABLES

    ## Set stitching functions
    chi_rho = 0.5 * (1 + np.tanh((rho_space_final[irho_minus:irho_plus+1] - rho_stitch) / width))
    chi_rho = np.tile(chi_rho, (new_shape[1] - num_new_pts_T, 1))

    chi_temp = 0.5 * (1 + np.tanh((T_space_final[iT_minus:iT_plus+1] - T_stitch) / width))
    chi_temp = np.tile(chi_temp, (new_shape[2], 1)).T

    ## Loop over Y_e
    for iYe, valYe in enumerate(tab_Ye):
        print("Calculating Ye = {:.4f}".format(valYe))

        ## Prepare input arrays for Helmholtz EOS
        helm_rhospace = 10**rho_space_final
        ## Helmholtz EOS has maximum rho * Y_e = 1e15. This will only be violated in the APREOS region or the forbidden high density, low temp region
        for irho in range(new_shape[2]):
            if helm_rhospace[irho] * valYe > 1e15:
                helm_rhospace[irho] = (1 - 1e-3) * 1e15 / valYe
        helm_rhospace = np.tile(helm_rhospace, (new_shape[1], 1))

        helm_Tspace = np.tile(10**T_space_final / K_TO_MEV, (new_shape[2], 1)).T

        # Create Helmholtz 2D table
        h = helmholtz.helmeos(dens=helm_rhospace, temp=helm_Tspace, abar=Abar[iYe,:,:], zbar=Abar[iYe,:,:]*valYe) # Y_e = zbar / abar 

        # Extend tables
        P = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, P, tab_P, h.ptot)
        S = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, S, tab_S, ENTFAC * h.stot)
        cs2 = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, cs2, tab_cs2, h.cs**2)
        gamma = extend_std_table_at_Ye(iYe, num_new_pts_T, num_new_pts_rho, irho_plus, irho_minus, iT_plus, iT_minus, chi_rho, chi_temp, gamma, tab_gamma, h.gam1)

        # eps needs special treatment due to eps_nuc (see Hayashi et al 2023)
        for iT in range(new_shape[1]):
            iT_star = max(0, iT - num_new_pts_T)
            for irho in range(new_shape[2]):
                irho_star = max(0, irho - num_new_pts_rho)
                if np.isnan(eps[iYe, iT, irho]):
                    eps_nuc = tab_eps[iYe, iT_star, irho_star] - h.etot[max(iT, num_new_pts_T), max(irho, num_new_pts_rho)]
                    eps[iYe, iT, irho] = h.etot[iT, irho] + eps_nuc

    ## End of Ye Loop

    # CALCULATE DERIVATIVES
    print("Calculating Derivatives")
    dT = 10**T_space_final
    drho = 10**rho_space_final

    dedT = np.gradient(eps, dT, axis=1)
    dPdT = np.gradient(P, dT, axis=1)
    dPde = dPdT / dedT
    dPdrho = np.gradient(P, drho, axis=2)
    dedrho = np.gradient(eps, drho, axis=2)
    dPdrhoe = dPdrho - dPdT * dedrho / dedT

    # Default to original values when possible, even if they are no longer accurate in the stitching regions. May enhance C2P? Needs to be tested.

    ## Load in original data
    tab_dedt = np.array(eos_file["dedt"])[:]
    tab_dPde = np.array(eos_file["dpderho"])[:]
    tab_dPdrhoe = np.array(eos_file["dpdrhoe"])[:]

    dedT[:,num_new_pts_T:,num_new_pts_rho:] = tab_dedt
    dPde[:,num_new_pts_T:,num_new_pts_rho:] = tab_dPde
    dPdrhoe[:,num_new_pts_T:,num_new_pts_rho:] = tab_dPdrhoe

    # Update new HDF5
    print("Outputting to file")
    out_file.create_dataset("pointsrho", data=np.array([len(rho_space_final)]), shape=(1,))
    out_file.create_dataset("pointstemp", data=np.array([len(T_space_final)]), shape=(1,))
    out_file.create_dataset("pointsye", data=np.array([len(tab_Ye)]), shape=(1,))
    out_file.create_dataset("logrho", data=rho_space_final)
    out_file.create_dataset("logtemp", data=T_space_final)
    out_file.create_dataset("ye", data=tab_Ye)
    out_file.create_dataset("logpress", data=np.log10(P))
    out_file.create_dataset("entropy", data=S)
    out_file.create_dataset("logenergy", data=np.log10(eps + tab_shift))
    out_file.create_dataset("energy_shift", data=np.array([tab_shift]), shape=(1,))
    out_file.create_dataset("muhat", data=muhat)
    out_file.create_dataset("mu_n", data=mu_n)
    out_file.create_dataset("mu_p", data=mu_p)
    out_file.create_dataset("mu_e", data=mu_e)
    out_file.create_dataset("munu", data=munu)
    out_file.create_dataset("dedt", data=dedT)
    out_file.create_dataset("dpdrhoe", data=dPdrhoe)
    out_file.create_dataset("dpderho", data=dPde)
    out_file.create_dataset("cs2", data=cs2)
    out_file.create_dataset("gamma", data=gamma)
    out_file.create_dataset("Xn", data=Xn)
    out_file.create_dataset("Xp", data=Xp)
    out_file.create_dataset("Xa", data=Xa)
    out_file.create_dataset("Xh", data=Xh)
    out_file.create_dataset("Abar", data=Abar)
    out_file.create_dataset("Zbar", data=Zbar)

    # Close files
    out_file.close()
    eos_file.close()

    return
 
if __name__ == '__main__':
    from sys import argv

    if len(argv) != 2:
        print(f"Correct usage: python3 {argv[0]} <Table Path>")
        exit(1)

    main(argv[1])
    print("Done :)")

