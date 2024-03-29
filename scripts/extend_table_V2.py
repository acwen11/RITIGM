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

    # CREATE HELMHOLTZ TABLE
    ## Prepare input arrays for Helmholtz EOS
    helm_rhospace = np.tile(10**rho_space_final, (new_shape[1], 1))
    helm_rhospace = np.tile(helm_rhospace, (new_shape[0], 1, 1))

    helm_Tspace = np.tile(10**T_space_final / K_TO_MEV, (new_shape[2], 1)).T
    helm_Tspace = np.tile(helm_Tspace, (new_shape[0], 1, 1))

    helm_yespace = np.tile(tab_Ye, (new_shape[2], 1)).T
    helm_yespace = np.tile(helm_yespace, (new_shape[1], 1, 1))
    helm_yespace = np.swapaxes(helm_yespace, 1, 0)

    ## Helmholtz EOS has maximum rho * Y_e = 1e15. This will only be violated in the APREOS region or the forbidden high density, low temp region
    for iYe in range(new_shape[0]):
        for iT in range(new_shape[1]):
            for irho in range(new_shape[2]):
                if helm_rhospace[iYe, iT, irho] * helm_yespace[iYe, iT, irho] > 1e15:
                    helm_rhospace[iYe, iT, irho] = (1 - 1e-3) * 1e15 / helm_yespace[iYe, iT, irho]

    h = helmholtz.helmeos(dens=helm_rhospace, temp=helm_Tspace, abar=Abar, zbar=Abar*helm_yespace) # Y_e = zbar / abar 

    # Initialize temporary tables
    P = np.empty(new_shape)
    eps = np.empty(new_shape)
    S = np.empty(new_shape)
    cs2 = np.empty(new_shape)
    gamma = np.empty(new_shape)




    for iYe, valYe in enumerate(tab_Ye):
        print("Calculating Ye = {:.4f}".format(valYe))

        # Loop over temperature range covered by original table
        print("\tOriginal temp range loop.")
        for iT in range(len(T_space_final) - 1, num_new_pts_T - 1, -1):
            tab_iT = iT - num_new_pts_T
            tempK = 10**T_space_final[iT] / 8.61733326214518e-11 # convert logT to K
            #print("\tDebug: calc temp {:.4f} K".format(tempK))

            # Get Abar from APR table, assume mean mass number does not change in low density region (see Hayashi et al 2023)
            abar_approx = tab_Abar[iYe, tab_iT, 0]

            # Set Helmholtz slice where needed.
            h = helmholtz.helmeos(dens=10**rho_space_final[:irho_plus+1], temp=tempK, abar=abar_approx, zbar=abar_approx*valYe) # Y_e = zbar / abar 

            # In low density region, fill temp tables with Helmholtz values 
            for irho in range(0, irho_minus):
                irho_star = max(0, irho - num_new_pts_rho)

                P[iYe, iT, irho] = h.ptot[irho]
                ## Set eps offset due to nuclear binding energy in extended table region (see Hayashi et al 2023)
                eps_nuc = tab_eps[iYe, tab_iT, irho_star] - h.etot[max(irho, num_new_pts_rho)]
                eps[iYe, iT, irho] = h.etot[irho] + eps_nuc
                #print("\t\teps calc: eps extended = {:.4e}; eps Helm = {:.4e}; eps APR = {:.4e}; eps_nuc = {:.4e}".format(eps[iYe,iT,irho], h.etot[max(irho, num_new_pts_rho)], tab_eps[iYe, tab_iT, irho_star], eps_nuc))
                S[iYe, iT, irho] = ENTFAC * h.stot[irho]
                cs2[iYe, iT, irho] = h.cs[irho]**2
                gamma[iYe, iT, irho] = h.gam1[irho]

                # Assume composition is frozen in at low density
                muhat[iYe, iT, irho] = tab_muhat[iYe, tab_iT, irho_star]
                mu_n[iYe, iT, irho] = tab_mu_n[iYe, tab_iT, irho_star]
                mu_p[iYe, iT, irho] = tab_mu_p[iYe, tab_iT, irho_star]
                mu_e[iYe, iT, irho] = tab_mu_e[iYe, tab_iT, irho_star]
                munu[iYe, iT, irho] = tab_munu[iYe, tab_iT, irho_star]
                Xp[iYe, iT, irho] = tab_Xp[iYe, tab_iT, irho_star]
                Xn[iYe, iT, irho] = tab_Xn[iYe, tab_iT, irho_star]
                Xa[iYe, iT, irho] = tab_Xa[iYe, tab_iT, irho_star]
                Xh[iYe, iT, irho] = tab_Xh[iYe, tab_iT, irho_star]
                Abar[iYe, iT, irho] = tab_Abar[iYe, tab_iT, irho_star]
                Zbar[iYe, iT, irho] = tab_Zbar[iYe, tab_iT, irho_star]

            # In transition region, apply interpolation
            for irho in range(irho_minus, irho_plus):
                chi = 0.5 * (1 + np.tanh((rho_space_final[irho] - rho_stitch) / width))
                tab_irho = irho - num_new_pts_rho

                P[iYe, iT, irho] = stitch(chi, tab_P[iYe, tab_iT, tab_irho], h.ptot[irho])
                # Copy eps from table. eps_nuc treatment from above at least guarantees continuity
                eps[iYe, iT, irho] = tab_eps[iYe, tab_iT, tab_irho]
                S[iYe, iT, irho] = stitch(chi, tab_S[iYe, tab_iT, tab_irho], ENTFAC * h.stot[irho])
                cs2[iYe, iT, irho] = stitch(chi, tab_cs2[iYe, tab_iT, tab_irho], h.cs[irho]**2)
                gamma[iYe, iT, irho] = stitch(chi, tab_gamma[iYe, tab_iT, tab_irho], h.gam1[irho])

                # Assume composition is frozen in at low density
                muhat[iYe, iT, irho] = tab_muhat[iYe, tab_iT, tab_irho]
                mu_n[iYe, iT, irho] = tab_mu_n[iYe, tab_iT, tab_irho]
                mu_p[iYe, iT, irho] = tab_mu_p[iYe, tab_iT, tab_irho]
                mu_e[iYe, iT, irho] = tab_mu_e[iYe, tab_iT, tab_irho]
                munu[iYe, iT, irho] = tab_munu[iYe, tab_iT, tab_irho]
                Xp[iYe, iT, irho] = tab_Xp[iYe, tab_iT, tab_irho]
                Xn[iYe, iT, irho] = tab_Xn[iYe, tab_iT, tab_irho]
                Xa[iYe, iT, irho] = tab_Xa[iYe, tab_iT, tab_irho]
                Xh[iYe, iT, irho] = tab_Xh[iYe, tab_iT, tab_irho]
                Abar[iYe, iT, irho] = tab_Abar[iYe, tab_iT, tab_irho]
                Zbar[iYe, iT, irho] = tab_Zbar[iYe, tab_iT, tab_irho]

            # In high density region, copy from original table
            for irho in range(irho_plus, len(rho_space_final)):
                tab_irho = irho - num_new_pts_rho

                P[iYe, iT, irho] = tab_P[iYe, tab_iT, tab_irho]
                eps[iYe, iT, irho] = tab_eps[iYe, tab_iT, tab_irho]
                S[iYe, iT, irho] = tab_S[iYe, tab_iT, tab_irho]
                cs2[iYe, iT, irho] = tab_cs2[iYe, tab_iT, tab_irho]
                gamma[iYe, iT, irho] = tab_gamma[iYe, tab_iT, tab_irho]
                muhat[iYe, iT, irho] = tab_muhat[iYe, tab_iT, tab_irho]
                mu_n[iYe, iT, irho] = tab_mu_n[iYe, tab_iT, tab_irho]
                mu_p[iYe, iT, irho] = tab_mu_p[iYe, tab_iT, tab_irho]
                mu_e[iYe, iT, irho] = tab_mu_e[iYe, tab_iT, tab_irho]
                munu[iYe, iT, irho] = tab_munu[iYe, tab_iT, tab_irho]
                Xp[iYe, iT, irho] = tab_Xp[iYe, tab_iT, tab_irho]
                Xn[iYe, iT, irho] = tab_Xn[iYe, tab_iT, tab_irho]
                Xa[iYe, iT, irho] = tab_Xa[iYe, tab_iT, tab_irho]
                Xh[iYe, iT, irho] = tab_Xh[iYe, tab_iT, tab_irho]
                Abar[iYe, iT, irho] = tab_Abar[iYe, tab_iT, tab_irho]
                Zbar[iYe, iT, irho] = tab_Zbar[iYe, tab_iT, tab_irho]
            del h

        # Loop over extended temperature range 
        print("\tExtended temp range loop.")

        ## Update points in transition region
        for iT in range(iT_plus, iT_minus, -1):
            tempK = 10**T_space_final[iT] / 8.61733326214518e-11 # convert logT to K
            #print("\tDebug: calc temp {:.4f} K".format(tempK))

            # From looking at the APR table, Abar ~ 74 at low temperatures and densities. This is not true at high densities but we hope that region of table space is not reached.
            h = helmholtz.helmeos(dens=10**rho_space_final[:helm_irho_max+1], temp=tempK, abar=74, zbar=74*valYe) # Y_e = zbar / abar 

            chi = 0.5 * (1 + np.tanh((T_space_final[iT] - T_stitch) / width))
            for irho in range(len(rho_space_final)):
                # Enforce Density limit on Helmholtz EOS. These values should not be used in a simulation.
                helm_irho_cap = min(helm_irho_max, irho)

                P[iYe, iT, irho] = stitch(chi, P[iYe, iT, irho], h.ptot[helm_irho_cap])
                # Don't update eps. eps_nuc treatment at least guarantees continuity
                S[iYe, iT, irho] = stitch(chi, S[iYe, iT, irho], ENTFAC * h.stot[helm_irho_cap])
                cs2[iYe, iT, irho] = stitch(chi, cs2[iYe, iT, irho], h.cs[helm_irho_cap]**2)
                gamma[iYe, iT, irho] = stitch(chi, gamma[iYe, iT, irho], h.gam1[helm_irho_cap])
            del h
        
        ## Copy values from Helm EOS in low temp range
        for iT in range(iT_minus, -1, -1):
            tempK = 10**T_space_final[iT] / 8.61733326214518e-11 # convert logT to K
            #print("\tDebug: calc temp {:.4f} K".format(tempK))
            h = helmholtz.helmeos(dens=10**rho_space_final[:helm_irho_max+1], temp=tempK, abar=74, zbar=74*valYe) # Y_e = zbar / abar 
            iT_star = max(0, iT - num_new_pts_T)
            for irho in range(len(rho_space_final)):
                helm_irho_cap = min(helm_irho_max, irho)

                P[iYe, iT, irho] = h.ptot[helm_irho_cap]
                ## Set eps offset due to nuclear binding energy in extended table region (see Hayashi et al 2023)
                irho_star = max(0, irho - num_new_pts_rho)
                eps_nuc = tab_eps[iYe, iT_star, irho_star] - h.etot[max(helm_irho_cap, num_new_pts_rho)]
                eps[iYe, iT, irho] = h.etot[helm_irho_cap] + eps_nuc
                S[iYe, iT, irho] = ENTFAC * h.stot[helm_irho_cap]
                cs2[iYe, iT, irho] = h.cs[helm_irho_cap]**2
                gamma[iYe, iT, irho] = h.gam1[helm_irho_cap]

                # Assume composition is frozen in at low density, temperature
                muhat[iYe, iT, irho] = tab_muhat[iYe, iT_star, irho_star]
                mu_n[iYe, iT, irho] = tab_mu_n[iYe, iT_star, irho_star]
                mu_p[iYe, iT, irho] = tab_mu_p[iYe, iT_star, irho_star]
                mu_e[iYe, iT, irho] = tab_mu_e[iYe, iT_star, irho_star]
                munu[iYe, iT, irho] = tab_munu[iYe, iT_star, irho_star]
                Xp[iYe, iT, irho] = tab_Xp[iYe, iT_star, irho_star]
                Xn[iYe, iT, irho] = tab_Xn[iYe, iT_star, irho_star]
                Xa[iYe, iT, irho] = tab_Xa[iYe, iT_star, irho_star]
                Xh[iYe, iT, irho] = tab_Xh[iYe, iT_star, irho_star]
                Abar[iYe, iT, irho] = tab_Abar[iYe, iT_star, irho_star]
                Zbar[iYe, iT, irho] = tab_Zbar[iYe, iT_star, irho_star]
            del h

        print("Done with Ye = {:.4f}".format(valYe))

    # Calculate Derivatives.
    print("Calculating Derivatives")
    dT = 10**tab_logtemp
    drho = 10**tab_logrho

    dedT = np.gradient(tab_eps, dT, axis=1)
    dPdT = np.gradient(tab_P, dT, axis=1)
    dPde = dPdT / dedT
    dPdrho = np.gradient(tab_P, drho, axis=2)
    dedrho = np.gradient(tab_eps, drho, axis=2)
    dPdrhoe = dPdrho - dPdT * dedrho / dedT

    # Default to original values when possible, even if they are no longer accurate in the stitching regions. May enhance C2P

    ## Load in original data
    orig_dedt = np.array(ref_file["dedt"])[:]
    orig_dPde = np.array(ref_file["dpderho"])[:]
    orig_dPdrhoe = np.array(ref_file["dpdrhoe"])[:]

    dedT[:,num_new_pts_temp:,num_new_pts_rho:] = orig_dedt
    dPde[:,num_new_pts_temp:,num_new_pts_rho:] = orig_dPde
    dPdrhoe[:,num_new_pts_temp:,num_new_pts_rho:] = orig_dPdrhoe

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

