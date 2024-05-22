# Load necessary Python modules
import os
import sys                                   # Multiplatform OS specific functions
from bisect import bisect_left                  # Bisection algorithm
# Support for large, multi-dimensional arrays and matrices and mathematical functions
import h5py as h5                               # HDF5 file system support
from collections import namedtuple              # C-like struct functionality
# 1d and 2d interpolating functions
from scipy.interpolate import interp1d, interp2d
import bisect                                   # Bisection algorithms
import astropy.constants as constants
from numpy import empty, array, zeros, log, exp, loadtxt


def initialize_eos_struct(eos_table_path=os.path.join('IllinoisGRMHD', 'doc', 'eostable.h5'), idx_Ye=0, idx_T=0):
    # Unit conversion: cgs to dimensionless
    # Set physical constants
    G = constants.G.cgs.value     # cm^3/g/s^2
    c = constants.c.cgs.value     # cm/s
    Msun = constants.M_sun.cgs.value  # g

    # Set the length factor:
    # .----------------------.
    # | [Msun * G / c^2] = L |
    # .----------------------.
    units_of_length = Msun * G / c**2

    # Set the time factor:
    # .----------------------.
    # | [Msun * G / c^3] = T |
    # .----------------------.
    units_of_time = units_of_length/c

    # Pressure has units of M * L^(-1) * T^(-2), so
    units_of_pressure = Msun * units_of_length**(-1) * units_of_time**(-2)

    # Density has units of M * L^(-3), so
    units_of_density = Msun * units_of_length**(-3)

    # Specific internal energy has units of L^(2) * T^(-2),so
    units_of_spec_int_energy = c**2

    # Load the relevant quantities from the EOS table
    eos_table = h5.File(eos_table_path, 'r')
    tab_lp = log(10**array(eos_table["logpress"][:]))
    tab_lr = log(10**array(eos_table["logrho"][:]))
    tab_le = log(10**array(eos_table["logenergy"][:]))
    tab_lp_of_lr = tab_lp[idx_Ye, idx_T, :]
    tab_le_of_lr = tab_le[idx_Ye, idx_T, :]

    # Define the eos_struct named tuple
    eos_struct = namedtuple(
        "eos_struct", "tab_P_of_rho tab_eps_of_rho tab_rho")

    # Set the eos "struct" with dimensionless quantities
    eos = eos_struct(tab_P_of_rho/units_of_pressure,
                     tab_eps_of_rho/units_of_spec_int_energy,
                     tab_rho/units_of_density)

    return eos


def Read_EOS_table_and_sliced_table(eos_table_path, sliced_table_path):
    # Read in the EOS table
    eos_table = h5.File(eos_table_path, 'r')
    logrho = array(eos_table["logrho"])

    # Read in the sliced table values
    tempindex, rhoindex, yeindex, P, S, dummy, dummy, dummy, dummy, eps = loadtxt(
        os.path.join(sliced_table_path, "eos.thermo")).T
    yeb = loadtxt(os.path.join(sliced_table_path, "eos.yqb")).T

    # Convert P to cgs units
    # pressconv = 17.980714384531204 # math.log10(m_n * 1.0E-3 * 1.0E-4 / mev_si)  #useful constant for pressure conversion
    # P = 10**(log10(P) + logrho + pressconv)

    # Set rho
    rho = 10**logrho

    # Check that the length of all quantities match
    length_check = len(logrho)
    entries = [rhoindex, yeindex, P, rho, yeb]
    for entry in entries:
        if (len(entry) != length_check):
            print("Error! Incorrect length found!")
            sys.exit(1)
    print("Loaded EOS quantities successfully")

    # Convert from cgs to dimensionless units
    P, rho, eps = Convert_table_to_code_units(P, rho, eps)

    print("Converted P,rho,eps from cgs to dimensionless units")

    # Set the eos namedtuple
    eos_tuple = namedtuple("eos_tuple", "rho P yeindex yeb eps")
    eos = eos_tuple(rho, P, yeindex, yeb, eps)

    print("Finished without errors")

    return eos


def read_and_slice_EOS_table_at_given_temperature(eos_file_path, T_in_MeV):
    eos_file = h5.File(eos_file_path, "r")

    # Read in the ye, P, eps, T, and rho
    tab_Ye = array(eos_file['ye'])[:]
    tab_lp = log(10**array(eos_file['logpress'])[:])
    tab_le = log(10**array(eos_file['logenergy'])[:])
    tab_lt = log(10**array(eos_file['logtemp'])[:])
    tab_lr = log(10**array(eos_file['logrho'])[:])
    tab_ent = array(eos_file['entropy'])[:]

    # Read in the chemical potentials
#     munu     =     array(eos_file['munu'   ])[:]
    mu_e = array(eos_file['mu_e'])[:]
    mu_hat = array(eos_file['muhat'])[:]
    mu_l = mu_e - mu_hat

    # Find temperature index
    T_idx = bisect.bisect_left(tab_lt, log(T_in_MeV))
    print("Slicing table at index %d. T[%d] = %.3lf MeV | Input temperature: %.3lf MeV" % (
        T_idx, T_idx, exp(tab_lt[T_idx]), T_in_MeV))

    # Set 2d interpolators for P and eps
    # as functions of rho and Ye
    lp_interp = interp2d(tab_lr, tab_Ye, tab_lp[:, T_idx, :])
    le_interp = interp2d(tab_lr, tab_Ye, tab_le[:, T_idx, :])
    ent_interp = interp2d(tab_lr, tab_Ye, tab_ent[:, T_idx, :])

    # Now, loop over the density and set the array
    Ye_be = empty(len(tab_lr))
    lp_be = empty(len(tab_lr))
    le_be = empty(len(tab_lr))
    ent_be = empty(len(tab_lr))
    Y_e_min = 1.1*tab_Ye.min()
    for rho_idx in range(len(tab_lr)):
        # Set the interpolator
        #         Ye_ito_munu_of_T_and_rho = interp1d(munu[:,T_idx,rho_idx],tab_Ye)
        Ye_ito_munu_of_T_and_rho = interp1d(mu_l[:, T_idx, rho_idx], tab_Ye)
        # Then find the value of Ye which sets munu to zero
        try:
            YeL = Ye_ito_munu_of_T_and_rho(0.0)
        except:
            YeL = Y_e_min
        # Set Ye in beta equilibrium
        Ye_be[rho_idx] = YeL
        # Finally, determine P and eps in beta equilibrium
        lp_be[rho_idx] = lp_interp(tab_lr[rho_idx], Ye_be[rho_idx])
        le_be[rho_idx] = le_interp(tab_lr[rho_idx], Ye_be[rho_idx])
        ent_be[rho_idx] = ent_interp(tab_lr[rho_idx], Ye_be[rho_idx])

    # Convert table to code units
    # Set physical constants in cgs units
    G = 6.6738480e-8   # cm^3/g/s^2
    c = 2.99792458e10  # cm/s
    Msun = 1.9884e33      # g

    # Set the length factor:
    # .----------------------.
    # | [Msun * G / c^2] = L |
    # .----------------------.
    units_of_length = Msun * G / c**2

    # Set the time factor:
    # .----------------------.
    # | [Msun * G / c^3] = T |
    # .----------------------.
    units_of_time = units_of_length/c

    # Pressure has units of M * L^(-1) * T^(-2), so
    units_of_pressure = Msun * units_of_length**(-1) * units_of_time**(-2)

    # Density has units of M * L^(-3), so
    units_of_density = Msun * units_of_length**(-3)

    # Specific internal energy has units of L^(2) * T^(-2),so
    units_of_spec_int_energy = c**2

    # Convert the table to dimensionless units
    lp_be = log(exp(lp_be)/units_of_pressure)
    le_be = log(exp(le_be)/units_of_spec_int_energy)
    tab_lr = log(exp(tab_lr)/units_of_density)

    # Functions to interpolate P and eps
    lr_of_lp = interp1d(lp_be, tab_lr)
    lp_of_lr = interp1d(tab_lr, lp_be)
    le_of_lr = interp1d(tab_lr, le_be)
    ye_of_lr = interp1d(tab_lr, Ye_be)
    lent_of_lr = interp1d(tab_lr, log(ent_be))

    def rho_of_P_func(P):
        return exp(lr_of_lp(log(P)))

    def P_of_rho_func(rho):
        return exp(lp_of_lr(log(rho)))

    def eps_of_rho_func(rho):
        return exp(le_of_lr(log(rho)))

    def ent_of_rho_func(rho):
        return exp(lent_of_lr(log(rho)))

    def Ye_of_rho_func(rho):
        return ye_of_lr(log(rho))

    rho_of_P = rho_of_P_func
    P_of_rho = P_of_rho_func
    eps_of_rho = eps_of_rho_func
    Ye_of_rho = Ye_of_rho_func
    ent_of_rho = ent_of_rho_func

    # Set the EOS tuple
    eos_tuple = namedtuple(
        "eos_tuple", "rho Ye T Temperature rho_of_P_interpolator P_of_rho_interpolator eps_of_rho_interpolator Ye_of_rho_interpolator units_of_density units_of_pressure units_of_spec_int_energy ent_of_rho_interpolator")
    beta_equilibrated_eos_slice = eos_tuple(exp(tab_lr), tab_Ye, exp(tab_lt),
                                            exp(tab_lt[T_idx]),
                                            rho_of_P,
                                            P_of_rho,
                                            eps_of_rho,
                                            Ye_of_rho,
                                            units_of_density,
                                            units_of_pressure,
                                            units_of_spec_int_energy,
                                            ent_of_rho)

    return beta_equilibrated_eos_slice
