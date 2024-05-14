#!/usr/bin/env python
#
# EOS_Tables: manages EOS tables for THC
# Loosely based on Filippo Galeazzi's Matlab scripts
# Copyright (C) 2016, David Radice <dradice@caltech.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import h5py
import numpy as np
import sys
import os

#sys.path.append(os.path.dirname(sys.argv[0]) + "/../modules")
print(os.getcwd() + "/../modules")
sys.path.append(os.getcwd() + "/../modules")
import betaeq
from betaeq import myinterp1d
from eos_bar_table import COLWIDTH, EOS_Table
import unitconv as ut

INFO = "INFO".ljust(COLWIDTH)

# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# -----------------------------------------------------------------------------
parser.add_argument("-y", "--hydro", dest="hydro", required=True,
    help="Hydro part of the EOS table")
parser.add_argument("-w", "--weak", dest="weak", required=True,
    help="Weak part of the EOS table")
parser.add_argument("-o", "--output", dest="output", required=True,
    help="Output base file name (will be overwritten)")
parser.add_argument("-t", "--temperature", dest="temp", type=float,
    help="Slice temperature in MeV (defaults to the minimum on the table)")
parser.add_argument("-s", "--entropy", dest="entropy", type=float,
    help="Make slice at a given entropy as opposed to a fixed temperature")
parser.add_argument("--rm-radiation", dest="rm_radiation", action="store_true",
    help="Remove the radiation pressure part of the EOS")
parser.add_argument("--attach-poly", dest="attach_poly", action="store_true",
    help="Attach a polytrope at low density")
parser.add_argument("--resample", dest="resample", type=int, default=-1,
    help="Resample the EOS table, if negative don't resample (default: -1)")
parser.add_argument("--hdf5", dest="hdf5", action="store_true",
    help="Output the EOS slice in HDF5 format")
parser.add_argument("--ye-type", dest="ye_type", default="beta",
    help="Choose electron fraction (beta, given, or table)")
parser.add_argument("--ye-value", dest="ye_value", type=float,
    help="Electron fraction (only used if --ye=\"given\")")
parser.add_argument("--ye-table", dest="ye_table",
    help="Path to electron fraction table (only used if --ye=\"table\")")
parser.add_argument("-l", "--lorene", dest="lorene", action="store_true",
    help="Output the EOS slice in LORENE format")
parser.add_argument("-p", "--pizza", dest="pizza", action="store_true",
    help="Output the EOS slice in PIZZA format")
parser.add_argument("-r", "--rns", dest="rns", action="store_true",
    help="Output the EOS slice in RNS format")
parser.add_argument("--tovl", dest="tovl", action="store_true",
    help="Output the EOS slice in TOVL format")
parser.add_argument("--rho-min", dest="rho_min", type=float,
    help="Minimum density reached in the coldslice")
parser.add_argument("--rho-max", dest="rho_max", type=float,
    help="Maximum density reached in the coldslice")
args = parser.parse_args()

if args.temp is not None and args.entropy is not None:
    msg = "Only one of temperature or entropy can be specified!"
    raise argparse.ArgumentTypeError(msg)
if args.ye_type not in ["beta", "given", "table"]:
    msg = "Unknown electron fraction specification \"{}\"".format(args.ye_type)
    raise argparse.ArgumentTypeError(msg)
if args.ye_type == "given":
    if args.ye_value is None:
        msg = "Ye is set to \"given\", but no electron fraction has been specified"
        raise argparse.ArgumentTypeError(msg)
if args.ye_type == "table":
    if args.ye_table is None:
        msg = "Ye is set to \"table\", but no electron fraction table has been specified"
        raise argparse.ArgumentTypeError(msg)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Read the table into memory, compute quantities of interest
# -----------------------------------------------------------------------------
table = {}
print(INFO + "hydro table: {}".format(args.hydro))
dfile = h5py.File(args.hydro, "r")
table["rho"] = 10**np.array(dfile["logrho"])
table["temp"] = 10**np.array(dfile["logtemp"])
table["ye"]  = np.array(dfile["ye"])
table["cs2"] = np.array(dfile["cs2"])
table["eps"] = 10**np.array(dfile["logenergy"]) - dfile['energy_shift'][0]
table["press"] = 10**np.array(dfile["logpress"])
table["entropy"] = np.array(dfile["entropy"])
# Mass factor is the free neutron mass in the APR table
table["mass_factor"] = 939.6 #MeV
del dfile

print(INFO + "weak table: {}".format(args.weak))
dfile = h5py.File(args.weak, "r")
table["mu_e"] = np.array(dfile["mu_e"])
table["mu_p"] = np.array(dfile["mu_p"])
table["mu_n"] = np.array(dfile["mu_n"])
del dfile

if args.rho_min is not None:
  print(INFO + "rho min: {}".format(float(args.rho_min)))
  idx = table["rho"] > float(args.rho_min)
  table["rho"] = table["rho"][idx]
  for vname in table.keys():
    if vname not in ["rho", "temp", "ye", "mass_factor"]:
      table[vname] = table[vname][:,:,idx]

if args.rho_max is not None:
  print(INFO + "rho max: {}".format(float(args.rho_max)))
  idx = table["rho"] < float(args.rho_max)
  table["rho"] = table["rho"][idx]
  for vname in table.keys():
    if vname not in ["rho", "temp", "ye", "mass_factor"]:
      table[vname] = table[vname][:,:,idx]

mass_factor_cgs = table["mass_factor"]*(ut.MEV_CGS/(ut.C_CGS**2))
nb = table["rho"]/(mass_factor_cgs/(ut.FM_CGS**3))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Compute beta equilibrium slice
# -----------------------------------------------------------------------------
temp_slice = np.ones_like(table["rho"])
ye_slice = np.zeros_like(table["rho"])

if args.ye_type == "given":
    ye_slice[:] = float(args.ye_value)
    if args.entropy is not None:
        for inb in range(nb.shape[0]):
            t_of_s = betaeq.find_temp_given_ent(table["temp"], table["ye"],
                    table["entropy"][:,:,inb], args.entropy)
            temp_slice[inb] = myinterp1d(table["ye"], t_of_s)(ye_slice[inb])
        print(INFO + "making slice at Ye = {}, S = {} kb".format(ye_slice[0], args.entropy))
    else:
        if args.temp is not None:
            temp_slice[:] = float(args.temp)
        else:
            temp_slice[:] = table["temp"].min()
        print(INFO + "making slice at Ye = {}, T = {}".format(ye_slice[0], temp_slice[0]))
elif args.ye_type == "table":
    rho_ye_tab = np.loadtxt(args.ye_table, unpack=True)
    ye_slice = np.interp(np.log10(table["rho"]), rho_ye_tab[0], rho_ye_tab[1])
    if args.entropy is not None:
        for inb in range(nb.shape[0]):
            t_of_s = betaeq.find_temp_given_ent(table["temp"], table["ye"],
                    table["entropy"][:,:,inb], args.entropy)
            temp_slice[inb] = myinterp1d(table["ye"], t_of_s)(ye_slice[inb])
        print(INFO + "making slice at S = {} kb with tabulated Ye".format(args.entropy))
    else:
        if args.temp is not None:
            temp_slice[:] = float(args.temp)
        else:
            temp_slice[:] = table["temp"].min()
        print(INFO + "making slice at T = {} with tabulated Ye".format(temp_slice[0]))
else:
    if args.entropy is not None:
        print(INFO + "making beta equilibrium slice at S = {} kb".format(args.entropy))
        for inb in range(nb.shape[0]):
            t_of_s = betaeq.find_temp_given_ent(table["temp"], table["ye"],
                    table["entropy"][:,:,inb], args.entropy)
            mu_e_1d = np.zeros_like(table["ye"])
            mu_p_1d = np.zeros_like(table["ye"])
            mu_n_1d = np.zeros_like(table["ye"])
            for iye in range(mu_e_1d.shape[0]):
                mu_e_1d[iye] = myinterp1d(table["temp"], table["mu_e"][iye,:,inb])(t_of_s[iye])
                mu_n_1d[iye] = myinterp1d(table["temp"], table["mu_n"][iye,:,inb])(t_of_s[iye])
                mu_p_1d[iye] = myinterp1d(table["temp"], table["mu_p"][iye,:,inb])(t_of_s[iye])
            ye_slice[inb] = betaeq.find_beta_eq(table["ye"], mu_e_1d, mu_n_1d, mu_p_1d)
            temp_slice[inb] = myinterp1d(table["ye"], t_of_s)(ye_slice[inb])
    else:
        if args.temp is not None:
            mytemp = args.temp
        else:
            mytemp = table["temp"].min()
        print(INFO + "making beta-equilibrium slice at T = {} MeV".format(mytemp))
        temp_slice[:] = mytemp

        for inb in range(nb.shape[0]):
            mu_e_1d = myinterp1d(table["temp"], table["mu_e"][:,:,inb], axis=1)(temp_slice[inb])
            mu_n_1d = myinterp1d(table["temp"], table["mu_n"][:,:,inb], axis=1)(temp_slice[inb])
            mu_p_1d = myinterp1d(table["temp"], table["mu_p"][:,:,inb], axis=1)(temp_slice[inb])
            ye_slice[inb] = betaeq.find_beta_eq(table["ye"], mu_e_1d, mu_n_1d, mu_p_1d)
        assert np.all(np.isfinite(ye_slice))

rho0_beta = table["rho"]
cs2_beta = np.empty_like(ye_slice)
eps_beta = np.empty_like(ye_slice)
press_beta = np.empty_like(ye_slice)
for inb in range(nb.shape[0]):
    # Interpolate quantities to the right temperature
    cs2_1d = myinterp1d(table["temp"], table["cs2"][:,:,inb], axis=1)(temp_slice[inb])
    eps_1d = myinterp1d(table["temp"], table["eps"][:,:,inb], axis=1)(temp_slice[inb])
    press_1d = myinterp1d(table["temp"], table["press"][:,:,inb], axis=1)(temp_slice[inb])
    # Interpolate quantities to the right Ye
    cs2_beta[inb] = myinterp1d(table["ye"], cs2_1d)(ye_slice[inb])
    eps_beta[inb] = myinterp1d(table["ye"], eps_1d)(ye_slice[inb])
    press_beta[inb] = myinterp1d(table["ye"], press_1d)(ye_slice[inb])
assert np.all(np.isfinite(cs2_beta))
assert np.all(np.isfinite(eps_beta))
assert np.all(np.isfinite(press_beta))

# Clean soundspeed
cs2_beta[cs2_beta < 0.0] = 0.0
cs2_beta[cs2_beta > 1.0] = 1.0

# Total energy density in g/cc
rho_beta = nb * (mass_factor_cgs/(ut.FM_CGS**3)) * \
    (1.0 + (eps_beta/(ut.C_CGS**2)))
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Remove the radiation pressure / energy
# -----------------------------------------------------------------------------
if args.rm_radiation:
    print(INFO + "removing radiation pressure")
    rad_press = 1.0/3.0 * ut.RAD_CGS * (temp_slice * ut.MEV_CGS/ut.KB_CGS)**4
    rad_rho = rad_press * 3.0 / (ut.C_CGS**2)

    press_beta -= rad_press
    rho_beta -= rad_rho

    eps_beta = (rho_beta - rho0_beta)/rho0_beta * ut.C_CGS**2

    assert np.all(press_beta > 0)
    assert np.all(rho_beta > 0)
    assert np.all(eps_beta > 0)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Create barotropic EOS table object
# -----------------------------------------------------------------------------
mb_SI = (table["mass_factor"]*ut.MEV_SI)/(ut.C_SI**2)
mb_PU = mb_SI / ut.PIZZA_UNITS.mass
c_PU  = ut.C_SI / ut.PIZZA_UNITS.velocity
uc    = ut.CGS_UNITS / ut.PIZZA_UNITS
eos_slice = EOS_Table(
    rho0_beta * uc.density,
    eps_beta/(ut.C_CGS**2),
    press_beta * uc.pressure,
    efr  = ye_slice,
    temp = temp_slice,
    mbar = mb_PU,
    csnd = np.sqrt(cs2_beta)*c_PU,
    name = args.output)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Attach polytrope at low density and resample
# -----------------------------------------------------------------------------
if args.attach_poly:
    print(INFO + "attaching polytrope")
    eos_slice = eos_slice.make_restmass_natural(3)
    eos_slice = eos_slice.make_poly_compatible(3)
    if args.resample > 0:
        npoints = args.resample
    else:
        npoints = eos_slice.rmd.shape[0] + 50
    eos_slice = eos_slice.attach_poly(100*uc.density, npoints)
    eos_slice = eos_slice.make_adiabatic(remove_unphys_csnd=False)
    eos_slice.csnd = np.minimum(eos_slice.csnd, 0.999)
elif args.resample > 0:
    print(INFO + "resampling the EOS slice")
    eos_slice = eos_slice.resample_geom(args.resample)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Ensure that the pressure is monotonically increasing
# -----------------------------------------------------------------------------
eos_slice = eos_slice.remove_unphys_points()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Output the data in ASCII format
# -----------------------------------------------------------------------------
print(INFO + "exporting tables")
if args.hdf5:
    eos_slice.save("{}.hdf5".format(args.output))
if args.lorene:
    eos_slice.save_lorene("{}.lorene".format(args.output))
if args.pizza:
    eos_slice.save_pizza("{}.pizza".format(args.output))
if args.rns:
    eos_slice.save_rns("{}.rns".format(args.output))
if args.tovl:
    eos_slice.save_tovl("{}.dat".format(args.output))
# -----------------------------------------------------------------------------
