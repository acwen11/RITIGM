# This code was was initially developed by Leo Werneck and further edited by Tanmayee Gupte

# Import needed modules
import os
#import slicetable
from tabulated_eos_fcts import read_and_slice_EOS_table_at_given_temperature
import numpy as np
# Set table path and the slicing temperature
tablepath = os.path.join(".","SLy4_3335_rho391_temp163_ye66_adjusted.h5")
slicetemp = 0.01 # in MeV
# This command will provide the 1D interpolators we need, i.e.
# P(rho), eps(rho), and Y_{e}(rho) in beta-equilibrium at
# T = slicetemp
eos = read_and_slice_EOS_table_at_given_temperature(tablepath,slicetemp)
# Then output to file as needed
rho_code_units = eos.rho
P_code_units   = eos.P_of_rho_interpolator(rho_code_units)
eps_code_units = eos.eps_of_rho_interpolator(rho_code_units)
Ye  = eos.Ye_of_rho_interpolator(eos.rho)
T   = eos.Temperature*np.ones(len(eos.rho))
m_b = 1.639274e-24

#Getting proper units
rho_cgs = rho_code_units*eos.units_of_density
P_cgs   = P_code_units*eos.units_of_pressure
eps_cgs = eps_code_units*eos.units_of_spec_int_energy
dens_code_units = (1+eps_code_units)*rho_code_units
dens_cgs = ((1+eps_code_units)*rho_code_units)*eos.units_of_density
rho_fm3 = rho_cgs/(m_b*10**39)


#Getting the format ready to be read in Lorene
#First column is dummy index
#Second column is baryon number density in fm^-3 units
#Third column is total energy density in cgs
#fourth column is Pressure in cgs
x = []

for i in range(len(rho_code_units)):
    x.append(41)

x = np.array(x)

with open("./test.txt",'w') as outfile:
    print("Writing results on",outfile.name,"...")
    header= "#"
    no_of_lines = len(rho_code_units)
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    outfile.write(str(no_of_lines) + "\n")
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    outfile.write(header + "\n")
    fmt2 = '%d  ' + '%1.15e  ' + '%1.15e  ' + '%1.15e'
    DAT =  np.column_stack((x, rho_fm3, dens_cgs, P_cgs))
    np.savetxt('test.txt', DAT, header='#######\n#\n\n',fmt=fmt2) 
    

