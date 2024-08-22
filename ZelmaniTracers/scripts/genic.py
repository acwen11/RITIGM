#!/usr/bin/env python

# Generates a GIZMO IC file from a ZelmaniTracer snapshot
# 1. Select which particles to keep based on different criteria
# 2. Apply symmetries

import h5py
import numpy
import optparse

from math import pi
from numpy import inf

# Useful routines for applying symmetries transformations
def double_copy_array(var):
    siz = var.shape[0]
    out = numpy.empty(2*siz)
    out[0:siz] = var
    out[siz:]  = var
    return out
def double_mirror_array(var):
    siz = var.shape[0]
    out = numpy.empty(2*siz)
    out[0:siz] = var
    out[siz:]  = -var
    return out

# Reflection symmetry across the equatorial plane
reflecting_z_scal = double_copy_array
def reflecting_z_vec(vec_x, vec_y, vec_z):
    out_x = double_copy_array(vec_x)
    out_y = double_copy_array(vec_y)
    out_z = double_mirror_array(vec_z)
    return out_x, out_y, out_z

# pi-symmetry
rotating_180_scal = double_copy_array
def rotating_180_vec(vec_x, vec_y, vec_z):
    out_x = double_mirror_array(vec_x)
    out_y = double_mirror_array(vec_y)
    out_z = double_copy_array(vec_z)
    return out_x, out_y, out_z

usage = "%prog [options] input output"
parser = optparse.OptionParser(usage)
parser.add_option("-b", "--bernoulli", dest="bernoulli", default=-inf,
        type="float", help="Remove particles with Bernoulli value" + \
                "smaller than a given value")
parser.add_option("-e", "--energy", dest="energy", default=-inf,
        type="float", help="Remove particles with specific kinetic" + \
                "energy at infinity smallar than a given value")
parser.add_option("--convert-temp", dest="conv_temp", default=False,
        action="store_true", help="Convert the temperature to Kelvin")
parser.add_option("--density-min", dest="rho_min", default=-inf,
        type="float", help="Enforce min density")
parser.add_option("--density-max", dest="rho_max", default=inf,
        type="float", help="Enforce max density")
parser.add_option("--mass-ratio", dest="mass_ratio", default=-inf, type="float",
        help="Select only particles with mass over max mass ratio " +\
                "larger than a given value")
parser.add_option("--reset-Abar", dest="Abar", default=None,
        type="float", help="Reset the composition to the given value of Abar")
parser.add_option("--reset-Neff", dest="Neff", default=None,
        type="float", help="Recompute the smoothing length from the " +\
                "given value of Neff")
parser.add_option("--rotating_180", dest="rotating_180", action="store_true",
        default=False, help="Apply pi-symmetry")
parser.add_option("--reflecting_z", dest="reflecting_z", action="store_true",
        default=False, help="Apply reflection symmetry across the z-plane")
parser.add_option("--temp-min", dest="temp_min", default=-inf,
        type="float", help="Enforce minimum temperature")
parser.add_option("--temp-max", dest="temp_max", default=inf,
        type="float", help="Enforce maximum temperature")
parser.add_option("--ye-min", dest="ye_min", default=-inf,
        type="float", help="Enforce minimum Ye")
parser.add_option("--ye-max", dest="ye_max", default=inf,
        type="float", help="Enforce maximum Ye")

options, args = parser.parse_args()

if len(args) != 2:
    parser.error("You must specify the input and output files")

# Tracer data filters
var_filters = {}

def filter_Density(data):
    rho = data.copy()
    rho[rho < options.rho_min] = options.rho_min
    rho[rho > options.rho_max] = options.rho_max
    return rho
var_filters['Density'] = filter_Density

def filter_Ye(data):
    ye = data.copy()
    ye[ye < options.ye_min] = options.ye_min
    ye[ye > options.ye_max] = options.ye_max
    return ye
var_filters['Ye'] = filter_Ye

def filter_Abar(data):
    if options.Abar is None:
        return data.copy()
    else:
        return options.Abar * numpy.ones(data.shape)
var_filters['Abar'] = filter_Abar

def filter_Temperature(data):
    if options.conv_temp:
        MeV = 1.16045e10 # Kelvin
        temp = data * MeV
    else:
        temp = data.copy()
    temp[temp < options.temp_min] = options.temp_min
    temp[temp > options.temp_max] = options.temp_max
    return temp
var_filters['Temperature'] = filter_Temperature

def filter_variable(name, data):
    if(var_filters.has_key(name)):
        return var_filters[name](data)
    else:
        return data.copy()

ifile = h5py.File(args[0],  'r')
ofile = h5py.File(args[1], 'w')

# First step: copy the header
ofile_header = ofile.create_group('/Header')
for hname in ifile['/Header'].attrs:
    ofile_header.attrs.create(hname, ifile['/Header'].attrs[hname])

# Second step: decide which particles to keep
eninf = numpy.array(ifile['/PartType0']['KineticEnergyAtInfinity'])
idx = eninf > options.energy

rho      = numpy.array(ifile['/PartType0']['Density'])
eps      = numpy.array(ifile['/PartType0']['InternalEnergy'])
mass     = numpy.array(ifile['/PartType0']['Masses'])
press    = numpy.array(ifile['/PartType0']['Pressure'])
enthalpy = 1.0 + eps + press/rho

B = enthalpy*(eninf + 1.0) - 1.0
idx = numpy.logical_and(idx, B > options.bernoulli)

mass_ratio = mass/mass.max()
idx = numpy.logical_and(idx, mass_ratio > options.mass_ratio)

siz = sum(idx)
if(options.reflecting_z):
    siz = siz*2
if(options.rotating_180):
    siz = siz*2
ofile_header.attrs['NumPart_ThisFile'] = numpy.array([siz, 0, 0, 0, 0, 0])
ofile_header.attrs['NumPart_Total']    = numpy.array([siz, 0, 0, 0, 0, 0])

# Third step: apply symmetries
ofile_ptype0 = ofile.create_group('/PartType0')

for name, dset in ifile['/PartType0'].iteritems():
    if name == "ParticleIDs":
        continue
    if options.Neff is not None and name == "SmoothingLength":
        continue
    old_data = numpy.array(dset)[idx]
    new_data = filter_variable(name, old_data)
    if len(new_data.shape) == 1:
        if options.reflecting_z:
            new_data = reflecting_z_scal(new_data)
        if options.rotating_180:
            new_data = rotating_180_scal(new_data)
    elif len(new_data.shape) == 2:
        data_x, data_y, data_z = tuple(new_data.transpose())
        if options.reflecting_z:
            data_x, data_y, data_z = reflecting_z_vec(data_x, data_y, data_z)
        if options.rotating_180:
            data_x, data_y, data_z = rotating_180_vec(data_x, data_y, data_z)
        new_data = numpy.array((data_x, data_y, data_z)).transpose()
    ofile_ptype0.create_dataset(name, data=new_data)
ofile_ptype0.create_dataset("ParticleIDs", data=numpy.arange(siz))

if options.Neff is not None:
    rho  = numpy.array(ofile_ptype0['Density'])
    mass = numpy.array(ofile_ptype0['Masses'])
    hsml = ((options.Neff * mass) / (4.0/3.0*pi*rho))**(1.0/3.0)
    ofile_ptype0.create_dataset("SmoothingLength", data=hsml)

ofile.close()
ifile.close()
