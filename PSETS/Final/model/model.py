import numpy as np
import matplotlib.pyplot as plt
import openmc
import openmc.model
import os
import sys
import shutil

# --------------------------------------------------------------------------- #
# Cross section data
cwd = os.getcwd()

if cwd[:6] == "/Users":
    os.environ['OPENMC_CROSS_SECTIONS'] = "/Users/lorenzomazzocco/openmc/endfb80_hdf5/cross_sections.xml"
    machine ="laptop"
if cwd[:13] == "/home/l_mazzo":
    os.environ['OPENMC_CROSS_SECTIONS'] = "/home/l_mazzo/nuclear_data/endfb80_hdf5/cross_sections.xml"
    machine ="cluster"
if cwd[:21] == "/home/lorenzomazzocco":
    os.environ['OPENMC_CROSS_SECTIONS'] = "/home/lorenzomazzocco/Desktop/nuclear_data/endfb80_hdf5/cross_sections.xml"
    machine ="office"

print(machine)





##################################################################
#                        MAIN PARAMETERS                         #
##################################################################

######################## GEOMETRY ########################

small_pitch = 1.8796
big_pitch = (3/2)*small_pitch

coolant_channel_OR = 0.7 / 2
fuel_compact_OR = 1.27 / 2

h = 2


######################## TEMPERATURES ########################

fuel_temperature = 1000 #K
other_temperature = 900 #K



##################################################################
#                      DEFINE MATERIALS                          #
##################################################################

####################### GRAPHITES ########################

block_graphite = openmc.Material(name='block_graphite')
block_graphite.set_density('g/cm3', 1.65)
block_graphite.add_element('C',    0.99998985, percent_type='wo')
block_graphite.add_nuclide('B10',  0.00000015, percent_type='wo')
block_graphite.add_nuclide('N14',  0.00001000, percent_type='wo')
block_graphite.add_s_alpha_beta('c_Graphite')
block_graphite.temperature = other_temperature


compact_graphite = openmc.Material(name='compact_graphite')
compact_graphite.set_density('g/cm3', 1.65)
compact_graphite.add_element('C',    0.99998985, percent_type='wo')
compact_graphite.add_nuclide('B10',  0.00000015, percent_type='wo')
compact_graphite.add_nuclide('N14',  0.00001000, percent_type='wo')
compact_graphite.add_s_alpha_beta('c_Graphite')
compact_graphite.temperature = fuel_temperature



####################### SALTS ########################

NaF = openmc.Material(name='NaFZrF4')
NaF.set_density('g/cm3', 2.96)
NaF.add_nuclide('F19',  0.4537, percent_type='wo')
NaF.add_nuclide('Na23', 0.1475, percent_type='wo')
NaF.add_element('Zr',   0.3988, percent_type='wo')
NaF.temperature = other_temperature


FLiBe_nat = openmc.Material(name='FLiBe_nat')
FLiBe_nat.set_density('g/cm3', 1.925)
FLiBe_nat.add_element('Li',  0.1414, percent_type='wo')
FLiBe_nat.add_nuclide('Be9', 0.0909, percent_type='wo')
FLiBe_nat.add_nuclide('F19', 0.7677, percent_type='wo')
FLiBe_nat.temperature = other_temperature


FLiBe_enr1 = openmc.Material(name='FLiBe_enr1')
FLiBe_enr1.set_density('g/cm3', 1.925)
FLiBe_enr1.add_element('Li',  0.1414, percent_type='wo', enrichment=99.95, enrichment_target='Li7', enrichment_type='ao')
FLiBe_enr1.add_nuclide('Be9', 0.0909, percent_type='wo')
FLiBe_enr1.add_nuclide('F19', 0.7677, percent_type='wo')
FLiBe_enr1.temperature = other_temperature


FLiBe_enr2 = openmc.Material(name='FLiBe_enr2')
FLiBe_enr2.set_density('g/cm3', 1.925)
FLiBe_enr2.add_element('Li',  0.1414, percent_type='wo', enrichment=99.995, enrichment_target='Li7', enrichment_type='ao')
FLiBe_enr2.add_nuclide('Be9', 0.0909, percent_type='wo')
FLiBe_enr2.add_nuclide('F19', 0.7677, percent_type='wo')
FLiBe_enr2.temperature = other_temperature

coolant = NaF



##################################################################
#                       DEFINE GEOMETRY                          #
##################################################################


####################### FUEL COMPACTS ########################

compact_cylinder = openmc.model.ZCylinder(r=fuel_compact_OR)
compact_cell = openmc.Cell(fill=compact_graphite, region= -compact_cylinder)
compact_outer_cell = openmc.Cell(fill=block_graphite, region= +compact_cylinder)
compact_univ = openmc.Universe(cells=[compact_cell, compact_outer_cell])

if machine=="laptop":
    # front view
    compact_univ.plot(origin=(0.0, 0.0, 0.0), width=(2, 2), pixels=(600, 600), basis='xy', color_by='material')
    plt.savefig('images/compact.png')
    plt.clf()



####################### COOLANT CHANNEL ########################

cc_cylinder = openmc.model.ZCylinder(r=coolant_channel_OR)
cc_cell = openmc.Cell(fill=coolant, region= -cc_cylinder)
cc_outer_cell = openmc.Cell(fill=block_graphite, region= +cc_cylinder)
cc_univ = openmc.Universe(cells=[cc_cell, cc_outer_cell])

if machine=="laptop":
    # front view
    cc_univ.plot(origin=(0.0, 0.0, 0.0), width=(2, 2), pixels=(600, 600), basis='xy', color_by='material')
    plt.savefig('images/coolant_channel.png')
    plt.clf()



####################### ASSEMBLY ########################

assembly_region = openmc.model.hexagonal_prism(edge_length=big_pitch, orientation='x', origin=(0.0, 0.0), boundary_type='periodic')

assembly_inner_cell = openmc.Cell(fill=block_graphite, region=assembly_region)
assembly_inner_univ = openmc.Universe(cells=[assembly_inner_cell])

a = compact_univ
b = cc_univ

assembly_lattice = openmc.HexLattice()
assembly_lattice.orientation='x'
assembly_lattice.center = (0, 0)
assembly_lattice.pitch = [small_pitch]
assembly_lattice.universes = \
[
[b]+[a]+[b]+[a]+[b]+[a], # 6 - Ring 1
[b] # 1 - Ring 0
]
assembly_lattice.outer = assembly_inner_univ

assembly_cell = openmc.Cell(fill=assembly_lattice, region=assembly_region)
assembly_univ = openmc.Universe(cells=[assembly_cell])


if machine=="laptop":
    # front view
    assembly_univ.plot(origin=(0.0, 0.0, 0.0), width=(5.6388, 5.6388), pixels=(600, 600), basis='xy', color_by='material', colors={block_graphite: 'yellow', coolant:'green', compact_graphite:'purple'})
    plt.savefig('images/assembly.png')
    plt.clf()