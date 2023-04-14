import numpy as np
import matplotlib.pyplot as plt
import openmc
import openmc.lib
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



def model(k_eff_filepath, coolant='NaF', mod_density=1.65, fuel_temp=1000, coolant_temp=900, mod_temp=900, boron_equivalent_impurities=0.15, inactive=1, batches=30, particles=100, activate_tallies=False, simulate_photons=False, statepoint_filepath=""):
                                                                                                         # boron_equivalent_impurities given in ppm

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

    fuel_temperature = fuel_temp #(ref: 1000K)
    coolant_temperature = coolant_temp #(ref: 900K)
    mod_temperature = mod_temp #(ref: 900K)



    ##################################################################
    #                      DEFINE MATERIALS                          #
    ##################################################################

    ####################### TRISO MATERIALS ########################

    kernel = openmc.Material(name='kernel')
    kernel.set_density('g/cm3', 10.5)
    kernel.add_nuclide('U238', 0.7502, percent_type='wo')
    kernel.add_nuclide('U235', 0.1376, percent_type='wo')
    kernel.add_nuclide('O16',  0.0897, percent_type='wo')
    kernel.add_element('C',    0.0224, percent_type='wo')
    kernel.temperature = fuel_temperature

    buffer = openmc.Material(name='buffer')
    buffer.set_density('g/cm3', 1)
    buffer.add_element('C', 1.0, percent_type='wo')
    buffer.temperature = fuel_temperature

    iPyC = openmc.Material(name='iPyC')
    iPyC.set_density('g/cm3', 1.9)
    iPyC.add_element('C', 1.0, percent_type='wo')
    iPyC.temperature = fuel_temperature

    SiC = openmc.Material(name='SiC')
    SiC.set_density('g/cm3', 3.2)
    SiC.add_element('C', 0.5)
    SiC.add_element('Si', 0.5)
    SiC.add_s_alpha_beta('c_Si_in_SiC')
    SiC.add_s_alpha_beta('c_C_in_SiC')
    SiC.temperature = fuel_temperature

    oPyC = openmc.Material(name='oPyC')
    oPyC.set_density('g/cm3', 1.9)
    oPyC.add_element('C', 1.0, percent_type='wo')
    oPyC.temperature = fuel_temperature





    ####################### GRAPHITE ########################

    graphite = openmc.Material(name='graphite')
    graphite.set_density('g/cm3', mod_density)
    graphite.add_element('C',    0.99998985, percent_type='wo')
    graphite.add_nuclide('B10',  boron_equivalent_impurities*1e-6, percent_type='wo')
    graphite.add_nuclide('N14',  0.00001000, percent_type='wo')
    graphite.add_s_alpha_beta('c_Graphite')
    graphite.temperature = mod_temperature


    ####################### SALTS ########################

    NaF = openmc.Material(name='NaFZrF4')
    NaF.set_density('g/cm3', 2.96)
    NaF.add_nuclide('F19',  0.4537, percent_type='wo')
    NaF.add_nuclide('Na23', 0.1475, percent_type='wo')
    NaF.add_element('Zr',   0.3988, percent_type='wo')
    NaF.temperature = coolant_temp


    FLiBe_nat = openmc.Material(name='FLiBe_nat')
    FLiBe_nat.set_density('g/cm3', 1.925)
    FLiBe_nat.add_element('Li',  0.1414, percent_type='wo')
    FLiBe_nat.add_nuclide('Be9', 0.0909, percent_type='wo')
    FLiBe_nat.add_nuclide('F19', 0.7677, percent_type='wo')
    FLiBe_nat.temperature = coolant_temp


    FLiBe_enr1 = openmc.Material(name='FLiBe_enr1')
    FLiBe_enr1.set_density('g/cm3', 1.925)
    FLiBe_enr1.add_element('Li',  0.1414, percent_type='wo', enrichment=99.95, enrichment_target='Li7', enrichment_type='ao')
    FLiBe_enr1.add_nuclide('Be9', 0.0909, percent_type='wo')
    FLiBe_enr1.add_nuclide('F19', 0.7677, percent_type='wo')
    FLiBe_enr1.temperature = coolant_temp


    FLiBe_enr2 = openmc.Material(name='FLiBe_enr2')
    FLiBe_enr2.set_density('g/cm3', 1.925)
    FLiBe_enr2.add_element('Li',  0.1414, percent_type='wo', enrichment=99.995, enrichment_target='Li7', enrichment_type='ao')
    FLiBe_enr2.add_nuclide('Be9', 0.0909, percent_type='wo')
    FLiBe_enr2.add_nuclide('F19', 0.7677, percent_type='wo')
    FLiBe_enr2.temperature = coolant_temp

    if coolant=="NaF":
        coolant = NaF
    elif coolant=="FLiBe_nat":
        coolant = FLiBe_nat
    elif coolant=="FLiBe_enr1":
        coolant = FLiBe_enr1
    elif coolant=="FLiBe_enr2":
        coolant = FLiBe_enr2
    else:
        sys.exit("ERROR: Select an acceptable coolant!")



    ##################################################################
    #                       DEFINE GEOMETRY                          #
    ##################################################################

    top = openmc.ZPlane(z0=h/2, boundary_type='reflective')
    bottom = openmc.ZPlane(z0=-h/2, boundary_type='reflective')


    ####################### TRISO PARTICLE ########################
    spheres = [openmc.Sphere(r=1e-4*r)
            for r in [212.5, 312.5, 347.5, 382.5, 422.5]]
    cells = [openmc.Cell(fill=kernel,  region=-spheres[0]),
            openmc.Cell(fill=buffer,  region=+spheres[0] & -spheres[1]),
            openmc.Cell(fill=iPyC,    region=+spheres[1] & -spheres[2]),
            openmc.Cell(fill=SiC,     region=+spheres[2] & -spheres[3]),
            openmc.Cell(fill=oPyC,    region=+spheres[3] & -spheres[4])]
    triso_univ = openmc.Universe(cells=cells)


    ####################### FUEL COMPACTS ########################

    TRISO_cylinder = openmc.model.ZCylinder(r=fuel_compact_OR-0.02)
    compact_cylinder = openmc.model.ZCylinder(r=fuel_compact_OR)
    TRISO_top = openmc.ZPlane(z0=h/2-0.02, boundary_type='reflective')
    TRISO_bottom = openmc.ZPlane(z0=-h/2+0.02, boundary_type='reflective')

    TRISO_region = -TRISO_cylinder & -TRISO_top & +TRISO_bottom



    centers = openmc.model.pack_spheres(radius=spheres[4].r, region=TRISO_region, pf=0.35)
    trisos = [openmc.model.TRISO(spheres[4].r, triso_univ, center) for center in centers]

    print("Number of TRISO particles: {}".format(len(trisos)))

    compact_cell = openmc.Cell(fill=graphite, region= -compact_cylinder)
    compact_outer_cell = openmc.Cell(fill=graphite, region= +compact_cylinder)


    lower_left, upper_right = compact_cell.region.bounding_box
    lower_left[2] = -h/2
    upper_right[2] = h/2
    shape = (3, 3, 3)
    pitch = (upper_right - lower_left)/shape
    lattice = openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, graphite)
    compact_cell.fill = lattice

    compact_univ = openmc.Universe(cells=[compact_cell, compact_outer_cell])


    # if machine=="laptop":
    #     # front view
    #     compact_univ.plot(origin=(0.0, 0.0, 0.0), width=(2, 2), pixels=(600, 600), basis='xy', color_by='material')
    #     plt.savefig('images/compact.png')
    #     plt.clf()



    ####################### COOLANT CHANNEL ########################

    cc_cylinder = openmc.model.ZCylinder(r=coolant_channel_OR)
    cc_cell = openmc.Cell(fill=coolant, region= -cc_cylinder)
    cc_outer_cell = openmc.Cell(fill=graphite, region= +cc_cylinder)
    cc_univ = openmc.Universe(cells=[cc_cell, cc_outer_cell])

    # if machine=="laptop":
    #     # front view
    #     cc_univ.plot(origin=(0.0, 0.0, 0.0), width=(2, 2), pixels=(600, 600), basis='xy', color_by='material')
    #     plt.savefig('images/coolant_channel.png')
    #     plt.clf()



    ####################### ASSEMBLY ########################

    assembly_region = openmc.model.hexagonal_prism(edge_length=big_pitch, orientation='x', origin=(0.0, 0.0), boundary_type='periodic')

    assembly_inner_cell = openmc.Cell(fill=graphite, region=assembly_region)
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

    assembly_cell = openmc.Cell(fill=assembly_lattice, region=assembly_region & +bottom & -top)
    assembly_univ = openmc.Universe(cells=[assembly_cell])

    geom = openmc.Geometry(assembly_univ)
    geom.export_to_xml()
    mats = list(geom.get_all_materials().values())
    openmc.Materials(mats).export_to_xml()


    # if machine=="laptop":
    #     # front view
    #     assembly_univ.plot(origin=(0.0, 0.0, 0.0), width=(5.6388, 5.6388), pixels=(600, 600), basis='xy', color_by='material', colors={block_graphite: 'yellow', coolant:'green', compact_graphite:'purple'})
    #     plt.savefig('images/assembly_top.png')
    #     plt.clf()

    #     # side view
    #     assembly_univ.plot(origin=(0.0, 0.0, 0.0), width=(5.6388, 5.6388), pixels=(600, 600), basis='xz', color_by='material', colors={block_graphite: 'yellow', coolant:'green', compact_graphite:'purple'})
    #     plt.savefig('images/assembly_side.png')
    #     plt.clf()





    ###########################################################
    #                          TALLIES                        #
    ###########################################################

    if activate_tallies:
        #FILTERS
        #particle
        particle_filter = openmc.ParticleFilter('neutron')

        #material
        coolant_filter = openmc.MaterialFilter(coolant)
        graphite_filter = openmc.MaterialFilter(graphite)
        kernel_filter = openmc.MaterialFilter(kernel)

        #cells
        fuel_cell_filter = openmc.CellFilter(bins=2586)
        coolant_cell_filter = openmc.CellFilter(bins=6438)

        #energy
        energy_groups_fine = np.logspace(np.log(1e-3), np.log(20e6), base=np.e, num=300)
        energy_groups_fine[0] = 0
        energy_groups_2 = openmc.mgxs.GROUP_STRUCTURES["CASMO-2"]

        energy_filter_fine = openmc.EnergyFilter(energy_groups_fine)
        energy_filter_2 = openmc.EnergyFilter(energy_groups_2)
        energy_out_filter_2 = openmc.EnergyoutFilter(energy_groups_2)


        #TALLIES
        tallies = openmc.Tallies()

        # coolant flux spectrum
        t_flux_coolant = openmc.Tally(name='tally_flux_coolant')
        t_flux_coolant.filters = [particle_filter, coolant_filter, energy_filter_fine]
        t_flux_coolant.scores = ['flux']
        tallies.append(t_flux_coolant)

        # moderator flux spectrum
        t_flux_moderator = openmc.Tally(name='tally_flux_moderator')
        t_flux_moderator.filters = [particle_filter, graphite_filter, energy_filter_fine]
        t_flux_moderator.scores = ['flux']
        tallies.append(t_flux_moderator)

        # moderator flux spectrum
        t_flux_fuel = openmc.Tally(name='tally_flux_fuel')
        t_flux_fuel.filters = [particle_filter, kernel_filter, energy_filter_fine]
        t_flux_fuel.scores = ['flux']
        tallies.append(t_flux_fuel)


        # fast fission tally
        t_fast_fission = openmc.Tally(name='tally_fast_fission')
        t_fast_fission.filters = [particle_filter, energy_filter_2]
        t_fast_fission.scores = ['fission', 'nu-fission']
        tallies.append(t_fast_fission)

        # resonance escape tally
        t_res_escape = openmc.Tally(name='tally_res_escape')
        t_res_escape.filters = [particle_filter, energy_out_filter_2, energy_filter_2]
        t_res_escape.scores = ['scatter']
        tallies.append(t_res_escape)

        # heat deposition fuel
        t_heat_fuel = openmc.Tally(name='tally_heat_fuel')
        t_heat_fuel.filters = [kernel_filter]
        t_heat_fuel.scores = ['heating']
        tallies.append(t_heat_fuel)

        # heat deposition coolant
        t_heat_coolant = openmc.Tally(name='tally_heat_coolant')
        t_heat_coolant.filters = [coolant_filter]
        t_heat_coolant.scores = ['heating']
        tallies.append(t_heat_coolant)

        # heat deposition moderator
        t_heat_moderator = openmc.Tally(name='tally_heat_moderator')
        t_heat_moderator.filters = [graphite_filter]
        t_heat_moderator.scores = ['heating']
        tallies.append(t_heat_moderator)


        # power normalization tally
        t_power_normalization = openmc.Tally(name='tally_power_normalization')
        t_power_normalization.scores = ['heating', 'scatter']
        tallies.append(t_power_normalization)

        tallies.export_to_xml()


    ###########################################################
    #                         SETTINGS                        #
    ###########################################################

    settings = openmc.Settings()
    settings.inactive = inactive
    settings.batches = batches
    settings.particles = particles 
    settings.temperature={'method': 'interpolation'}
    settings.output = {
            'tallies': False,
            }

    # Shannon Entropy
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = [-big_pitch, -big_pitch, -h/2]
    entropy_mesh.upper_right = [big_pitch, big_pitch, h/2]
    entropy_mesh.dimension = (8, 8, 1)
    settings.entropy_mesh = entropy_mesh

    # Source initialization (point)
    point = openmc.stats.Point((-1.7, 0, 0))
    src = openmc.Source(space=point)
    settings.source = src
    settings.photon_transport = simulate_photons

    settings.export_to_xml()



    ##################################################################
    #                           RUN OPENMC                           #
    ##################################################################

    openmc.lib.init()
    openmc.lib.run()
    k_eff = openmc.lib.keff()
    openmc.lib.finalize()

    np.savetxt(k_eff_filepath, k_eff, delimiter=",")



    ##################################################################
    #                        ORDER FOLDER                            #
    ##################################################################

    os.remove("materials.xml")
    os.remove("geometry.xml")
    os.remove("settings.xml")

    if statepoint_filepath=="":
        os.remove("statepoint.{}.h5".format(batches))
    else:
        shutil.move("statepoint.{}.h5".format(batches), "{}".format(statepoint_filepath, batches))




    
