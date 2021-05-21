"""
    This is an example script to automatically generate TOUGH input files from
    GEOS outputs. It performs 5 tasks:
        1 - Call MeshMaker to build TOUGH mesh
        2 - Insert GEOS outputs into the fracture cells
        3 - Calculate heterogeneous fracture permeability
        4 - Build secondary fracture (SF) network between hydraulic fractures (HF)
        5 - Write INPUT file for TOUGH

    There are several requirements on the test problem to be set up:
        1 - The domain must be rectangle with Cartesian grids
        2 - The HF must be planar and it must spread in the x-z plane
        3 - The horizontal well must be along the y-axis (perpendicular to HF)
        4 - The numerical code must be pTOUGH+OGB with aqueous, oil and gas phases

    Author : Zhi Li, Lawrence Berkeley National Laboratory, 05-10-2021
"""
import numpy as np
from mm_input import *
import subprocess as sp

"""
    Top-level options
"""
#   Embeds GEOS outputs or not
use_geos = True
#   Uses MINC or not
use_minc = True
#   Adds secondary fractures or not
use_secf = True
#   Use closed aperture or not
use_clos = False
#   Builds fracture-conforming mesh or not
remesh = True

"""
    Information of the present test problem
"""
#   Simulation ID
#   This ID allows the user to write problem-specific settings. To build a new
#   test problem, the user should define a new ID
sim_id = '09L'
#   Fracture ID
#   This ID determines which fracture (among the 5 fractures in a typical GEOS
#   simulation) to be modeled. It ranges from 1 to 5.
ifrac = 2
#   Mean and standard deviation of the target aperture field [m]
#   If these values are not None, the GEOS aperture will be scaled to obtain
#   certain mean and standard deviation. For example, to simulate a homogeneous
#   aperture field, the user should set std_aper = 0
avg_aper = None
std_aper = None
#   Blockage ID
#   This ID allows the user to set up specific instructions to remove isolated
#   fracture regions for a given fracture.
blk_key = 'up09l2'
#   Default domain dimensions
#   NX, NY, NZ of the model domain.
dim = [117, 35, 45]
#   Default grid resolutions [m]
#   dx, dy, dz of the model domain.
delta = [6.0, 0.002, 4.0]
#   y-coordinate of the HF plane [m]
yfrac = 8.0
#   Time info of the simulation [sec], [simulation length, initial dt, max dt]
timestep = [1.296e8, 1.0, 8.64e5]
#   Time to save model outputs [sec]
t_save = [8.64e4, 8.64e5, 8.64e6, 8.64e7, 1.296e8]

"""
    Information on the media
"""
#   Information of the horizontal well
#   Note : This also shows how some options can be modified for specific test
#           problems by changing sim_id.
if sim_id == '09L':
    boundary = [1, ['HWELL', [354.0, 360.0, 0.0, 16.0, -92.0, -88.0,
                    ['properties', 2.6e3, 1.0, 1e-4, 0e-8, 6, 5]]]
                ]
    dim = [117, 35, 45]
elif sim_id == '09H' or sim_id == 'SMH':
    boundary = [1, ['HWELL', [354.0, 360.0, 0.0, 16.0, -252.0, -248.0,
                    ['properties', 2.6e3, 1.0, 1e-4, 0e-8, 6, 5]]]
                ]
    dim = [117, 35, 85]
#   Medium properties
#   These information will be written to the ROCKS block of TOUGH
#   The structure of 'region' is as follows:
#   region = [number of regions,
#            [region name, [x_min, x_max, y_min, y_max, z_min, z_max],
#                    ['properties', density, porosity, permeability, compressibility, relperm func, cap func],
#                    ['dualporosity', density, porosity, permeability, compressibility, relperm func, cap func]],
#            [the next region...]   ]
#   Note : To use MINC for a medium, the 'dualporosity' row should be present and
#           the region name should starts with 'F'.
region = [4, ['MATRX', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        ['properties', 2.6e3, 0.07, 1e-21, 1e-8, 8, 8],
                        ['dualporosity', 2.6e3, 0.07, 1e-21, 1e-8, 8, 8]]
             ],
            ['FSRV1', [20.0, 680.0, 1.0, yfrac-0.5*delta[1], -(dim[2]-1)*delta[2], -delta[2],
                        ['properties', 2.6e3, 0.07, 1e-21, 1e-8, 8, 8],
                        ['dualporosity', 2.6e3, 0.9, 1e-21, 1e-8, 8, 8]]
            ],
            ['FSRV2',   [20.0, 680.0, yfrac+0.5*delta[1], 15.0, -(dim[2]-1)*delta[2], -delta[2],
                        ['properties', 2.6e3, 0.07, 1e-21, 1e-8, 8, 8],
                        ['dualporosity', 2.6e3, 0.9, 1e-21, 1e-8, 8, 8]]
            ],
            ['MHFRA', [20.0, 680.0, yfrac-0.5*delta[1], yfrac+0.5*delta[1], -(dim[2]-1)*delta[2], -delta[2],
                        ['properties', 2.6e3, 0.07, 1e-21, 1e-8, -8, -8]]
            ]
          ]
#   Information on MINC regions
#   minc = [number_of_media, number_of_fractured_media, number_of_MINCs, number_of_specified_VolFractions, matrix_to_matrix_flow,
#      [fractured_medium_number, [?, ?, ?, ?], volume_fraction],
#      [the next MINC zone...] ]
minc = [4, 2, 2, 1, 'A',
            [2, [5.0, 5.0, 5.0, 5.0], 0.01],
            [3, [5.0, 5.0, 5.0, 5.0], 0.01]]
#   The medium ID in MESH for the HF and the fractures in the SRV
#   When using MINC, the ID must be numbers.
if use_minc:
    frac_zone = ['    7']
    srv_frac_zone = ['    3', '    5']
else:
    frac_zone = ['MHFRA']
    srv_frac_zone = ['FSRV1', 'FSRV2']
#   Information on the secondary fractures
#   The secondary fractures are assumed parallel to the HF with permeability kx = kz != ky.
#   secf = [[origin_x, origin_y], [L, W], max_aper, min_aper, kx/ky, perm_cutoff]
#       origin = Location with the largest aperture. Often set to the intersection
#               between HF and the well.
#       [L, W] = Length, width of the secf zone in [m].
#       max_aper = Maximum aperture near the well [m].
#       min_aper = Minimum aperture at L from the well [m].
#       kx/ky = Anisotropy ratio for [HF, SRV].
#       min_perm = Minimum permeability a secondary fracture can get [m^2].
secf = [[357.0, 8.0], [100.0, 7.0], 1e-7, 1e-10, [1e-21, 1000.0], 1e-21]

"""
    Information on domain discretization
"""
#   discrete = {
#               'disc_x':   [number_of_zones, [number_of_grids_in_this_zone, mode, grid_resolution], [the next zone...]]
#               }
#   'Mode' can be 'Equal', 'Uniform' or 'Log'. This follows the available options
#   in the new MeshMaker.
discret = {
        'disc_x':   [1, [dim[0], 'Equal', delta[0]]],
        'disc_y':   [9, [1, 'Uniform', 1.0], [14, 'Log', 1.0, 7.9, 1.0],
                        [1, 'Uniform', 0.078], [1, 'Uniform', 0.021+0.5*(0.002-delta[1])],
                        [1, 'Uniform', delta[1]],
                        [1, 'Uniform', 0.021+0.5*(0.002-delta[1])], [1, 'Uniform', 0.078],
                        [14, 'Log', 1.0, 15.0, 0.1865], [1, 'Uniform', 1.0]],
        'disc_z':   [1, [dim[2], 'Equal', delta[2]]],
}

"""
    Information on initial and boundary conditions
"""
#   Default initial conditions
#   This will be the INDOM block of the TOUGH INPUT
#   Each row of indom has the format : ['AqO', P, GOR, So, T] or ['AOG', P, Sg, So, T]
indom = {
            'Default'   :   ['AqO', 2.5e7, 3.0e3, 0.55, 65.0],
            'HWELL'     :   [1, 'AOG', 7.5e6, 0.4, 0.4, 65.0]
}

#   Initial conditions used in the INCON file
#   Each key contains a list of 3 elements, which represent [Matrix, HF, Well]
#   The SRV (if exists) has the same initial condition as the matrix
#   Note that when use_geos=True, the HF porosity and permeability will be
#   overwritten with values from GEOS.
incon = {
            'state'         :   ['AqO', 'AqO', 'AOG'],
            'porosity'      :   [0.07, 0.9, 1.0],
            'permeability'  :   [1e-21, 1e-21, 1e-4],
            'pressure'      :   [2.5e7, 2.5e7, 5.5e6],
            'gor'           :   [3.0e3, 3.0e3, 0.4],
            'saturation'    :   [0.55, 0.55, 0.4],
            'temperature'   :   [65.0, 65.0, 65.0]
}

"""
    Information on the GEOS-TOUGH coupling scheme
"""
#   Just a list of keys in the incon dict
ikeys = ['porosity', 'permeability', 'pressure', 'gor', 'saturation', 'temperature']
#   This determines (in sequence) if a field in incon will be imported from GEOS outputs
#   Normally only porosity and permeability are imported, but GEOS outputs also
#   contain pressure.
import_from_geos = [True, True, False, False, False, False]

#   This is the difference bewteen the y-coordinates of the HFs in GEOS and TOUGH
#   In GEOS, the HFs are often symmetric about y=0 (e.g. -32m, -16m, 0m, 16m, 32m)
#   In TOUGH, the y coordinates are all positive (e.g. 8m, 24m, 40m, 56m, 72m)
#   ygeos characterizes the difference between the two.
ygeos = [None, 40.0, 24.0, 8.0, -8.0, -24.0]

#   Information on the GEOS outputs
#   fgeos = name of the GEOS aperture file to be loaded
#   offset = differences between GEOS and TOUGH coordinates in [x, y, z]
if avg_aper != None and std_aper != None:
    fprop = 'aper'
elif use_clos:
    fprop = 'clos_aper'
else:
    fprop = 'prop_aper'


if sim_id == '09L':
    fgeos = 'GEOS09_upslo_'+fprop+'.csv'
    offset = [354.0, ygeos[ifrac], -120.0]
elif sim_id == '09H':
    fgeos = 'GEOS09_upshi_'+fprop+'.csv'
    offset = [354.0, ygeos[ifrac], -280.0]
elif sim_id == 'SMH':
    fgeos = 'GEOS03_samepump25_'+fprop+'.csv'
    offset = [354.0, ygeos[ifrac], -280.0]

#   Number of downscaled HFs in a swarm. This value should be consistent with
#   the value GEOS uses.
n_swarm = 8.0
#   The range within which the HF is located in the TOUGH domain
xlim = [20.0, 680.0, 7.9, 8.1, -(dim[2]-1)*delta[2], 0.0]
#   Minimum aperture allowed in TOUGH to avoid non-convergence
cutoff = 1e-4
#   This is the empirical parameter 'gamma' used for estimating aperture in the
#   isolated regions that are created by the stress-shadowing effect.
iso_correction = 1.2
#   blk_range is the instruction for removing the isolated regions.
#   blk_range = {
#        blk_key : [[x_start, x_end, z_start, z_end, direction],
#                   [next isolated region...]]
#       }
#   For each isolated region, the code searches for the minimum aperture within
#   range [x_start, x_end, z_start, z_end] and identifies the region in the
#   'direction' of the minimum aperture as the isolated region.
blk_range = {
    'up09l1'   :    [[180.0, 270.0, -180.0, -60.0, 'x-']],
    'up09l2'   :    [[240.0, 330.0, -180.0, -60.0, 'x-'],
                        [420.0, 510.0, -180.0, -60.0, 'x+']],
    'up09l3'   :    [[510.0, 570.0, -180.0, -60.0, 'x+']],
    'up09h1'   :    [[180.0, 660.0, -200.0, -140.0, 'z+'],
                        [510.0, 570.0, -240.0, -180.0, 'x+']],
    'up09h2'   :    [[180.0, 660.0, -200.0, -120.0, 'z+']],
    'up09h3'   :    [[180.0, 660.0, -200.0, -80.0, 'z+']],
    'up03h3'   :    [[180.0, 660.0, -200.0, -80.0, 'z+'],
                    [330.0, 390.0, -264.0, -248.0, 'z-'],
                    [180.0, 300.0, -320.0, -200.0, 'x-']]
}


""" ###########################################################################
                        Execution - No user change below
    ####################################################################### """

mesh = Mesh(dim, delta, region, boundary, discret, [])
mesh.write_input(use_minc)
mesh.build_mesh()

if use_minc:
    mesh = Mesh(dim, delta, region, boundary, discret, minc)
    mesh.write_input(use_minc)
    mesh.build_mesh()
    elem, conn = mesh.load_mesh('MINC')
else:
    elem, conn = mesh.load_mesh('MESH')

map = mesh.neighbor_map(elem, conn)
conx = mesh.get_boundary_id(elem, conn, map)

if use_geos:
    # Initialize dict for incon variables
    incon_var = {}
    for key in ikeys:
        incon_var[key] = []
    # Load GEOS results
    geos = Geos(fgeos, offset, xlim, cutoff)
    geos_aper = geos.custom_range()
    for ii in range(np.shape(geos_aper)[0]):
        if geos_aper[ii,3] < cutoff:
            geos_aper[ii,3] = cutoff

    if sim_id == 'SMH' and ifrac == 3 and blk_key == 'up03h3':
        print(' ----------> Special treatment for SMH-3!!!')
        for ii in range(np.shape(geos_aper)[0]):
            if geos_aper[ii,2] == -238.0:
                if geos_aper[ii,0] >= 39*6 and geos_aper[ii,0] <= 67*6:
                    geos_aper[ii,3] = 6e-4
            if abs(geos_aper[ii,2]-(-184)) <= 2.0 and abs(geos_aper[ii,0]-288) <= 3.0:
                geos_aper[ii,3] = 1.5e-4

    if avg_aper != None and std_aper != None:
        # scale to match target mean / std
        geos_aper = geos.scale(geos_aper, [avg_aper, std_aper], cutoff)

    if blk_key != None:
        geos_aper = geos.block_checking(geos_aper, blk_key, blk_range, iso_correction)

    if avg_aper == None and std_aper != None:
        # scale to match target mean / std
        geos_aper = geos.scale(geos_aper, [avg_aper, std_aper], cutoff)

    if sim_id == 'SMH' and ifrac == 3 and blk_key == 'up03h3':
        for ii in range(np.shape(geos_aper)[0]):
            if geos_aper[ii,2] > -166.0:
                geos_aper[ii,3] = np.nan
            if geos_aper[ii,2] < -262.0:
                geos_aper[ii,3] = np.nan
            if geos_aper[ii,0] < 231.0:
                geos_aper[ii,3] = np.nan

    # remove extreme aperture values
    for ii in range(np.shape(geos_aper)[0]):
        if geos_aper[ii,3] < cutoff:
            geos_aper[ii,3] = cutoff
        if geos_aper[ii,3] > delta[1]:
            print(geos_aper[ii,:])
            geos_aper[ii,3] = delta[1]*0.99

    # Calculate perm and poro
    geos_perm = []
    geos_poro = []
    phi = incon['porosity'][1]
    for ii in range(np.shape(geos_aper)[0]):
        geos_perm.append([geos_aper[ii,0], geos_aper[ii,1], geos_aper[ii,2],
            geos_aper[ii,3]**2.0 / (12.0 * n_swarm**2.0)])
        geos_poro.append([geos_aper[ii,0], geos_aper[ii,1], geos_aper[ii,2], geos_aper[ii,3]/delta[1]])

    incon_var['porosity'] = np.array(geos_poro)
    incon_var['permeability'] = np.array(geos_perm)
    print('GEOS data range: x --> [',np.amin(geos_aper[:,0]),',',np.amax(geos_aper[:,0]),']')
    print('GEOS data range: z --> [',np.amin(geos_aper[:,2]),',',np.amax(geos_aper[:,2]),']')
    print('GEOS in TOUGH domain: max, min, mean , std aperture = ', np.nanmax(geos_aper[:,3]), np.nanmin(geos_aper[:,3]),
        np.nanmean(geos_aper[:,3]), np.nanstd(geos_aper[:,3]))
    gmap, isfrac = mesh.geos_map_minc(elem, conn, geos_aper, frac_zone)
    if remesh:
        elem, conn = mesh.adjust_mesh_poro(elem, conn, gmap, map, geos_poro)
        mesh.update_mesh(elem, conn)

#
#   Insert GEOS aperture
#
if use_geos:
    data = np.zeros((mesh.Nele, len(ikeys)), dtype=float)
    for ii in range(len(ikeys)):
        variable = Incon(incon[ikeys[ii]][1], len(elem), import_from_geos[ii])
        out1d = variable.insert_geos(incon_var[ikeys[ii]], gmap)
        data[:,ii] = out1d

#
#   Calculate spatial-variable MINC permeability
#
if use_secf:
    data = variable.secf_model(data, elem, geos_aper, incon, secf, srv_frac_zone)

#
#   Write INPUT file
#
input = Input(elem, conn, use_minc)
if use_geos:
    if use_secf:
        incon_id = input.write_incon(data, 'AqO', frac_zone+srv_frac_zone, isfrac, secf)
    else:
        incon_id = input.write_incon(data, 'AqO', frac_zone, isfrac, secf)
input.write_memory('INPUT_OIL_TEMPLATE')
input.write_rocks(region, boundary)
input.write_param(timestep, indom)
input.write_indom(indom, use_minc)
input.write_times(t_save)
input.write_conx(conx)


#
#   Unify medias
#
if use_minc:
    mesh = Mesh(dim, delta, region, boundary, discret, [])
    elem, conn = mesh.load_mesh('MESH')
    for ii in range(len(elem)):
        if elem[ii][1] in frac_zone:
            if not elem[ii][0] in incon_id:
                elem[ii][1] = '    1'
    mesh.update_mesh(elem, conn)
