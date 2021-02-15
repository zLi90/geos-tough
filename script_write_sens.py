""" Test the mm_input scripts """
import numpy as np
from mm_input import *
import subprocess as sp

use_geos = True
use_minc = True
remesh = True

#
#   Dimension of model domain
#
min_ratio = 2.0
dim = [117, 35, 70]
delta = [6.0, 0.002, 4.0]
yfrac = 8.0
#
#   Time control of simulation
#
timestep = [1.296e8, 1.0, 8.64e5]
tVec = [8.64e4, 8.64e5, 8.64e6, 8.64e7, 1.296e8]

#
#   Mediums and regionsr
#
region = [4, ['MMATX', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        ['properties', 2.6e3, 0.07, 1e-20, 1e-8, 8, 8]]
             ],
            ['FSRV1', [20.0, 680.0, 6.0, yfrac-0.5*delta[1], -272.0, -8.0,
                        ['properties', 2.6e3, 0.07, 1e-20, 1e-8, 8, 8],
                        ['dualporosity', 2.6e3, 0.9, 1e-18, 1e-8, 8, 8]]
            ],
            ['FSRV2',   [20.0, 680.0, yfrac+0.5*delta[1], 10.0, -272.0, -8.0,
                        ['properties', 2.6e3, 0.07, 1e-20, 1e-8, 8, 8],
                        ['dualporosity', 2.6e3, 0.9, 1e-18, 1e-8, 8, 8]]
            ],
            ['MHFRA', [20.0, 680.0, yfrac-0.5*delta[1], yfrac+0.5*delta[1], -272.0, -8.0,
                        ['properties', 2.6e3, 0.07, 1e-18, 1e-8, 8, 8]]
            ]
          ]
# high landing point
# boundary = [1, ['HWELL', [354.0, 360.0, 0.0, 16.0, -188.0, -184.0,
#                 ['properties', 2.6e3, 1.0, 1e-4, 0e-8, 6, 5]]]
#             ]
# low landing point
# boundary = [1, ['HWELL', [354.0, 360.0, 0.0, 16.0, -196.0, -192.0,
#                 ['properties', 2.6e3, 1.0, 1e-4, 0e-8, 6, 5]]]
#             ]
# GEOS09 low landing point
boundary = [1, ['HWELL', [354.0, 360.0, 0.0, 16.0, -212.0, -208.0,
                ['properties', 2.6e3, 1.0, 1e-4, 0e-8, 6, 5]]]
            ]


#
#   MINC settings
#
#[number_of_media, numger_of_fractured_media, number_of_MINCs,
#    number_of_specified_VolFractions, matrix_to_matrix_flow,
#      [fractured_medium_number, [], volume_fraction] ]
minc = [4, 2, 2, 1, 'A',
            [2, [5.0, 5.0, 5.0, 5.0], 0.01],
            [3, [5.0, 5.0, 5.0, 5.0], 0.01]]

#
#   Discretization
#
discret = {
        'disc_x':   [1, [dim[0], 'Equal', delta[0]]],
        'disc_y':   [9, [1, 'Uniform', 1.0], [14, 'Log', 1.0, 7.9, 1.0],
                        [1, 'Uniform', 0.078], [1, 'Uniform', 0.021+0.5*(0.002-delta[1])],
                        [1, 'Uniform', delta[1]],
                        [1, 'Uniform', 0.021+0.5*(0.002-delta[1])], [1, 'Uniform', 0.078],
                        [14, 'Log', 1.0, 15.0, 0.1865], [1, 'Uniform', 1.0]],
        'disc_z':   [1, [dim[2], 'Equal', delta[2]]],
}


# Initial conditions
indom = {
            'Default'   :   ['AqO', 2.5e7, 6.0e2, 0.53, 65.0],
            'HWELL'     :   [1, 'AOG', 5.5e6, 0.4, 0.4, 65.0]
}
# index of fracture zone in incon
idx_frac = 1
incon = {
            'state'         :   ['AqO', 'AqO', 'AOG'],
            'porosity'      :   [0.07, 0.07, 1.0],
            'permeability'  :   [1e-18, 1e-18, 1e-4],
            'pressure'      :   [2.5e7, 2.5e7, 5.5e6],
            'gor'           :   [6.0e2, 6.0e2, 0.4],
            'saturation'    :   [0.53, 0.53, 0.4],
            'temperature'   :   [65.0, 65.0, 65.0]
}
# Keys in the incon dict
ikeys = ['porosity', 'permeability', 'pressure', 'gor', 'saturation', 'temperature']
# Whether or not import incon from geos
import_from_geos = [True, True, False, False, False, False]

# Medium ID for hydraulic fracture
if use_minc:
    frac_zone = ['    7']
else:
    frac_zone = ['MHFRA']

# GEOS related settings
# fgeos = 'GEOS08_upshi_prop_aper.csv'
# offset = [354.0, 40.0, -120.0]

fgeos = 'GEOS09_upslo_aper.csv'
offset = [354.0, 40.0, -240.0]

n_swarm = 8.0
xlim = [20.0, 680.0, 7.9, 8.1, -270.0, 0.0]
avg_aper = 4.16e-4
std_aper = 2.91e-4
cutoff = 1e-4
cutoff2 = 1e-4
iso_correction = 2.0
blk_key = None
blk_range = {
    'upslo1'   :    [[450.0, 570.0, -280.0, -160.0, 'x+']],
    'upslo2'   :    [[120.0, 240.0, -280.0, -160.0, 'x-']],
    'upslo3'   :    [[450.0, 570.0, -280.0, -160.0, 'x+']],
    'upshi1'   :    [[175.0, 270.0, -280.0, -180.0, 'x-'],
                        [330.0, 390.0, -200.0, -160.0, 'z+']],
    'upshi2'   :    [[450.0, 570.0, -280.0, -180.0, 'x+'],
                        [270.0, 420.0, -180.0, -160.0, 'z+']],
    'upshi3'   :    [[180.0, 270.0, -280.0, -180.0, 'x-'],
                        [300.0, 480.0, -180.0, -160.0, 'z+']]
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
    if blk_key != None:
        geos_aper = geos.block_checking(geos_aper, blk_key, blk_range[blk_key], iso_correction)

    # scale to match target mean / std
    geos_aper = geos.scale(geos_aper, [avg_aper, std_aper], cutoff)


    # remove small aperture
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
        variable = Incon(incon[ikeys[ii]][idx_frac], len(elem), import_from_geos[ii])
        out1d = variable.insert_geos(incon_var[ikeys[ii]], gmap)
        data[:,ii] = out1d

#
#   Write INPUT file
#
input = Input(elem, conn, use_minc)
if use_geos:
    input.write_incon(data, 'AqO', frac_zone, isfrac, incon['permeability'][0])
input.write_memory('INPUT_OIL_TEMPLATE')
input.write_rocks(region, boundary)
input.write_param(timestep, indom)
input.write_indom(indom, use_minc)
input.write_times(tVec)
input.write_conx(conx)

#
#   Put inputs into a folder
#
# sp.call(['mv', 'INPUT', 'INPUT0'])
