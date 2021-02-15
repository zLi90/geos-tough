""" Test MeshMaker """
import os
import subprocess as sp
import numpy as np
from BuildMesh import *
from Input import *
from Geos import *
from InconVariable import *
from utils import *

''' -------------------- BEGIN USER SETTINGS -------------------- '''
simname = '4Mbase'
ikeys = ['Porosity', 'Perm', 'Pressure', 'Saturation', 'SatOil', 'Temperature']

timestep = [3.024e8, 1.0, 8.64e5]
tVec = [8.64e4, 8.64e5, 8.64e6, 8.64e7, 2.592e8]

use_geos = True
nswarm = 8
r_adj = 16.0

dimension = {
            'delta'     :   [6.0, 1.0, 4.0],
            'dim'       :   [118, 135, 46],
            'offset'    :   [[354.0, 48.0, -10.0]],
            # 'dim'       :   [118, 135, 75],
            # 'offset'    :   [[354.0, 48.0, -132.0]],
            'ratio'     :   [1.0, 0.25, 1.0],
            'yfrac'     :   [-32.0,-16.0,0.0,16.0,32.0],
            'bound'     :   []
}
#   Settings for generating mesh using MeshMaker
fmesh = {
                'MaxElem'       :   dimension['dim'][0]*dimension['dim'][1]*dimension['dim'][2],
                'Longest'       :   200,
                'Regions'       :   ['MATRX','FRAC1','WELL1'],
                'NX'            :   [dimension['dim'][0], dimension['delta'][0]],
                'NY'            :   [dimension['dim'][1], [[10, 0.8]]],
                'NZ'            :   [dimension['dim'][2], dimension['delta'][2]]
}
for ii in range(5):
    fmesh['NY'][1].append([11, 0.025, 7.99, '-'])
    fmesh['NY'][1].append([1, 0.025])
    fmesh['NY'][1].append([11, 0.025, 7.985, '+'])
    # fmesh['NY'][1].append([1, 0.0055])
    # fmesh['NY'][1].append([11, 0.025, 7.985+0.0195, '+'])

fmesh['NY'][1].append([10, 0.8])

rocks = {
            'MATRX'     :   [['range', []],
                            ['properties', [2.6e3, 0.08, 1e-18, 1e-8, 8, 8]]],
            'FRAC1'     :   [['range'],
                            ['properties', [2.6e3, 0.8, 3.95e-8, 1e-8, 8, 8]]],
            'WELL1'     :   [['range', [354.0, 360.0, 0.0, 96.0, -86.0, -82.1]],
                            ['properties', [2.6e3, 1.00, 1e-6, 0e-8, 6, 5]]]
            # 'WELL1'     :   [['range', [354.0, 360.0, 0.0, 96.0, -200.0, -196.1]],
            #                 ['properties', [2.6e3, 1.00, 1e-6, 0e-8, 6, 5]]]
}
yfrac = 15.99
for ii in range(10):
    # low
    rocks['FRAC1'][0].append([30.0, 696.0, yfrac, yfrac+0.025, -168.0, -16.0])
    # high
    # rocks['FRAC1'][0].append([30.0, 696.0, yfrac, yfrac+0.025, -284.0, -16.0])
    # rectangular low
    # rocks['FRAC1'][0].append([165.0, 597.0, yfrac, yfrac+0.0055, -162.0, -26.0])
    yfrac += 16.0

ic = {
            'Default'   :   ['AqO', 3e7, 7.0e2, 0.5, 70.0],
            'FRAC1'     :   ['AqO', 3e7, 7.0e2, 0.5, 70.0],
            'WELL1'     :   ['AOG', 4e6, 0.4, 0.4, 70.0]
}

isource = {
                'state'         :   'AqO',
                'Porosity'      :   ['binary', 'none', [0.08,0.8], 'none'],
                # 'Perm'          :   ['binary', 'none', [1e-18,1e-10], 'none'],
                'Perm'          :   ['geos', 'GEOS_A5frac_upslo', [1e-18,1e-8], 'aper2perm'],
                'Pressure'      :   ['binary', 'none', [3e7,3e7], 'none'],
                'Saturation'    :   ['const', 'none', [7.0e2,7.0e2], 'none'],
                'SatOil'        :   ['binary', 'none', [0.5,0.5], 'none'],
                'Temperature'   :   ['const', 'none', [70.0,70.0], 'none']
}

bc = {
            'useBC'     :   0,
            'name'      :   ['FRAC1','WELL1'],
            'state'     :   ['AOG'],
            'isheet'    :   [0, 0, None, None, 1],
            'period'    :   ['20160115', '20160120'],
            'col'       :   [[0], [2], [None], [None], [2]],
            'units'     :   ['d2sec', 'psi2Pa', None, None, 'F2C']
}
wells = ['WELL1']
#   Names of the aperture data files
fgeos = {
        'Perm'      :   ['geos_5cluster_upslo.csv',[['GEOS_A5frac_upslo',None,None]],[]]
        # 'Pressure'  :   ['../../../data/GEOS_P.csv',[['GEOS_P'+simname+'.csv',0.625]],[88,100,44]]
}
#   Name of production data file
fdata = '../../../data/Production_by_well.xlsx'
''' -------------------- END OF USER SETTINGS -------------------- '''

if use_geos:
    geos = Geos(fgeos['Perm'], dimension)
    geos.getRange()
    data2D_interp = geos.interp2D()
    data1D = geos.from2Dto1D(data2D_interp)
    data1D = geos.duplicate(data1D)
    data_out = geos.scale(data1D, fgeos['Perm'][1][0])
    np.savetxt(fgeos['Perm'][1][0][0], data_out, delimiter=',')

#   Create mesh
mesh = Mesh(fmesh, rocks, dimension['delta'], dimension['dim'], dimension['bound'], bc, wells)
fmeshin = mesh.genMeshInput()
mesh.genMesh(fmeshin)
elem, conn = mesh.readMesh()
elem, flag_rewrite = mesh.coordNearFrac(elem,[16.0,32.0,48.0,64.0,80.0])
isbound, iswell = mesh.markBound(elem)
if flag_rewrite == 1:
    mesh.writeNewMesh(elem, conn, isbound)
    elem, conn = mesh.readMesh()
    isbound, iswell = mesh.markBound(elem)
#   Adjust mesh to fit fractures
map = mesh.getNeighbor(elem, conn)
if use_geos:
    gmap, isfrac, = mesh.buildMap(elem, conn, np.genfromtxt(fgeos['Perm'][1][0][0], delimiter=','), ['FRAC1'])
    elem, conn = mesh.adjustMesh(elem, conn, gmap, map, r_adj)
    elem = mesh.fracMask(elem, isfrac, map, ['FRAC1'], 'MATRX')
#   Add I labels to the mesh
conx = mesh.getWellInterface(elem, conn, map, iswell)
mask = mesh.wellMask(elem, ['MATRX'])
mesh.writeNewMesh(elem, conn, isbound)

#   Create INCON from GEOS aperture
if use_geos:
    data = np.zeros((mesh.Nele, len(ikeys)))
    for ii in range(0,len(ikeys)):
        variable = InconVariable(isource[ikeys[ii]], mesh.Nele, isfrac, True, nswarm)
        data1D = variable.processField(gmap)
        data[:,ii] = data1D[:,0]
        print("Variable "+ikeys[ii]+" has been processed!")
    writeInconOil(elem, isource['state'], data, mask, isource['Perm'][2][0], map, wells, ic)

#   Create boundary condition from Production data
databc = []
for ii in range(len(bc['name'])):
    one_bc = []
    one_bc.append(bc['name'][ii])
    for jj in range(len(ic[bc['name'][ii]])):
        one_bc.append(ic[bc['name'][ii]][jj])
    databc.append(one_bc)

#   Write new INPUT file
input = Input(elem, conn)
input.writeMEMORY('INPUT_OIL_TEMPLATE')
input.writeROCKS(rocks)
input.writePARAM(timestep, ic)
input.writeINDOM(databc)
input.writeTIMES(tVec)
input.writeCONX(conx['WELL1'])
