""" Create aperture for TOUGH from GEOS output """

import numpy as np
import silo_parser as sp
import matplotlib.pyplot as plt
from geos_convert import *

savedat = False
useprop = True
#
#   Select the GEOS results to be used
#
# ''' < upscaling with high landing point > '''
# id = 'upshi'
# fname = 'FiveFrac_base_652236'
# fdir = 'GEOS_upscaling202008/HighLandingPoint_b/sub'
# offset = [350.0, 40.0, -120.0]

# ''' < upscaling with low landing point > '''
# id = 'upslo'
# fname = 'FiveFrac_base_652236'
# fdir = 'GEOS_upscaling202008/LowLandingPoint_b/sub'
# offset = [350.0, 40.0, -120.0]

# ''' < upscaling with low landing point + 3D stress > '''
# id = 'upslo'
# fname = 'FiveFrac_base_571491'
# fdir = 'GEOS202009/4SM_LowLandingPoint/sub'
# offset = [350.0, 40.0, -120.0]

''' < upscaling with low landing point + 3D stress > '''
id = 'upshi'
fname = 'FiveFrac_base_562605'
fdir = 'GEOS202009/4SM_HighLandingPoint/sub'
offset = [350.0, 40.0, -240.0]

# ''' < upscaling with well 4SU > '''
# id = 'ups4u'
# fname = 'FiveFrac_base_788235'
# fdir = '4SU_40/sub'
# offset = [350.0, 40.0, -240.0]


savename = 'GEOS09_'+id+'_prop_aper.csv'
out_fields = ['Aperture','ProppantVolumeFraction']
# a list of all available field names:
# ['Aperture','FaceArea','Pressure','ProppantVolumeFraction','Volume','birthTime',
# 'permeability','stressNOnFace','totalLeakedThickness']


#
#   TOUGH domain settings
#
y_frac = [8.0, 24.0, 40.0, 56.0, 72.0]
dx = [6.0, 16.0, 4.0]
xlim = [0.0, 700.0, 0.0, 80.0, -276.0, 0.0]

#
#   Figure settings
#
y_plot = [0, 1, 2]
fs = 12
clim = [0.0, 0.02, 0.0, 0.6, 0.0, 0.002]

#
#   Execution
#

geos = geos_convert(fname, fdir, y_frac)

# extract required fields
data1d = {}
data1d_raw = {}
for field in out_fields:
    data1d_raw[field] = geos.get_one_field(field)
coord_geos = geos.get_coordinates(data1d_raw['Aperture'])

# interpolation for non-uniform grid resolution
# (currently non-uniform resolution is only found in z direction)
newz = geos.update_axis(coord_geos['zz'], dx[2])
for field in out_fields:
    data1d[field] = geos.interp_z(data1d_raw[field], coord_geos['zz'], newz)
coord_geos = geos.get_coordinates(data1d['Aperture'])

map = geos.get_actv_map(data1d['Aperture'])
aper = geos.apply_actv_map(data1d['Aperture'], map)
prop = geos.apply_actv_map(data1d['ProppantVolumeFraction'], map)

# calculate propped aperture
prop_aper = geos.prop_aper(aper, prop)
print('   >>> ------------')
print('Mean propped aperture = ',np.nanmean(prop_aper[:,3]))
print('Std propped aperture = ',np.nanstd(prop_aper[:,3]))
print('Max propped aperture = ',np.nanmax(prop_aper[:,3]))
print('Min propped aperture = ',np.nanmin(prop_aper[:,3]))
print('   >>> ------------')


# convert to 3D data for plotting
aper3D, coord = geos.convert_1D_to_3D(aper, xlim, dx, offset)
prop3D, coord = geos.convert_1D_to_3D(prop, xlim, dx, offset)
prop_aper3D, coord = geos.convert_1D_to_3D(prop_aper, xlim, dx, offset)

# print the coordinates
print('   >>> GEOS X: ',coord_geos['xx'])
print('   >>> GEOS Z: ',coord_geos['zz'])
print('   >>> TOUGH X: ',coord['xx'])
print('   >>> TOUGH Z: ',coord['zz'])

#
#   Save data
#
if savedat:
    if useprop:
        np.savetxt(savename, prop_aper, delimiter=',')
    else:
        np.savetxt(savename, aper, delimiter=',')

#   Grouping
levels = [5e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4]
slice = prop_aper3D[:,2,:]
dim = np.shape(slice)
slice_level = len(levels) * np.ones(dim)
for ii in range(dim[0]):
    for jj in range(dim[1]):
        kk = 0
        while True:
            if slice[ii,jj] == 0 or np.isnan(slice[ii,jj]):
                slice_level[ii,jj] = np.nan
                break
            else:
                if slice[ii,jj] <= levels[kk]:
                    slice_level[ii,jj] = kk
                    break
                else:
                    kk += 1
                if kk >= len(levels):
                    break
        if np.isnan(slice[ii,jj]):
            slice_level[ii,jj] = np.nan
#    Calculate area reduction
for ii in range(len(levels)):
    aa = slice_level == ii
    print('Cutoff = ', levels[ii], '  Total area = ',np.nansum(aa))
aa = np.isnan(slice_level)
print('Total fracture area = ', dim[0]*dim[1] - np.nansum(aa))

#
#   Make plot
#

plt.figure(1)

ifig = 1
for ii in range(len(y_plot)):
    plt.subplot(len(y_plot), 3, ifig)
    plt.imshow(np.flipud(np.transpose(aper3D[:,y_plot[ii],:])),
        vmin=clim[0], vmax=clim[1], interpolation='bilinear', cmap='jet')
    plt.colorbar()
    ifig += 1

    plt.subplot(len(y_plot), 3, ifig)
    plt.imshow(np.flipud(np.transpose(prop3D[:,y_plot[ii],:])),
        vmin=clim[2], vmax=clim[3], interpolation='bilinear', cmap='jet')
    plt.colorbar()
    ifig += 1

    plt.subplot(len(y_plot), 3, ifig)
    plt.imshow(np.flipud(np.transpose(prop_aper3D[:,y_plot[ii],:])),
        vmin=clim[4], vmax=clim[5], interpolation='bilinear', cmap='jet')
    plt.colorbar()
    ifig += 1


plt.figure(2)
plt.imshow(np.flipud(np.transpose(slice_level)), cmap='jet')
plt.colorbar()


plt.show()
