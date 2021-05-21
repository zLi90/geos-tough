""" Create aperture for TOUGH from GEOS output """

import numpy as np
import silo_parser as sp
import matplotlib.pyplot as plt
from geos_convert import *

savedat = True
useprop = True
getstress = True
hf_clos = True
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

''' < upscaling with low landing point + 3D stress > '''
id = 'upslo'
fname = 'FiveFrac_base_571491'
fdir = 'GEOS202009/4SM_LowLandingPoint/sub'
offset = [350.0, 40.0, -120.0]

# ''' < upscaling with high landing point + 3D stress > '''
# id = 'upshi'
# fname = 'FiveFrac_base_562605'
# fdir = 'GEOS202009/4SM_HighLandingPoint/sub'
# offset = [350.0, 40.0, -280.0]

# ''' < upscaling with well 4SU > '''
# id = 'ups4u'
# fname = 'FiveFrac_base_788235'
# fdir = '4SU_40/sub'
# offset = [350.0, 40.0, -240.0]

# ''' < high landing point + early flush > '''
# id = 'hiflu1'
# fname = 'FiveFrac_base_333123'
# fdir = 'GEOS_proppants202103'
# offset = [350.0, 40.0, -280.0]

# ''' < high landing point + late flush > '''
# id = 'hiflu2'
# fname = 'FiveFrac_base_344050'
# fdir = 'GEOS_proppants202103'
# offset = [350.0, 40.0, -280.0]

# ''' < stage25 with same pumping schedule (HLP) > '''
# id = 'samepump25'
# fname = 'FiveFrac_base_595388'
# fdir = 'GEOS_same_pumping202103/stage25_20210322/sub/'
# offset = [350.0, 40.0, -280.0]

# ''' < stage37 with same pumping schedule (LLP) > '''
# id = 'samepump37'
# fname = 'FiveFrac_base_604792'
# fdir = 'GEOS_same_pumping202103/stage37_20210321/sub/'
# offset = [350.0, 40.0, -120.0]

savename = 'GEOS09_'+id+'_aper.csv'
if useprop:
    savename = 'GEOS09_'+id+'_clos_aper.csv'
if getstress:
    strename = 'GEOS09_'+id+'_stre.csv'
out_fields = ['Aperture','ProppantVolumeFraction','stressNOnFace']
# a list of all available field names:
# ['Aperture','FaceArea','Pressure','ProppantVolumeFraction','Volume','birthTime',
# 'permeability','stressNOnFace','totalLeakedThickness']


#
#   TOUGH domain settings
#
y_frac = [8.0, 24.0, 40.0, 56.0, 72.0]
dx = [6.0, 16.0, 4.0]
xlim = [0.0, 700.0, 0.0, 80.0, -340.0, 0.0]
# xlim = [0.0, 700.0, 0.0, 80.0, -276.0, 0.0]

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
    # data1d[field] = data1d_raw[field]
    data1d[field] = geos.interp_z(data1d_raw[field], coord_geos['zz'], newz, coord_geos['xx'], coord_geos['yy'], 154.0)
coord_geos = geos.get_coordinates(data1d['Aperture'])

map = geos.get_actv_map(data1d['Aperture'])
aper = geos.apply_actv_map(data1d['Aperture'], map)
prop = geos.apply_actv_map(data1d['ProppantVolumeFraction'], map)
stre = geos.apply_actv_map(data1d['stressNOnFace'], map)

# calculate propped aperture
prop_aper = geos.prop_aper(aper, prop, hf_clos)
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

stre3D, coord = geos.convert_1D_to_3D(stre, xlim, dx, offset)

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
        # for kk in range(np.shape(prop_aper)[0]):
        #     if prop_aper[kk,2] < -220.0+280.0 and prop_aper[kk,3] < 2e-4:
        #         if prop_aper[kk,0] > 300.0-354.0 and prop_aper[kk,0] < 450.0-354.0:
        #             prop_aper[kk,3] = np.amax([2e-4, prop_aper[kk,3]])

        np.savetxt(savename, prop_aper, delimiter=',')
        if getstress:
            np.savetxt(strename, stre, delimiter=',')
    else:
        np.savetxt(savename, aper, delimiter=',')
        if getstress:
            np.savetxt(strename, stre, delimiter=',')

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

    # plt.subplot(len(y_plot), 3, ifig)
    # plt.imshow(np.flipud(np.transpose(prop_aper3D[:,y_plot[ii],:])),
    #     vmin=clim[4], vmax=clim[5], interpolation='bilinear', cmap='jet')
    # plt.colorbar()
    # ifig += 1

    plt.subplot(len(y_plot), 3, ifig)
    plt.imshow(np.flipud(np.transpose(-stre3D[:,y_plot[ii],:])),
        vmin=0.0, vmax=4e7, interpolation='bilinear', cmap='jet')
    plt.colorbar()
    ifig += 1


plt.figure(2)
plt.imshow(np.flipud(np.transpose(slice_level)), cmap='jet')
plt.colorbar()


plt.show()
