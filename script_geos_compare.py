""" Plot the GEOS aperture """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
from matplotlib.patches import Circle
import matplotlib
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mm_input import *

testid = '09'

if testid == 'LLP':
    fname = ['GEOS08_upslo_prop_aper.csv','GEOS08_upslo_prop_aper.csv','GEOS08_upslo_aper.csv']
    blk_key = ['upslo1','upslo2','upslo3']
    blk_check = [False, True, True]
    prop_scale = [[[None, None], [None, None], [3.9e-4, 3.4e-4]],
                    [[None, None], [None, None], [4.0e-4, 3.7e-4]],
                    [[None, None], [None, None], [3.8e-4, 3.4e-4]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
elif testid == 'HLP':
    fname = ['GEOS08_upshi_prop_aper.csv','GEOS08_upshi_prop_aper.csv','GEOS08_upshi_aper.csv']
    blk_key = ['upshi1','upshi2','upshi3']
    blk_check = [False, True, True]
    prop_scale = [[[None, None], [None, None], [2.8e-4, 2.5e-4]],
                    [[None, None], [None, None], [2.8e-4, 2.8e-4]],
                    [[None, None], [None, None], [2.75e-4, 2.8e-4]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
elif testid == '4U':
    fname = ['GEOS01_ups4u_prop_aper.csv','GEOS01_ups4u_prop_aper.csv','GEOS01_ups4u_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    prop_scale = [[[None, None], [None, None], [3.08e-4, 1.86e-4]],
                    [[None, None], [None, None], [3.07e-4, 2.13e-4]],
                    [[None, None], [None, None], [2.96e-4, 2.17e-4]]]
    offset = [[354.0, 40.0, -240.0],[354.0, 24.0, -240.0],[354.0, 8.0, -240.0]]
elif testid == '09':
    fname = ['GEOS09_upslo_prop_aper.csv','GEOS09_upslo_prop_aper.csv','GEOS09_upslo_prop_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    prop_scale = [[[None, None], [None, None], [3.08e-4, 1.86e-4]],
                    [[None, None], [None, None], [3.07e-4, 2.13e-4]],
                    [[None, None], [None, None], [2.96e-4, 2.17e-4]]]
    offset = [[354.0, 40.0, -240.0],[354.0, 24.0, -240.0],[354.0, 8.0, -240.0]]

savedata = False
savefig = False
max_val = 1.5

yfrac = ['16 m','32 m','48 m','64 m','80 m']
case = ['propped','propped+ROI','unpropped+ROI']
figind = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

cutoff = 1e-4
xlim = [0.0, 702.0, 7.9, 8.1, -280.0, 0.0]

iso_coef = 2.0
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

fs = 7

def get3Ddata(data):
    N = np.shape(data)[0]
    yy = np.unique(data[:,1])
    xx = np.linspace(3.0, 699.0, 117)
    zz = np.linspace(-2.0, -278.0, 70)
    dim = [len(xx),len(yy),len(zz)]
    geos3d = np.nan * np.ones((dim))
    print('Dimension of GEOS data = ',dim)

    for ll in range(N):
        ii = np.where(xx == data[ll,0])
        jj = np.where(yy == data[ll,1])
        kk = np.where(zz == data[ll,2])
        geos3d[ii,jj,kk] = data[ll,3]
    return geos3d, xx, yy, zz


aper3D = {}
xx = {}
yy = {}
zz = {}
ind = 0
for ii in range(len(fname)):
    for jj in range(3):
        key = str(ind)
        ind += 1
        geos = Geos(fname[ii], offset[jj], xlim, cutoff)
        geos_aper = geos.custom_range()
        for rr in range(np.shape(geos_aper)[0]):
            if geos_aper[rr,3] < cutoff:
                geos_aper[rr,3] = cutoff
        if blk_check[ii]:
            aper = geos.block_checking(geos_aper, blk_key[jj], blk_range[blk_key[jj]], iso_coef)
            aper_sc = geos.scale(aper, prop_scale[jj][ii], cutoff)
        else:
            aper_sc = geos.scale(geos_aper, prop_scale[jj][ii], cutoff)
        aper3D[key], xx[key], yy[key], zz[key] = get3Ddata(aper_sc)
        print('  >>>>> Mean, Std = ',np.nanmean(aper3D[key]),np.nanstd(aper3D[key]))


# Make plot
font = {'family' : 'Arial',
        'size'   : fs}
matplotlib.rc('font', **font)
cm = 1.0 / 2.54
plt.figure(1,figsize=(18*cm,12*cm))
ifig = 0
for col in range(3):
    for ii in range(len(fname)):
        # key = str(ii)
        key = str(ifig)
        ax = plt.subplot(3,len(fname),ifig+1)
        pos1 = ax.get_position()
        pos2 = [pos1.x0 + 0.035*(ii-1), pos1.y0 + 0.02*(1-col), pos1.width*1.1, pos1.height]
        # ax.set_position(pos2)
        # slice = np.transpose(geos_union[key][:,col,:])
        slice = np.transpose(aper3D[key][:,0,:])
        # aa = slice < (np.nanmax(slice)/5.0)
        aa = slice <= 1e-4
        slice[aa] = 1e-4


        # img=ax.imshow(np.transpose(slice)*1e3, cmap='jet', vmin=0.0, vmax=max_val, interpolation='bilinear')
        img=ax.imshow(slice*1e3, cmap='jet', vmin=0.0, vmax=max_val, interpolation='bilinear')
        # img=ax.imshow(slice*1e3, cmap='jet', vmin=0.0, vmax=max_val)


        if col == 2:
            plt.xlabel('X [m]',fontsize=fs)
            plt.xticks([0, 50, 100],['0','300','600'])
        else:
            plt.xticks([],[])
        if ii == 0:
            plt.ylabel('Z [m]',fontsize=fs)
        plt.yticks([0, 23, 47, 70],['0','-93','-186','-280'])
        if testid == 'LLP':
            well = plt.Circle((59,48),1.0,color='k')
        elif testid == 'HLP':
            well = plt.Circle((59,48),1.0,color='k')
        elif testid == '4U':
            well = plt.Circle((59,52),1.0,color='k')
        elif testid == '09':
            well = plt.Circle((59,52),1.0,color='k')
        plt.title(figind[ifig] + ' ' +testid+'-'+str(ii+1)+' ('+(case[col])+')',fontsize=fs)
        ax.add_artist(well)

        # plt.colorbar(cax=ax)
        if ii == 2:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.02)
            #
            plt.colorbar(img, cax=cax)

        ifig += 1
plt.savefig('temp.eps', format='eps')
if savefig:
    plt.savefig('Figure_5.eps', format='eps')




plt.show()

# Interpolate to desired range and resolution
#
