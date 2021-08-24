""" Plot the GEOS aperture """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
from matplotlib.patches import Circle
import matplotlib
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mm_input import *

testid = '09L'

if testid == 'LLP':
    fname = ['GEOS08_upslo_prop_aper.csv','GEOS08_upslo_prop_aper.csv','GEOS08_upslo_aper.csv']
    blk_key = ['upslo1','upslo2','upslo3']
    blk_check = [False, True, True]
    prop_scale = [[[None, None], [None, None], [3.9e-4, 3.4e-4]],
                    [[None, None], [None, None], [4.0e-4, 3.7e-4]],
                    [[None, None], [None, None], [3.8e-4, 3.4e-4]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -280.0, 0.0]
    dim = [117, 35, 70]
elif testid == 'HLP':
    fname = ['GEOS08_upshi_prop_aper.csv','GEOS08_upshi_prop_aper.csv','GEOS08_upshi_aper.csv']
    blk_key = ['upshi1','upshi2','upshi3']
    blk_check = [False, True, True]
    prop_scale = [[[None, None], [None, None], [2.8e-4, 2.5e-4]],
                    [[None, None], [None, None], [2.8e-4, 2.8e-4]],
                    [[None, None], [None, None], [2.75e-4, 2.8e-4]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -280.0, 0.0]
    dim = [117, 35, 70]
elif testid == '4U':
    fname = ['GEOS01_ups4u_prop_aper.csv','GEOS01_ups4u_prop_aper.csv','GEOS01_ups4u_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    prop_scale = [[[None, None], [None, None], [3.08e-4, 1.86e-4]],
                    [[None, None], [None, None], [3.07e-4, 2.13e-4]],
                    [[None, None], [None, None], [2.96e-4, 2.17e-4]]]
    offset = [[354.0, 40.0, -240.0],[354.0, 24.0, -240.0],[354.0, 8.0, -240.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -280.0, 0.0]
    dim = [117, 35, 70]
elif testid == '09L':
    testlab = 'LLP'
    fname = ['GEOS09_upslo_clos_aper.csv','GEOS09_upslo_prop_aper.csv','GEOS09_upslo_prop_aper.csv']
    blk_key = ['up09l1', 'up09l2', 'up09l3']
    blk_check = [False, False, True]
    blk_result = [[225.0, 0.146], [291.0, 0.278, 465.0, 0.234], [555.0, 0.132]]
    prop_scale = [[[3.94e-4, 2.91e-4], [None, None], [None, None]],
                    [[4.0e-4, 3.1e-4], [None, None], [None, None]],
                    [[3.80e-4, 2.99e-4], [None, None], [None, None]]]
    prop_scale = [[[None, None], [None, None], [None, None]],
                    [[None, None], [None, None], [None, None]],
                    [[None, None], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -180.0, 0.0]
    dim = [117, 35, 45]
elif testid == '09H':
    testlab = 'HLP'
    fname = ['GEOS09_upshi_aper.csv','GEOS09_upshi_prop_aper.csv','GEOS09_upshi_prop_aper.csv']
    blk_key = ['up09h1', 'up09h2', 'up09h3']
    blk_check = [False, False, True]
    blk_result = [[-154.0, 0.103], [-162.0, 0.113], [-166.0, 0.118]]
    prop_scale = [[[2.16e-4, 1.59e-4],[None, None], [None, None]],
                    [[2.10e-4, 1.32e-4],[None, None], [None, None]],
                    [[2.17e-4, 1.42e-4], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -280.0],[354.0, 24.0, -280.0],[354.0, 8.0, -280.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -340.0, 0.0]
    dim = [117, 35, 85]
elif testid == 'HFL1':
    testlab = 'HFL1'
    fname = ['GEOS09_hiflu1_aper.csv','GEOS09_hiflu1_prop_aper.csv','GEOS09_hiflu1_prop_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    blk_result = [[-154.0, 0.103], [-162.0, 0.113], [-166.0, 0.118]]
    prop_scale = [[[2.16e-4, 1.59e-4],[None, None], [None, None]],
                    [[2.10e-4, 1.32e-4],[None, None], [None, None]],
                    [[2.17e-4, 1.42e-4], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -280.0],[354.0, 24.0, -280.0],[354.0, 8.0, -280.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -340.0, 0.0]
    dim = [117, 35, 85]
elif testid == 'HFL2':
    testlab = 'HFL2'
    fname = ['GEOS09_hiflu2_aper.csv','GEOS09_hiflu2_prop_aper.csv','GEOS09_hiflu2_prop_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    blk_result = [[-154.0, 0.103], [-162.0, 0.113], [-166.0, 0.118]]
    prop_scale = [[[2.16e-4, 1.59e-4],[None, None], [None, None]],
                    [[2.10e-4, 1.32e-4],[None, None], [None, None]],
                    [[2.17e-4, 1.42e-4], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -280.0],[354.0, 24.0, -280.0],[354.0, 8.0, -280.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -340.0, 0.0]
    dim = [117, 35, 85]
elif testid == 'SML':
    testlab = 'SML'
    fname = ['GEOS03_samepump37_aper.csv','GEOS03_samepump37_prop_aper.csv','GEOS03_samepump37_prop_aper.csv']
    blk_key = [None, None, None]
    blk_check = [False, False, False]
    blk_result = [[-154.0, 0.103], [-162.0, 0.113], [-166.0, 0.118]]
    prop_scale = [[[2.16e-4, 1.59e-4],[None, None], [None, None]],
                    [[2.10e-4, 1.32e-4],[None, None], [None, None]],
                    [[2.17e-4, 1.42e-4], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -120.0],[354.0, 24.0, -120.0],[354.0, 8.0, -120.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -340.0, 0.0]
    dim = [117, 35, 85]
elif testid == 'SMH':
    testlab = 'SMH'
    fname = ['GEOS03_samepump25_aper.csv','GEOS03_samepump25_prop_aper.csv','GEOS03_samepump25_prop_aper.csv']
    blk_key = ['up03h1', 'up03h2', 'up03h3']
    blk_check = [False, False, True]
    blk_result = [False, False, [-166.0, 0.113, -262.0, 0.116, 231.0, 0.245]]
    prop_scale = [[[2.29e-4, 1.70e-4],[None, None], [None, None]],
                    [[2.28e-4, 1.54e-4],[None, None], [None, None]],
                    [[2.31e-4, 1.55e-4], [None, None], [None, None]]]
    offset = [[354.0, 40.0, -280.0],[354.0, 24.0, -280.0],[354.0, 8.0, -280.0]]
    xlim = [0.0, 702.0, 7.9, 8.1, -340.0, 0.0]
    dim = [117, 35, 85]

savedata = False
savefig = False
max_val = 1.2

yfrac = ['16 m','32 m','48 m','64 m','80 m']
case = ['unpropped','propped','propped+modified']
figind = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

cutoff = 1e-7


iso_coef = 1.0
blk_range = {
    'upslo1'   :    [[450.0, 570.0, -280.0, -160.0, 'x+']],
    'upslo2'   :    [[120.0, 240.0, -280.0, -160.0, 'x-']],
    'upslo3'   :    [[450.0, 570.0, -280.0, -160.0, 'x+']],
    'upshi1'   :    [[175.0, 270.0, -280.0, -180.0, 'x-'],
                        [330.0, 390.0, -200.0, -160.0, 'z+']],
    'upshi2'   :    [[450.0, 570.0, -280.0, -180.0, 'x+'],
                        [270.0, 420.0, -180.0, -160.0, 'z+']],
    'upshi3'   :    [[180.0, 270.0, -280.0, -180.0, 'x-'],
                        [300.0, 480.0, -180.0, -160.0, 'z+']],
    'up09l1'   :    [[180.0, 270.0, -180.0, -60.0, 'x-']],
    'up09l2'   :    [[240.0, 330.0, -180.0, -60.0, 'x-'],
                        [420.0, 510.0, -180.0, -60.0, 'x+']],
    'up09l3'   :    [[510.0, 570.0, -180.0, -60.0, 'x+']],
    'up09h1'   :    [[180.0, 660.0, -200.0, -140.0, 'z+']],
    'up09h2'   :    [[180.0, 660.0, -200.0, -120.0, 'z+']],
    'up09h3'   :    [[180.0, 660.0, -200.0, -80.0, 'z+']],
    'up03h1'   :    [[180.0, 660.0, -200.0, -140.0, 'z+']],
    'up03h2'   :    [[180.0, 660.0, -200.0, -120.0, 'z+']],
    'up03h3'   :    [[180.0, 660.0, -200.0, -80.0, 'z+'],
                    [330.0, 390.0, -264.0, -248.0, 'z-'],
                    [180.0, 300.0, -320.0, -200.0, 'x-']]
}

fs = 7

def get3Ddata(data, xvec, zvec):
    N = np.shape(data)[0]
    yvec = np.unique(data[:,1])
    dim3 = [len(xvec),len(yvec),len(zvec)]
    geos3d = np.nan * np.ones((dim3))
    print('Dimension of GEOS data = ',dim3)

    for ll in range(N):
        ii = np.where(xvec == data[ll,0])
        jj = np.where(yvec == data[ll,1])
        kk = np.where(zvec == data[ll,2])
        geos3d[ii,jj,kk] = data[ll,3]
    return geos3d, xvec, yvec, zvec


aper3D = {}
xx = {}
yy = {}
zz = {}
ind = 0
for jj in range(3):
    print('--------------------')
    for ii in range(len(fname)):
        key = str(ind)
        ind += 1
        geos = Geos(fname[ii], offset[jj], xlim, cutoff)
        geos_aper = geos.custom_range()
        for rr in range(np.shape(geos_aper)[0]):
            if geos_aper[rr,3] < cutoff:
                geos_aper[rr,3] = cutoff
        if blk_check[ii]:
            aper = geos.scale(geos_aper, prop_scale[jj][ii], cutoff)
            aper_sc = geos.block_checking(aper, blk_key[jj], blk_range, iso_coef)
            # for kk in range(np.shape(aper_sc)[0]):
                # if aper_sc[kk,2] < -220.0 and aper_sc[kk,2] > -255.0 and aper_sc[kk,3] < 4e-4:
                #     if aper_sc[kk,0] > 340.0 and aper_sc[kk,0] < 375.0:
                #         aper_sc[kk,3] = np.amax([4e-4, aper_sc[kk,3]])
                # if aper_sc[kk,3] > 6e-4:
                #     aper_sc[kk,3] = 6e-4
                # if aper_sc[kk,2] > -200.0 and aper_sc[kk,3] > 1.2e-4:
                #     aper_sc[kk,3] = 1.2e-4
                # if aper_sc[kk,2] > -150:
                #     aper_sc[kk,3] = np.nan
            slice = []
            for kk in range(np.shape(aper_sc)[0]):
                if aper_sc[kk,2] == -238.0:
                    if aper_sc[kk,0] >= 39*6 and aper_sc[kk,0] <= 67*6:
                        slice.append(aper_sc[kk,3])
                        aper_sc[kk,3] = 6e-4
                if abs(aper_sc[kk,2]-(-184)) <= 2.0 and abs(aper_sc[kk,0]-288) <= 3.0:
                    aper_sc[kk,3] = 1.5e-4
            print('SLICE :::')
            slice = np.array(slice)
            ns = len(slice)
            print(np.nanmean(slice))
            print(np.nansum(slice**3.0), np.nansum(ns*(np.nanmean(slice)**3.0)), np.nansum(ns*(6e-4**3.0)))




        else:
            aper_sc = geos.scale(geos_aper, prop_scale[jj][ii], cutoff)

        xvec = np.linspace(3.0, 699.0, dim[0])
        zvec = np.linspace(-2.0, xlim[4]+2.0, dim[2])
        aper3D[key], xx[key], yy[key], zz[key] = get3Ddata(aper_sc, xvec, zvec)
        print('  >>>>> Mean, Std = ',np.nanmean(aper3D[key]),np.nanstd(aper3D[key]))
print(xvec)
print(zvec)


# Make plot
font = {'family' : 'Arial',
        'size'   : fs}
matplotlib.rc('font', **font)
cm = 1.0 / 2.54
if testid == '09L':
    max_val = 1.2
    plt.figure(1,figsize=(18*cm,8*cm))
else:
    plt.figure(1,figsize=(18*cm,14*cm))
ifig = 0

for ii in range(len(fname)):
    for col in range(3):
        # key = str(ii)
        key = str(ifig)
        ax = plt.subplot(len(fname),3,ifig+1)
        pos1 = ax.get_position()
        pos2 = [pos1.x0 + 0.03*(col-1), pos1.y0 + 0.02*(1-ii), pos1.width*1.1, pos1.height]
        ax.set_position(pos2)
        # slice = np.transpose(geos_union[key][:,col,:])
        slice = np.transpose(aper3D[key][:,0,:])
        # aa = slice < (np.nanmax(slice)/5.0)
        # aa = slice <= 1e-4
        # slice[aa] = 1e-4


        # img=ax.imshow(np.transpose(slice)*1e3, cmap='jet', vmin=0.0, vmax=max_val, interpolation='bilinear')
        img=ax.imshow(slice*1e3, cmap='jet', vmin=0.0, vmax=max_val, interpolation='bilinear')
        # img=ax.imshow(slice*1e3, cmap='jet', vmin=0.0, vmax=max_val)

        # plot the modified zone
        if col == 2:
            if testid == '09L':
                plt.plot([blk_result[ii][0]/6.0, blk_result[ii][0]/6.0], [dim[2]-1, 2], 'r--')
                plt.text(blk_result[ii][0]/6.0-35, 10, str(blk_result[ii][1])+' mm', color='r', fontsize=fs)
                if len(blk_result[ii]) > 2:
                    plt.plot([blk_result[ii][2]/6.0, blk_result[ii][2]/6.0], [dim[2]-1, 2], 'r--')
                    plt.text(blk_result[ii][2]/6.0+5, 10, str(blk_result[ii][3])+' mm', color='r', fontsize=fs)
            elif testid == '09H':
                plt.plot([10, 100], [-blk_result[ii][0]/4.0, -blk_result[ii][0]/4.0], 'r--')
                plt.text(2, -blk_result[ii][0]/4.0-3, str(blk_result[ii][1])+' mm', color='r', fontsize=fs)
            elif testid == 'SMH':
                if ii == 2:
                    plt.plot([10, 100], [-blk_result[ii][0]/4.0, -blk_result[ii][0]/4.0], 'r--')
                    plt.text(2, -blk_result[ii][0]/4.0-3, str(blk_result[ii][1])+' mm', color='r', fontsize=fs)
                    plt.plot([10, 100], [-blk_result[ii][2]/4.0, -blk_result[ii][2]/4.0], 'r--')
                    plt.text(2, -blk_result[ii][2]/4.0+6, str(blk_result[ii][3])+' mm', color='r', fontsize=fs)
                    plt.plot([blk_result[ii][4]/6.0, blk_result[ii][4]/6.0], [dim[2]-1, 35], 'r--')
                    plt.text(blk_result[ii][4]/6.0-32, 82, str(blk_result[ii][5])+' mm', color='r', fontsize=fs)


        if ii == 2:
            plt.xlabel('X [m]',fontsize=fs)
            plt.xticks([0, 50, 100],['0','300','600'])
        else:
            plt.xticks([],[])
        if col == 0:
            plt.ylabel('Z [m]',fontsize=fs)

        if testid == '09L':
            plt.yticks([0, 23, 45],['0','-93','-180'])
        else:
            plt.yticks([0, 23, 47, 70],['0','-93','-186','-280'])

        if testid == 'LLP':
            well = plt.Circle((59,48),1.0,color='k')
        elif testid == 'HLP':
            well = plt.Circle((59,48),1.0,color='k')
        elif testid == '4U':
            well = plt.Circle((59,52),1.0,color='k')
        elif testid == '09L' or testid == 'SML':
            well = plt.Circle((59,23),1.0,color='k')
        elif testid == '09H' or testid == 'HFL1' or testid == 'HFL2' or testid == 'SMH':
            well = plt.Circle((59,62),1.0,color='k')
        plt.title(figind[ifig] + ' ' +testlab+'-'+str(ii+1)+' ('+(case[col])+')',fontsize=fs)
        ax.add_artist(well)

        # plt.colorbar(cax=ax)
        if col == 2:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.02)
            #
            plt.colorbar(img, cax=cax)

        ifig += 1

# plt.savefig('HLP.eps', format='eps')
if savefig:
    plt.savefig('Figure_5_SMH.eps', format='eps')




plt.show()

# Interpolate to desired range and resolution
#
