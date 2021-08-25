""" Create contour plots from Plot_Data_Elem """
import numpy as np
import matplotlib.pyplot as plt
import os.path
from tough_mapping import *
# from matplotlib.colors import LogNorm
from matplotlib import colors

#
# User settings
#
kernels = ['']

folders = ['/Users/zli/Documents/HFTS/results_fracture_network/upslo2_sf2d_100-100-k15_1500_monitor/',
        '/Users/zli/Documents/HFTS/results_fracture_network/upslo2_sf2d_100-100-k15_hf1500/']
# folders = ['sens_ntwk/upslo2mink_ksrv1e-20_1e-7_ani1000_L100_cf12/','sens_ntwk/upslo2mink_ksrv1e-20_5e-8_ani100_L100_cf12/']
flabels = ['No closure','HF closed']

savefig = False

dim = [117,35,45]
dim_order = [1,2,0]
T = [1.296e8]
Z = [17]

fields = ['P','S_aqu','S_org','S_gas']
# fields = ['P','S_aqu','T','S_gas']
figlabel = ['Pressure [MPa]','Water saturation','Oil saturation','Gas saturation']
fs = 12
ls = ['-','--',':','--',':',':']
co = ['b','r']
mk = ['o','+','x']
lb = ['at 1d','at 1000d']

#
# Load data
#
map = buildMap(dim,dim_order)
data3D = get3Ddata(folders, dim, dim_order, fields, T)


ifig = 1
plt.figure(ifig,figsize=(8,5))

subfig = 1
for ifolder in range(len(folders)):
    plt.subplot(len(folders), 1, subfig)
    slice = np.transpose(data3D[folders[ifolder]]['P'][:,22,:,0]) / 1e6
    plt.imshow(slice, vmin=22, vmax=25, cmap='jet')
    plt.colorbar()
    plt.ylabel('Y [m]',fontsize=fs)
    plt.xlabel('X [m]',fontsize=fs)
    xlabel = [-354,-177,0,177,354]
    xtick = np.linspace(0,dim[0]-1,5)
    # ylabel = [0,-46,-92,-138]
    ylabel = [0,8,16]
    ytick = np.linspace(0,dim[1]-1,3)
    plt.xticks(xtick,xlabel)
    plt.yticks(ytick,ylabel)
    # plt.xticks([],[])
    # plt.yticks([],[])
    plt.text(2,5,flabels[ifolder],{'color':'white','fontsize':fs})

    subfig += 1
if savefig:
    savename = 'fig_pressureXY.eps'
    plt.savefig(savename, format='eps')



for plottime in range(len(T)):
    ff = 0
    ifig += 1

    # crg = [[2e6,2.5e7],[0.4,0.9],[0.1,0.6],[0,0.4]]
    crg = [[4,32],[0.0,0.2],[0.8,1.0],[0,0.2]]
    ilabel = -1
    for field in fields:
        ilabel += 1
        # print(np.shape(data3D[folders[ifolder]]['P']))
        for ifolder in range(len(folders)):
            plt.figure(ifig,figsize=(8,2.8))
            # slice = np.transpose(data3D[folders[ifolder]][field][:,30,:,plottime])
            slice = np.transpose(data3D[folders[ifolder]][field][:,:,Z[0],plottime])
            print(np.shape(data3D[folders[ifolder]][field]))
            if field == 'P':
                slice = slice / 1e6
            #     slice0 = data3D[folders[ifolder]][field][:,:,21,plottime]
            #     # slice0 = data3D[folders[ifolder]][field][:,:,40,plottime]
            #     slice = np.zeros((np.shape(slice0)[0]*2, np.shape(slice0)[1]))
            #     for ii in range(np.shape(slice0)[0]):
            #         for jj in range(np.shape(slice0)[1]):
            #             slice[ii*2:(ii+1)*2,jj] = slice0[ii,jj]
            #     slice = np.transpose(slice)
            # else:
            #     slice = np.transpose(data3D[folders[ifolder]][field][Z[0],:,:,plottime])
            plt.imshow(slice, vmin=crg[ff][0], vmax=crg[ff][1], cmap='jet')
            # plt.imshow(slice, vmin=crg[ff][0], vmax=crg[ff][1], cmap='jet', interpolation='bilinear')
            # plt.contourf(slice, vmin=crg[ff][0], vmax=crg[ff][1], cmap='jet')
            plt.colorbar()

            # if field == 'P':
            #     plt.ylabel('X [m]',fontsize=fs)
            #     plt.xlabel('Y [m]',fontsize=fs)
            #     xlabel = [-48,-32,-16,0,16,32,48]
            #     xtick = np.linspace(0,dim[1]*2-1,7)
            #     ylabel = [-354,-177,0,177,354]
            #     ytick = np.linspace(0,dim[0]-1,5)
            #     plt.title(flabels[ifolder],fontsize=fs)
            # else:
            plt.ylabel('Z [m]',fontsize=fs)
            plt.xlabel('X [m]',fontsize=fs)
            xlabel = [-354,-177,0,177,354]
            xtick = np.linspace(0,dim[0]-1,5)
            ylabel = [0,-46,-92,-138,-184]
            ytick = np.linspace(0,dim[2]-1,5)
            plt.xticks(xtick,xlabel)
            plt.yticks(ytick,ylabel)
            plt.title(figlabel[ilabel]+' ('+flabels[ifolder]+')',fontsize=fs)
            # if ifolder == 0:
            #     if field == 'P':
            #         plt.text(5,20,'Pressure [Pa]',{'color':'white','fontsize':fs})
            #     else:
            #         plt.text(5,10,figlabel[ilabel],{'color':'white','fontsize':fs})
            ifig += 1
            if savefig:
                savename = 'fig_'+flabels[ifolder]+'_'+field+'.eps'
                plt.savefig(savename, format='eps')
        ff += 1
plt.show()
