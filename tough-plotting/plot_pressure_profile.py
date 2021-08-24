""" Create contour plots from Plot_Data_Elem """
import numpy as np
import matplotlib.pyplot as plt
import os.path
from tough_mapping import *
# from matplotlib.colors import LogNorm
from matplotlib import colors
import matplotlib

#
# User settings
#
id = '1'
kernels = ['']
# folders = ['sens_paper/upslo_geos_minc2/','sens_paper/upslo_geos_minc2s0/','sens_paper/upslo_geos_minc2np/','sens_paper/upslo_geos_minc2blk/']
# folders = ['sens_paper/upslo2/','sens_paper/upslo2s0/','sens_paper/upslo2np/']
folders = ['minc09/upslo'+id+'s0/','minc09/upslo'+id+'eq15/','minc09/upslo'+id+'np2/']
# folders = ['sens_paper/ups09h1/']

flabels = ['LLP-'+id+' Mean','LLP-'+id+' $\lambda=1.5$','LLP-'+id+' unpropped']

# fgeos = '../scripts/meshmaker/GEOS08_upshi_prop_aper.csv'
fgeos = '../scripts/meshmaker/GEOS09_upslo_prop_aper.csv'

savefig = False

dim = [117,35,45]
offset = [354.0, 40.0, -120.0]
# dim = [117,35,85]
# offset = [354.0, 24.0, -280.0]

dx = [6.0, 0.002, 4.0]
dim_order = [0,2,1]


T = [8.64e7]
Z = [87]


fields = ['P']
figlabel = ['Pressure [MPa]']
fs = 8
lw = 1.0
ls = ['-','--','--','--',':',':']
co = ['b','r','g','k']
mk = ['o','+','x']

#
#   Load TOUGH results
#
elem, coord = readEleme(folders[0]+'MESH', dim)
eid, cid = getSliceID(elem, coord, [0, dim[0]*dx[0], 7.99, 8.01, -112, -108])
# eid, cid = getSliceID(elem, coord, [350,354, 7.99, 8.01, -280, -1])
sliceout = getSlice(folders, eid, fields, dim, T)

#
#   Load GEOS aperture
#
aper = np.genfromtxt(fgeos, delimiter=',')
aper2D = np.flipud(np.transpose(from1Dto2D(aper, dim, dx, offset, 8.0)))
print(np.shape(aper2D))

#
#   Make plot
#
xvec = np.linspace(0.5*dx[0], dim[0]*dx[0]-0.5*dx[0], dim[0])
cm = 1.0 / 2.54
font = {'size'   : 8}
matplotlib.rc('font', **font)
xlabel = [0, 175, 351, 526, 702]
xtick = np.linspace(0,dim[0]-1,5)
ylabel = [0, -40, -80, -120, -160]
ytick = np.linspace(0,40,5)


fig, (ax0, ax1, ax2) = plt.subplots(nrows=3,figsize=(12*cm,11*cm),
                  gridspec_kw={"height_ratios":[0.05, 1, 1]})

# Aperture
img = ax1.imshow(aper2D[:,:], vmin=0.0, vmax=1.5e-3,
        interpolation='bilinear', cmap='jet')
ax1.set_ylabel('Z [m]')
ax1.set_xticks(xtick)
ax1.set_xticklabels(xlabel)
ax1.set_yticks(ytick)
ax1.set_yticklabels(ylabel)
ax1.annotate(text='', xy=(0,28), xytext=(117,28), arrowprops=dict(arrowstyle='<->',color='r'),color='r')


# Pressure
for ifolder in range(len(folders)):
    slice = np.array(sliceout[folders[ifolder]]['P'])
    ax2.plot(slice / 1e6, color=co[ifolder], linestyle=ls[ifolder], linewidth=lw)
ax2.set_xlim([0,117])
ax2.set_xticks(xtick)
ax2.set_xticklabels(xlabel)
ax2.set_xlabel('X [m]', fontsize=fs)
ax2.set_ylabel('Pressure [MPa]',fontsize=fs)
ax2.legend(flabels, fontsize=fs)


fig.colorbar(img, cax=ax0, orientation="horizontal")
ax0.set_title('Aperture [m]')


if savefig:
    savename = 'Figure_12a.eps'
    plt.savefig(savename, format='eps')


plt.show()
