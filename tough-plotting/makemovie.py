"""
    Create movie for TOUGH simulations, ZhiLi 20210512
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path
from tough_mapping import *
from matplotlib import colors
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

"""
    User settings for the simulation
"""
#   Directory of simulation results
directory = ['sens_ntwk2/upslo2_k1e-7_r1000_rp/']
#   Variable to be filmed
field = ['S_aqu']
#   Save times in the simulation
times = [8.64e4,8.64e5,8.64e6,8.64e7,1.296e8]
#   Dimension of computational domain and grid resolutions (x, y, z)
#   NOTE: Since dy is variable, adding y-axis labels could be troublesome
dim = [117,35,45]
dx = [6.0,0.002,4.0]
#   Rank the dimensions from shortest to longest
#   For example, if dim=[30,20,10], since 30>20>10, then dim_order=[2,1,0]
dim_order = [1,2,0]
#   Coordinate to extract output slice
#   The order of cutpoint is from the longest to the shortest (Apologize for the confusion!)
#   For example, if dim=[30,20,10], then cutpoint=[None,None,0] will extract the
#   x-y plane at z=0.
#   NOTE: Unused coordinate must be 'None'
#   For LLP and HLP scenarios, cutpoint=[x, z, y]
cutpoint = [None, None, 17]
#   Whether or not to save the output movie
saveani = True
savename = 'anime.mp4'

"""
    User settings for plotting
"""
#   Figure caption and variable label
caption = ['LLP-2']
flabel = ['$s_{g}$']
#   Range of colorbar
crg = [0.0,0.7]
# crg = [1e7,2.5e7]
#   Figure size in cm
figsz = [40,15]
#   Font size
fs = 16
#   Axis labels (x, y)
ax_labels = ['X [m]', 'Z [m]']
#   Add tick labels or not
#   NOTE: Adding y-labels should be careful due to non-uniform resolution
add_tick_labels = True
xtick = np.linspace(0, dim[0], 5)
ytick = np.linspace(0, dim[2], 5)
xtick_label = [0,177,352,527,702]
ytick_label = [0,-46,-92,-138,-184]


"""
    Executing code
"""
data3D = get3Ddata(directory, dim, dim_order, field, times)
xyz = ['x','y','z']
print(' >>> Order of output dimension is ',xyz[dim_order[2]],xyz[dim_order[1]],xyz[dim_order[0]])
print(' >>> Dimension of output = ',np.shape(data3D[directory[0]][field[0]]))

cm = 1.0 / 2.54
imgs = []
fig = plt.figure(1,figsize=(figsz[0]*cm,figsz[1]*cm))
for t in range(len(times)):
    if cutpoint[0] != None:
        slice = np.transpose(data3D[directory[0]][field[0]][cutpoint[0],:,:,t])
    elif cutpoint[1] != None:
        slice = np.transpose(data3D[directory[0]][field[0]][:,cutpoint[1],:,t])
    elif cutpoint[2] != None:
        slice = np.transpose(data3D[directory[0]][field[0]][:,:,cutpoint[2],t])
    img = plt.imshow(slice, vmin=crg[0], vmax=crg[1], cmap='jet', animated=True)
    txt = plt.annotate(str(times[t]/86400)+' Days', (5,5), color='w', fontsize=fs)
    imgs.append([img, txt])
plt.colorbar()
plt.title(caption[0]+' '+flabel[0],fontsize=fs)
plt.xlabel(ax_labels[0],fontsize=fs)
plt.ylabel(ax_labels[1],fontsize=fs)
if add_tick_labels:
    plt.xticks(xtick,xtick_label)
    plt.yticks(ytick,ytick_label)

ani = animation.ArtistAnimation(fig, imgs, interval=500, blit=True)

if saveani:
    ani.save(savename)

plt.show()
