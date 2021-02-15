""" Plot the GEOS aperture """
import numpy as np
import matplotlib.pyplot as plt

fname = 'GEOS_upshi_A.csv'
savename = 'geos_5cluster_upshi.csv'
mirror = False
geos = np.genfromtxt(fname, delimiter=',')

# Convert to 3D fields
N = np.shape(geos)[0]
for ii in range(N):
    for jj in range(3):
        geos[ii,jj] = np.round(geos[ii,jj])

xx = np.unique(geos[:,0])
yy = np.unique(geos[:,1])
zz = np.unique(geos[:,2])
print(xx)
print(yy)
print(zz)
dim = [len(xx),len(yy),len(zz)]
geos3d = np.zeros((dim))
print('Dimension of GEOS data = ',dim)


for ll in range(N):
    ii = np.where(xx == geos[ll,0])
    jj = np.where(yy == geos[ll,1])
    kk = np.where(zz == geos[ll,2])
    geos3d[ii,jj,kk] = geos[ll,3]

# Shrink to the real domain size
geosOut = geos3d
# geosOut = geos3d[30:120,:,5:45]
dimOut = np.shape(geosOut)
print('Final output domain size = ',np.shape(geosOut))

# Print maximum aperture
for jj in range(len(yy)):
    slice = geosOut[:,jj,:]
    print('Slice ',jj,' - Max aper = ',np.amax(slice))

# Get coordinates
geos1d = np.zeros((dimOut[0]*dimOut[1]*dimOut[2],4))
xVec = np.linspace(xx[0],xx[-1],len(xx))
yVec = yy
# zVec = np.linspace(-dimOut[2]*2.0+1,-1.0,dimOut[2]) - 10.0
zVec = np.linspace(zz[0],zz[-1],len(zz))
# zVec -= 2.0
ll = 0
for jj in range(len(yVec)):
    for ii in range(len(xVec)):
        for kk in range(len(zVec)):
            geos1d[ll,0] = xVec[ii]
            geos1d[ll,1] = yVec[jj]
            geos1d[ll,2] = zVec[kk]
            geos1d[ll,3] = geosOut[ii,jj,kk]
            ll += 1

# Mirror aperture
if mirror:
    N = np.shape(geos1d)[0]
    offset = [-np.amin(geos1d[:,0])+1.0, 0.0, -np.amin(geos1d[:,2])+1.0, 0.0]
    print(offset)
    out = np.zeros((2*N,4))
    # copy original aperture
    for ii in range(N):
        for jj in range(4):
            out[ii,jj] = geos1d[ii,jj] + offset[jj]
    # copy mirrored aperture
    for ii in range(N):
        for jj in range(4):
            out[ii+N,jj] = geos1d[ii,jj] + offset[jj]
            if jj == 2:
                out[ii+N,jj] *= -1
    # final adjust of offset
    out[:,2] -= offset[2]
    geos1d = out
    # get 3D geos field
    N = np.shape(geos1d)[0]
    xx = np.sort(np.unique(geos1d[:,0]))
    yy = np.sort(np.unique(geos1d[:,1]))
    zz = np.sort(np.unique(geos1d[:,2]))
    dim = [len(xx),len(yy),len(zz)]
    geos3d = np.zeros((dim))
    print(zz)
    print('Dimension of GEOS data = ',dim)
    for ll in range(N):
        ii = np.where(xx == geos1d[ll,0])
        jj = np.where(yy == geos1d[ll,1])
        kk = np.where(zz == geos1d[ll,2])
        geos3d[ii,jj,kk] = geos1d[ll,3]
    geosOut = geos3d

# Save geos aperture
# np.savetxt(savename, geos1d, delimiter=",")

# Make plot
plt.figure(1,figsize=(10,6))

for ifig in range(dim[1]):
    # plt.subplot(2,3,ifig+1)
    # if ifig < 3:
    plt.subplot(3,2,ifig+1)
    plt.imshow(np.flipud(np.transpose(geosOut[:,ifig,:])), cmap='jet', vmin=0.0, vmax=2e-2)
    plt.colorbar()
# plt.imshow(np.flipud(np.transpose(geosOut[:,2,:])), cmap='jet', vmin=0.0, vmax=2e-2)


plt.show()

# Interpolate to desired range and resolution
