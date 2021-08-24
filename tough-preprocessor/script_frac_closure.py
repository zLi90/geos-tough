""" Estimate fracture closure following Morris2016 """
import numpy as np

"""
    Input parameters
"""
#   GEOS output file
faper = 'GEOS09_upslo_prop_aper.csv'
fstre = 'GEOS09_upslo_stre.csv'

#   Output file name
save_output = True
fsave = 'GEOS09_upslo_clos_aper.csv'

#   Normal stress [Pa]
sigma = np.genfromtxt(fstre, delimiter=',')

#   Initial aperture [m]
b = np.genfromtxt(faper, delimiter=',')

#   Young's modulus [Pa]
E = 1e11

#   Initial guess of final apeture [m]
binit = 1e-4


"""
    Use iterative method to estimate fracture closure
"""
#   Loop to estimate fracture closure
aper = []
tol = 1e-7
print(' \n >>>>> Begin computation! Max, Min, Mean, Std aperture = ',
    np.nanmax(b[:,3]), np.nanmin(b[:,3]), np.nanmean(b[:,3]), np.nanstd(b[:,3]))
for ii in range(np.shape(sigma)[0]):
    if not np.isnan(b[ii,3]):
        b0 = binit
        eps = 1.0
        iter = 0
        while eps > tol:
            b1 = b[ii,3] - abs(sigma[ii,3])*b[ii,3]/(E*b0 + abs(sigma[ii,3]))
            eps = abs(b1 - b0)
            b0 = np.maximum(b1, 0.0)
            iter += 1
            if iter > 100:
                print(' >>> Iteration cannot converge for [x, y, z, b] = ',b[ii,:],', loop breaks!')
                break
            elif b0 == 0.0:
                print(' >>> Aperture goes negative for [x, y, z, b] = ',b[ii,:],', loop breaks!')
                b1 = 0.0
                break
        b1 = np.maximum(b1, 0.0)
        aper.append([b[ii,0], b[ii,1], b[ii,2], b1])
        if b1 > b[ii,3]:
            print(' >>> Aperture increases to ', b1,' for [x, y, z, b] = ',b[ii,:])
aper = np.array(aper)

print(' >>>>> Computation completed! Max, Min, Mean, Std aperture = ',
    np.nanmax(aper[:,3]), np.nanmin(aper[:,3]), np.nanmean(aper[:,3]), np.nanstd(aper[:,3]))

if save_output:
    np.savetxt(fsave, aper, delimiter=',')
