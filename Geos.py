""" This class is used to process GEOS output data to match TOUGH domain """
import numpy as np
import numpy.matlib as nplib
import matplotlib.pyplot as plt
import math
import copy
from scipy.interpolate import griddata


class Geos():
    """ Each instance of the Geos class represents one GEOS variable (e.g. Aperture,
        Pressure, etc.). The variable is originally store in an .csv file. The methods
        of the Geos class are used to read the .csv file, extract the data that fits
        the TOUGH domain, scale/smooth the data, and interpolate to different resolutions
        to meet requirements of any specific simulation scenario.

    Args:
        fgeos (dict)    : Input dict containing basic geos info. It has the format:
                        [fname, [savename, smooth, scale], bcutoff]
                        fname   : Name of the GEOS output file of the variable
                        savename: Name of the new GEOS file after processing
                        smooth  : If True, average the GEOS data to remove heterogeneity
                        scale   : Scale the magnitude of the GEOS data
                        bcutoff : Remove GEOS data outside the range of bcutoff
        geometry (dict) : Geometry information of the model domain. It has the keys:
                        delta   : Baseline element size in [x, y, z] directions
                        dim     : Number of elements in [x, y, z] directions
                        offset  : Distance between GEOS coordinates and TOUGH well
                                  locations. len(offset) = number of wells
                        ratio   : Scaling factor to interpolate GEOS data onto TOUGH mesh
                        yfrac   : Y-coordinates of the fractures in the GEOS data
                        bound   : Boundaries ('x+','z-',etc.) that will be marked 'I'

    Attributes:
        fname (str)     : Same as fgeos[variable][fname]
        delta (list)    : Same as geometry['delta']
        yfrac (list)    : Same as geometry['yfrac']
        offset (list)   : Same as geometry['offset']
        ratio (list)    : Same as geometry['ratio']
        geos (nparray)  : GEOS output read from 'fname'
        N (int)         : Length of GEOS data
        range (list)    : [min, max] of GEOS coordinates
        xVec (nparray)  : X-coordinates of GEOS data, the default GEOS grid resolution
                          is 4m
        zVec (nparray)  : Z-coordinates of GEOS data
        dim (list)      : Dimension of GEOS fracture plane, [len(xVec), len(zVec)]

    """
    def __init__(self, fgeos, geometry):
        self.fname = fgeos[0]
        self.delta = geometry['delta']
        self.yfrac = geometry['yfrac']
        self.offset = geometry['offset']
        self.ratio = geometry['ratio']
        self.geos = np.genfromtxt(self.fname, delimiter=',')
        if len(fgeos[2]) > 0:
            for ii in range(3):
                aa = abs(self.geos[:,ii]) < fgeos[2][ii]
                self.geos = self.geos[aa,:]
        self.N = self.geos.shape[0]

    def getRange(self):
        """ Get range of GEOS data in [xmin, xmax, zmin, zmax]"""
        # get range of GEOS data
        self.range = []
        self.range.append(np.amin(self.geos[:,0]))
        self.range.append(np.amax(self.geos[:,0]))
        self.range.append(np.amin(self.geos[:,2]))
        self.range.append(np.amax(self.geos[:,2]))
        # create coordinates along each axis
        self.xVec = np.linspace(self.range[0], self.range[1], int((self.range[1]-self.range[0])/2.0+1))
        self.zVec = np.linspace(self.range[2], self.range[3], int((self.range[3]-self.range[2])/2.0+1))
        # dimension of 2D GEOS domian
        self.dim = [len(self.xVec), len(self.zVec)]
        print(self.xVec)
        print(self.zVec)

    def from1Dto2D(self):
        """ Convert 1D GEOS output into 2D numpy matrix """
        data2D = np.nan * np.ones((self.dim[0],len(self.yfrac),self.dim[1],4))
        for ll in range(self.N):
            xx = np.where(abs(self.geos[ll,0]-self.xVec) == np.amin(abs(self.geos[ll,0]-self.xVec)))
            zz = np.where(abs(self.geos[ll,2]-self.zVec) == np.amin(abs(self.geos[ll,2]-self.zVec)))
            yy = np.where(abs(self.geos[ll,1]-self.yfrac) == np.amin(abs(self.geos[ll,1]-self.yfrac)))
            for ii in range(4):
                data2D[xx,yy,zz,ii] = self.geos[ll,ii]
        return data2D

    def interp2D(self):
        """ Interpolate to match required resolution

        NOTE: Since GEOS data contains NANs, direct 2D interpolation onto finer
              resolution could be troublesome. For simplicity, only piecewise
              constant 'interpolation' is performed so far.
        """
        dim = [math.floor(self.dim[0]/self.ratio[0]), len(self.yfrac), math.floor(self.dim[1]/self.ratio[2]), 4]
        out = np.nan * np.ones(dim)
        # loop over all fractures
        for ll in range(dim[1]):
            # get the GEOS data for this fracture
            aa = abs(self.geos[:,1]-self.yfrac[ll]) < self.delta[1]
            xVec = self.geos[aa,0].tolist()
            zVec = self.geos[aa,2].tolist()
            aper = self.geos[aa,3].tolist()
            coord = np.zeros((dim[0]*dim[2],2))
            # compute new coordinates
            xNew = np.zeros((dim[0],1))
            zNew = np.zeros((dim[2],1))
            kk = 0
            for ii in range(dim[0]):
                xNew[ii] = self.range[0] + ii*self.delta[0]
                for jj in range(dim[2]):
                    if ii == 0:
                        zNew[jj] = self.range[2] + jj*self.delta[2]
                    coord[kk,0] = self.range[0] + ii*self.delta[0]
                    coord[kk,1] = self.range[2] + jj*self.delta[2]
                    kk += 1
            # interpolate GEOS data to match new coordinates
            aper2D = griddata((xVec, zVec), aper, (coord[:,0],coord[:,1]), method='linear')
            out[:,ll,:,3] = np.reshape(aper2D,(dim[0],dim[2]))
            out[:,ll,:,0] = nplib.repmat(xNew,1,dim[2])
            out[:,ll,:,1] = self.yfrac[ll]
            out[:,ll,:,2] = np.transpose(nplib.repmat(zNew,1,dim[0]))

        return out

    def from2Dto1D(self, data2D):
        """ Convert 2D data back to 1D (N by 4) array """
        for ii in range(4):
            temp = np.reshape(data2D[:,:,:,ii],(np.size(data2D[:,:,:,ii])))
            if ii == 0:
                out1D = np.zeros((np.size(temp),4))
            out1D[:,ii] = temp
        for ii in range(np.shape(out1D)[0]):
            if np.isnan(out1D[ii,3]):
                out1D[ii,3] = 0.0
        return out1D

    def shift(self, out1D, offset):
        """ Shift the coordinates to match the TOUGH well location """
        temp = np.zeros((out1D.shape))
        for ii in range(3):
            temp[:,ii] = out1D[:,ii] + offset[ii]
        temp[:,3] = out1D[:,3]
        return temp

    def duplicate(self, out1D):
        """ Duplicate (and shift) to match multi-well scenarios """
        for kk in range(len(self.offset)):
            temp = self.shift(out1D, self.offset[kk])
            if kk == 0:
                out = temp
            else:
                out = np.vstack((out,temp))
        out[:,1] = np.round(out[:,1])
        return out

    def scale(self, data, scale, min_aper=1e-5, max_aper=5e-2):
        """ Scale geos data """
        targ_avg = scale[1]
        targ_std = scale[2]
        out = copy.deepcopy(data)
        # out[:,3] = data[:,3] * scale

        #   Adjust aperture to match desired mean and std
        if targ_avg != None:
            actv = out[:,3] > 0.0
            new_mean = np.mean(out[actv,3])
            iter = 0
            while abs(new_mean - targ_avg)/targ_avg > 0.05:
                for ii in range(np.shape(out)[0]):
                    if actv[ii] == 1:
                        out[ii,3] = data[ii,3] - (np.mean(data[actv,3]) - targ_avg)
                new_mean = np.mean(out[actv,3])
                if targ_std != None:
                    coeff = targ_std / np.std(data[actv,3])
                    out[actv,3] = (out[actv,3] - new_mean) * coeff + new_mean
                # Iterate to remove negative aperture
                for ii in range(np.shape(out)[0]):
                    if actv[ii] == 1 and out[ii,3] <= 0.0:
                        out[ii,3] = min_aper
                new_mean = np.mean(out[actv,3])
                print('Scaling GEOS: iter ',iter,' aperture has mean, std = ',np.mean(out[actv,3]),np.std(out[actv,3]))
                iter += 1
                if iter > 10:
                    break
        elif targ_std != None:
            actv = out[:,3] > 0.0
            new_mean = np.mean(out[actv,3])
            coeff = targ_std / np.std(data[actv,3])
            out[actv,3] = (out[actv,3] - new_mean) * coeff + new_mean
            for ii in range(np.shape(out)[0]):
                if actv[ii] == 1 and out[ii,3] <= 0.0:
                    out[ii,3] = min_aper
            new_mean = np.mean(out[actv,3])
            print('Scaling GEOS: aperture has mean, std = ',np.mean(out[actv,3]),np.std(out[actv,3]))
        #   remove large and small aperture
        aa = out[:,3] == 0.0
        out[aa,3] = np.nan
        aa = out[:,3] > max_aper
        out[aa,3] = max_aper
        aa = out[:,3] < min_aper
        out[aa,3] = min_aper
        return out

    def smooth(self, data):
        """ Use constant (mean) geos data """
        out = copy.deepcopy(data)
        out[:,3] = np.nanmean(data[:,3])
        return out
