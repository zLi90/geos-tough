""" Functions used to convert GEOS to TOUGH """

import numpy as np
import silo_parser as sp

class geos_convert():

    def __init__(self, fgeos, fdir, y_frac, nproc=16):
        self.fgeos = fgeos
        self.fdir = fdir
        self.yVec = y_frac
        self.fields = ['Aperture','FaceArea','Pressure','ProppantVolumeFraction',
                  'Volume','birthTime',
                  'permeability','stressNOnFace','totalLeakedThickness']
        self.data = sp.parse_file(self.fgeos, 'FaceFields', self.fields, self.fdir, nproc)

    def get_one_field(self, field):
        """Function to extract one field from .silo file
        It creates a .csv file with columns [x, y, z, value]
        """
        print('   >>> Extracting ', field)
        out = np.r_['1,2,0',self.data['FaceCenter[0]'],self.data['FaceCenter[1]']
            ,self.data['FaceCenter[2]'],self.data[field]]
        return np.array(out)

    def get_actv_map(self, data1d):
        """ Get a map that only include active elements (remove nans) """
        map = np.zeros(np.shape(data1d)[0], dtype=int)
        for ii in range(np.shape(data1d)[0]):
            if not np.isnan(data1d[ii,3]):
                map[ii] = 1
        return map

    def apply_actv_map(self, data1d, map):
        out = []
        for ii in range(len(map)):
            if map[ii] == 1:
                out.append([np.round(data1d[ii,0]), np.round(data1d[ii,1]),
                    np.round(data1d[ii,2]), data1d[ii,3]])
        return np.array(out)

    def get_coordinates(self, data1d):
        """ Get coordinates (x,y,z) of GEOS data """
        coord = {}
        coord['xx'] = np.unique(np.round(data1d[:,0]))
        coord['yy'] = np.unique(np.round(data1d[:,1]))
        coord['zz'] = np.unique(np.round(data1d[:,2]))
        return coord

    def update_axis(self, ax, dx):
        """ Get coordinates with fixed grid resolution """
        new_ax = [ax[0]]
        while True:
            new_ax.append(new_ax[-1] + dx)
            if new_ax[-1] > ax[-1]:
                break
        return np.array(new_ax)

    def interp_z(self, data, ax, new_ax):
        """ Interpolation from non-uniform to uniform axis """
        add_rows = []
        for ii in range(np.shape(data)[0]):
            z = data[ii,2]
            val = np.amin(np.abs(new_ax - z))
            if val > 1e-3:
                ind = (np.abs(ax - z)).argmin()
                z_curr = ax[ind]
                z_next = ax[ind]
                z_prev = ax[ind]
                if ind < len(ax)-1:
                    z_next = ax[ind+1]
                if ind > 0:
                    z_prev = ax[ind-1]
                # distribute coarse reso values onto fine reso grids
                for jj in range(len(new_ax)):
                    if new_ax[jj]-z_curr >= 0 and new_ax[jj]-z_curr <= 0.5*(z_next-z_curr):
                        add_rows.append([data[ii,0], data[ii,1], new_ax[jj], data[ii,3]])
                    elif new_ax[jj]-z_curr <= 0 and z_curr-new_ax[jj] < 0.5*(z_curr-z_prev):
                        add_rows.append([data[ii,0], data[ii,1], new_ax[jj], data[ii,3]])
            else:
                add_rows.append([data[ii,0], data[ii,1], data[ii,2], data[ii,3]])
        return np.array(add_rows)


    def prop_aper(self, aper, prop, min_aper=5e-5):
        """ Calculate propped aperture """
        max_prop = np.nanmax(prop[:,3])
        print("   >>> Maximum proppant volume fraction = ",max_prop)
        prop_aper = []
        for ii in range(np.shape(aper)[0]):
            val = prop[ii,3] * aper[ii,3] / max_prop
            if val < min_aper:
                prop_aper.append([aper[ii,0], aper[ii,1], aper[ii,2], min_aper])
            else:
                prop_aper.append([aper[ii,0], aper[ii,1], aper[ii,2], val])
        return np.array(prop_aper)

    def convert_1D_to_3D(self, data1d, xlim, dx, offset):
        """ Convert 1D GEOS output to 3D (x,y,z) data just for plotting """
        # get coordinates based on the user-defined TOUGH domain
        coord = {}
        coord['xx'] = np.linspace(xlim[0], xlim[1], int((xlim[1]-xlim[0])/dx[0]+1))
        coord['yy'] = self.yVec
        coord['zz'] = np.linspace(xlim[4], xlim[5], int((xlim[5]-xlim[4])/dx[2]+1))
        dim = [len(coord['xx']),len(coord['yy']),len(coord['zz'])]
        geos3d = np.nan * np.ones((dim))
        print('   >>> Dimension of 3D GEOS data = ',dim)

        N = np.shape(data1d)[0]
        keys = ['xx','yy','zz']
        for ii in range(N):
            ind = []
            for jj in range(3):
                loc = np.round(data1d[ii,jj]) + offset[jj]
                ind.append((np.abs(coord[keys[jj]] - loc)).argmin())
            geos3d[ind[0], ind[1], ind[2]] = data1d[ii,3]
        return geos3d, coord
