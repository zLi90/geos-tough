""" Python functions that write MeshMaker input """
import numpy as np
import subprocess as sp
import copy
import os
import math

"""
    ============================================================================
                    Class 'Mesh' for generating TOUGH MESH
    ============================================================================
"""
class Mesh():
    """ A 'Mesh' object controls the initialization and modification of a
        TOUGH mesh.

    Assumptions:
        - Cartesian mesh with 'old' output format
        - Length unit is 'm'
        - Origin at (0,0,0)
        - Media defined by number (to allow use of MINC)

    Args:
        dim (list)      : Number of elements in each direction, [nx, ny, nz]
        region (list)   : Defines the regions (mediums)
                region = [num_region, default_medium, other_medium]
                    other_medium = [[name,x_min,x_max,y_min,y_max,z_min,z_max],[],...]
        boundary (list) : Defines the boundaries (e.g. wells)
                boundary = [num_boundaries, boundaries]
                    boundaries = [name,x_min,x_max,y_min,y_max,z_min,z_max]
        discret (dict)  : Defines the discretizations
                discret = {disc_x, disc_y, disc_z}
                    disc_x = [num_subset, disc_x1, disc_x2, ...]
                        disc_x1 = {num_dx, opt, dx, dx_max, dx_log}
        minc (bool)     : Whether or not use MINC
        remesh (bool)   : Create original mesh or modify base on an existing mesh

    Attributes:


    """
    def __init__(self, dim, delta, region, boundary, discret, minc=[], remesh=False):
        self.nx = dim[0]
        self.ny = dim[1]
        self.nz = dim[2]
        self.dx = delta[0]
        self.dy = delta[1]
        self.dz = delta[2]
        self.region = region
        self.boundary = boundary
        self.disc = discret
        self.minc = minc
        self.remesh = remesh

    ''' <<<<<<<< Some useful utility functions >>>>>>>>>> '''

    def findInd(self, value, array, delta):
        """ Find the index of element in array who is closest to value

        Args:
            value (float): Value to be put into array for searching
            array (nparray): The array of all possible values

        Returns:
            Index of value in array. If value is out of the range of array, return nan.
        """
        if value >= np.amin(array)-0.5*delta and value <= np.amax(array)+0.5*delta:
            ind = np.where(abs(value - array) == np.amin(abs(value - array)))
            return np.where(abs(value - array) == np.amin(abs(value - array)))
        else:
            return np.nan


    def reformat(self, strn):
        """ Make python formatE starts with zero

        Args:
            strn (str): A number in Python E format, e.g., 1.2345E+01

        Returns:
            stro (str): Same number but starts with zero, e.g., 0.1234E+02
        """
        if float(strn) != 0.0:
            if strn[-3] == '+':
                if strn[-2] == '0':
                    expo = 'E+0'+str(int(strn[-1]) + 1)
                else:
                    expo = 'E+'+str(int(strn[-2:]) + 1)
            elif strn[-3] == '-':
                if strn[-2] == '0':
                    expo = 'E-0'+str(int(strn[-1]) - 1)
                else:
                    if strn[-1] == '0':
                        expo = 'E-'+str(int(strn[-2])-1)+'9'
                    else:
                        expo = 'E-'+str(int(strn[-2:]) - 1)
            if float(strn) > 0.0:
                sign = '0.'
            else:
                sign = '-.'
            digits = str(round(float(strn)/float('1'+strn[-4:])*1000))[-4:]
            stro = f"{sign}{digits}{expo}"
        else:
            stro = strn
        return stro

    ''' <<<<<<<< End of utility functions >>>>>>>>>> '''

    def build_mesh(self, fname='mesh_input'):
        """ Calls MeshMaker to generate 'MESH' from 'mesh_input' """
        with open(fname,'r') as fid:
            fmtstr = fid.read()
        proc = sp.Popen(['./mm3'], stdin=sp.PIPE)
        proc.communicate(str.encode(fmtstr))
        return 0

    def load_mesh(self, fname):
        """ Load an existing MESH file for modification

        Returns:
            elem (list): List of all elements in ELEME block
                         Format: [elem id, media, volume, area, x, y, z, isbound]
            conn (list): List of all connections in CONNE block
                         Format: [conn id1, id2, axis(x,y,z), l1, l2, area]
        """
        fm = open(fname, 'r')
        elem = []
        conn = []
        while True:
            line = fm.readline()
            if line[0:5] == 'CONNE':
                while True:
                    cline = fm.readline()
                    if cline[0:5] == '     ':
                        break
                    else:
                        # conn has structure (conn id1, id2, axis, d1, d2, area)
                        conn.append([cline[0:5], cline[5:10], int(cline[29]),
                            float(cline[30:40]), float(cline[40:50]), float(cline[50:60]), cline[60:70]])
                break
            elif line[0:5] != 'ELEME' and line[0:5] != '     ' and line[5] == ' ':
                # elem has the structure (element id, media, V, A, x, y, z, isbound)
                isbound = 0
                if line[81:82] == 'I':
                    isbound = 1
                elem.append([line[0:5], line[15:20], float(line[20:30]), float(line[30:40]),
                    float(line[50:60]), float(line[60:70]), float(line[70:80]), isbound])
        fm.close()
        self.Nele = len(elem)
        self.Ncon = len(conn)
        # read (x,y,z) from elem
        self.xVec = np.sort(np.unique(np.array([item[4] for item in elem])))
        self.yVec = np.sort(np.unique(np.array([item[5] for item in elem])))
        self.zVec = np.sort(np.unique(np.array([item[6] for item in elem])))
        self.zVec[::-1].sort()
        return elem, conn

    def get_boundary_id(self, elem, conn, map):
        """ Get the connection ID for the well """
        allid = []
        dir = ['iP','iM','kP','kM']
        for ii in range(self.Nele):
            if elem[ii][7] == 1:
                eid = elem[ii][0]
                for jj in range(self.Ncon):
                    if eid == conn[jj][0] and conn[jj][2] != 2:
                        allid.append(conn[jj][0]+conn[jj][1])
                        # for kk in range(len(dir)):
                        #     if conn[jj][1] == elem[map[dir[kk]][ii,0]][0] and elem[map[dir[kk]][ii,0]][1] != elem[ii][1]:
                        #         allid.append(conn[jj][0]+conn[jj][1])
                    elif eid == conn[jj][1] and conn[jj][2] != 2:
                        allid.append(conn[jj][0]+conn[jj][1])
                        # for kk in range(len(dir)):
                        #     if conn[jj][0] == elem[map[dir[kk]][ii,0]][0] and elem[map[dir[kk]][ii,0]][1] != elem[ii][1]:
                        #         allid.append(conn[jj][0]+conn[jj][1])
        return allid

    def write_input(self, use_minc):
        """ Write MeshMaker input file named 'mesh_input' """
        finput = 'mesh_input'
        fid = open(finput, 'w')
        # Write general info
        fid.write("Input file for MeshMaker\n")
        fid.write(">>>GENERAL_INFO\n")
        fid.write(f"&Grid_Specifications coordinate_system  = 'Cartesian',\n")
        fid.write(f"                     output_file_format = 'old',\n")
        fid.write(f"                     length_units        = 'm',\n")
        fid.write(f"                     grid_numbering_system = 'standard',\n")
        fid.write(f"                     ElemName_NumCharacters = 5,\n")
        fid.write(f"                     MaxNum_X_Subdivisions = {self.nx},\n")
        fid.write(f"                     MaxNum_Y_Subdivisions = {self.ny},\n")
        fid.write(f"                     MaxNum_Z_Subdivisions = {self.nz},\n")
        fid.write(f"                     AxesOrigin_X          = 0.0d0,\n")
        fid.write(f"                     AxesOrigin_Y          = 0.0d0,\n")
        fid.write(f"                     AxesOrigin_Z          = 0.0d0,\n")
        fid.write(f"                     inclination_angle     = 0.0,\n")
        fid.write(f"                     areas_for_HeatExch_Solution = .TRUE.,\n")
        if use_minc:
            fid.write(f"                     media_by_number       = .TRUE.,\n")
        else:
            fid.write(f"                     media_by_number       = .FALSE.,\n")
        fid.write(f"                     special_features      = .FALSE.\n")
        fid.write("      /\n")
        fid.write("<<<\n\n")
        # Write MINC block
        if len(self.minc) > 0:
            fid.write(">>>MINC Data Block\n")
            fid.write(f"&General_MINC_Data number_of_media                = {self.minc[0]},\n")
            fid.write(f"                   number_of_fractured_media      = {self.minc[1]},\n")
            fid.write(f"                   number_of_MINCs                = {self.minc[2]},\n")
            fid.write(f"                   number_of_specified_VolFractions = {self.minc[3]},\n")
            fid.write(f"                   fracture_is_1st_continuum      = .TRUE.,\n")
            fid.write(f"                   matrix_to_matrix_flow          = '{self.minc[4]}',\n")
            fid.write("      /\n")
            for ii in range(self.minc[1]):
                minc_medium = self.minc[5+ii]
                fid.write(f"&Medium_MINC_Data  fractured_medium_number    = {minc_medium[0]},\n")
                fid.write(f"                   proximity_function_type    = 'THRED',\n")
                fid.write(f"                   fracture_DistParameters    = ")
                for jj in range(len(minc_medium[1])):
                    fid.write(f"{minc_medium[1][jj]}, ")
                fid.write("\n")
                fid.write(f"                   VolFraction_data_format    = '*'\n / \n")
                fid.write(f"{minc_medium[2]} \n")
            fid.write("<<<\n\n")
        # Write regions
        fid.write(">>>REGIONS\n")
        num_reg = 0
        for ii in range(1, len(self.region)):
            num_reg += len(self.region[ii])-1
        fid.write(f"&Heterogeneous_Regions number_of_regions = {num_reg},\n")
        fid.write(f"                       dominant_medium = '{self.region[1][0]}' /\n")
        for ii in range(self.region[0]-1):
            reg_list = self.region[ii+2]
            for jj in range(1, len(reg_list)):
                fid.write(f"    &HetRegion_GeneralInfo region_name  = '{reg_list[0]}',\n")
                fid.write(f"                           region_shape = 'Rectangular',\n")
                fid.write(f"                           length_units = 'm',\n")
                fid.write("     /\n")
                fid.write(f"    &Rectangular_HetRegion X_min = {reg_list[jj][0]}, Y_min = {reg_list[jj][2]}, Z_min = {reg_list[jj][4]},\n")
                fid.write(f"                           X_max = {reg_list[jj][1]}, Y_max = {reg_list[jj][3]}, Z_max = {reg_list[jj][5]},\n")
                fid.write("     /\n")
        fid.write("<<<\n\n")
        # Write boundaries
        fid.write(">>>BOUNDARIES\n")
        fid.write(f"&Boundary_Regions number_of_boundaries = {self.boundary[0]}/\n")
        for ii in range(self.boundary[0]):
            bd_list = self.boundary[ii+1]
            fid.write(f"    &Boundary_GeneralInfo boundary_name  = '{bd_list[0]}',\n")
            fid.write(f"                          boundary_shape = 'Rectangular',\n")
            fid.write(f"                          boundary_type  = 'I',\n")
            fid.write(f"                          length_units   = 'm',\n")
            fid.write("     /\n")
            fid.write(f"    &Rectangular_Boundary  X_min = {bd_list[1][0]}, Y_min = {bd_list[1][2]}, Z_min = {bd_list[1][4]},\n")
            fid.write(f"                           X_max = {bd_list[1][1]}, Y_max = {bd_list[1][3]}, Z_max = {bd_list[1][5]},\n")
            fid.write("     /\n")
        fid.write("<<<\n\n")
        # Discretization
        fid.write(">>>DISCRETIZATION\n")
        fid.write(f":::>>>X-Discretization   &DX_Subsets num_DX_subsets = {self.disc['disc_x'][0]} /\n")
        for ii in range(self.disc['disc_x'][0]):
            disc = self.disc['disc_x'][ii+1]
            fid.write(f"    &DX_data number_of_DXs = {disc[0]},\n")
            fid.write(f"             option        = '{disc[1]}',\n")
            fid.write(f"             Delta_X       = {disc[2]},\n")
            if disc[1][0:3] == 'Log':
                fid.write(f"             X_max         = {disc[3]},\n")
                fid.write(f"             DX_log        = {disc[4]},\n")
            fid.write("     /\n")
        fid.write(":::<<<\n")
        fid.write(f":::>>>Y-Discretization   &DY_Subsets num_DY_subsets = {self.disc['disc_y'][0]} /\n")
        for ii in range(self.disc['disc_y'][0]):
            disc = self.disc['disc_y'][ii+1]
            fid.write(f"    &DY_data number_of_DYs = {disc[0]},\n")
            fid.write(f"             option        = '{disc[1]}',\n")
            fid.write(f"             Delta_Y       = {disc[2]},\n")
            if disc[1][0:3] == 'Log':
                fid.write(f"             Y_max         = {disc[3]},\n")
                fid.write(f"             DY_log        = {disc[4]},\n")
            fid.write("     /\n")
        fid.write(":::<<<\n")
        fid.write(f":::>>>Z-Discretization   &DZ_Subsets num_DZ_subsets = {self.disc['disc_z'][0]} /\n")
        for ii in range(self.disc['disc_z'][0]):
            disc = self.disc['disc_z'][ii+1]
            fid.write(f"    &DZ_data number_of_DZs = {disc[0]},\n")
            fid.write(f"             option        = '{disc[1]}',\n")
            fid.write(f"             Delta_Z       = {disc[2]},\n")
            if disc[1][0:3] == 'Log':
                fid.write(f"             Z_max         = {disc[3]},\n")
                fid.write(f"             DZ_log        = {disc[4]},\n")
            fid.write("     /\n")
        fid.write(":::<<<\n")
        # End of input file
        fid.write("<<<\n\n")
        fid.close()
        return finput


    def neighbor_map(self, elem, conn):
        """ Build map among TOUGH elements

        Args:
            elem, conn (list): Obtained from read_mesh()

        Returns:
            map (dict): Dict with keys 'iP'(x-plus), 'iM'(x-minus), 'jP', 'jM', 'kP', 'kM'
                        e.g., map['iP'] = [iP elem, i-iP connection]
        """
        # initialize data arrays
        map = {}
        keys = ['iP','iM','jP','jM','kP','kM']
        for ii in range(len(keys)):
            map[keys[ii]] = -1 * np.ones((self.Nele,2),dtype=int)
        # save eid into a dict
        eids = {}
        for ii in range(self.Nele):
            eids[elem[ii][0]] = ii
        # loop over all connections
        for ii in range(self.Ncon):
            cid1 = conn[ii][0]
            cid2 = conn[ii][1]
            axis = conn[ii][2]
            # find index in elem
            ind1 = eids[cid1]
            ind2 = eids[cid2]
            # save to map
            map[keys[(axis-1)*2]][ind1,0] = ind2
            map[keys[(axis-1)*2]][ind1,1] = ii
            map[keys[(axis-1)*2+1]][ind2,0] = ind1
            map[keys[(axis-1)*2+1]][ind2,1] = ii
        return map

    def geos_map(self, elem, conn, geos, zone):
        """ Build map between GEOS and TOUGH mesh

        Args:
            elem, conn (list): Obtained from readMesh()
            geos (nparray)   : Aperture data from GEOS, [x, y, z, aperture]
            zone (list)      : Media names for fractures
        Returns:
            gmap (nparray)   : gmap[row index in 'geos'] = row index in 'elem'
        """
        dim = np.shape(geos)
        self.aper = geos
        # initialize map arrays
        gmap = -1 * np.ones((dim[0],1), dtype=int)
        isfrac = np.zeros((self.Nele, 1), dtype=bool)
        # find corresponding geos element in tough
        for ii in range(dim[0]):
            rowx = self.findInd(self.aper[ii,0], self.xVec, self.dx)
            rowy = self.findInd(self.aper[ii,1], self.yVec, self.dy)
            rowz = self.findInd(self.aper[ii,2], self.zVec, self.dz)
            # Note that the order of elements in TOUGH seems depend on the
            # dimension of the domain. So here we sort nx, Ny and Nz to get
            # the correct element index.
            if any([isinstance(rowx,float),isinstance(rowy,float),isinstance(rowz,float)]):
                gmap[ii] = -1
            else:
                rowarray = np.array([[self.nz, rowz[0]],[self.ny, rowy[0]],[self.nx, rowx[0]]])
                rows = rowarray[rowarray[:,0].argsort()]
                row = rows[2,1] * rows[0,0] * rows[1,0] + rows[1,1] * rows[0,0] + rows[0,1]
                if len(row) > 1:
                    row = int(row[0])
                else:
                    row = int(row)
                if elem[row][1] in zone:
                    if not np.isnan(self.aper[ii,3]):
                        gmap[ii] = row
                        isfrac[row] = True
        return gmap, isfrac

    def geos_map_minc(self, elem, conn, geos_aper, zone):
        """ Similar to geos_map, but returns the element id in elem """
        dim = np.shape(geos_aper)
        self.aper = geos_aper
        # initialize map arrays
        gmap = -1 * np.ones((dim[0],1), dtype=int)
        isfrac = np.zeros((self.Nele, 1), dtype=bool)
        # gmap = []
        # isfrac = []
        # get full coordinates
        self.xVec_full = np.array([item[4] for item in elem])
        self.yVec_full = np.array([item[5] for item in elem])
        self.zVec_full = np.array([item[6] for item in elem])
        # find corresponding geos element in tough
        for ii in range(len(elem)):
            if elem[ii][1] in zone:
                elem_x = elem[ii][4]
                elem_y = elem[ii][5]
                elem_z = elem[ii][6]
                for jj in range(dim[0]):
                    match_x = abs(elem_x - self.aper[jj,0]) < self.dx
                    match_y = abs(elem_y - self.aper[jj,1]) < self.dy
                    match_z = abs(elem_z - self.aper[jj,2]) < self.dz
                    if all([match_x, match_y, match_z]):
                        if not np.isnan(geos_aper[jj,3]):
                            gmap[jj] = ii
                            isfrac[ii] = True
        return gmap, isfrac


    def adjust_mesh(self, elem, conn, gmap, map, min_ratio=20.0):
        """ Adjust mesh based on aperture

        Change the element volume, interface area and distance to interface
        based on the aperture data. After this adjustment, aperture is
        represented by exactly one element.

        Args:
            elem, conn (list): Obtained from load_mesh()
            gmap (np array)  : Obatined from geos_map()
            map (np array)   : Obtained from neighbor_map()
            min_ratio (double) : Threshold for aperture reduction.

        Returns:
            elem, conn (list): 'elem' and 'conn' with updated values
        """
        N = self.aper.shape[0]
        min_aper = np.amax([np.nanmax(self.aper[:,3])/min_ratio, np.nanmin(self.aper[:,3])])
        print('Minimum allowed aperture, volume = ',min_aper, min_aper*self.dx*self.dz, ' min_ratio = ',min_ratio)
        # loop over all apertures to adjust mesh size
        count = 0
        for ii in range(N):
            if self.aper[ii,3] > 0 and gmap[ii,0] >= 0:
                aperture = self.aper[ii,3]
                if aperture < min_aper:
                    aperture = min_aper
                eid = gmap[ii,0]
                count += 1
                # adjust volume
                if map['jP'][eid,0] >= 0:
                    elem[map['jP'][eid,0]][2] += (elem[eid][2] - aperture*self.dx*self.dz)/2.0
                if map['jM'][eid,0] >= 0:
                    elem[map['jM'][eid,0]][2] += (elem[eid][2] - aperture*self.dx*self.dz)/2.0
                elem[eid][2] = round(aperture * self.dx * self.dz, 4)
                # adjust distance and face area
                if map['jP'][eid,1] >= 0:
                    # distance to face
                    dlP = conn[map['jP'][eid,1]][3] - 0.5 * aperture
                    # conn[map['jP'][eid,1]][4] += dlP
                    conn[map['jP'][eid,1]][4] = (2.0*conn[map['jP'][eid,1]][4] + conn[map['jP'][eid,1]][3] - aperture/2.0)/2.0
                    conn[map['jP'][eid,1]][3] = 0.5 * aperture
                    # area of face
                    if map['kP'][eid,1] >= 0:
                        conn[map['kP'][eid,1]][5] = aperture * self.dx
                    if map['iP'][eid,1] >= 0:
                        conn[map['iP'][eid,1]][5] = aperture * self.dz
                if map['jM'][eid,1] >= 0:
                    dlM = conn[map['jM'][eid,1]][4] - 0.5 * aperture
                    conn[map['jM'][eid,1]][3] = (2.0*conn[map['jM'][eid,1]][3] + conn[map['jM'][eid,1]][4] - aperture/2.0)/2.0
                    conn[map['jM'][eid,1]][4] = 0.5 * aperture
        # enforce minimum area between 2 fracture connections
        for ii in range(N):
            if self.aper[ii,3] > 0 and gmap[ii,0] >= 0:
                aperture = self.aper[ii,3]
                if aperture < min_aper:
                    aperture = min_aper
                eid = gmap[ii,0]
                if map['kM'][eid,1] >= 0:
                    if aperture * self.dx < conn[map['kM'][eid,1]][5]:
                        conn[map['kM'][eid,1]][5] = aperture * self.dx
                if map['iM'][eid,1] >= 0:
                    if aperture * self.dz < conn[map['iM'][eid,1]][5]:
                        conn[map['iM'][eid,1]][5] = aperture * self.dz
            # check results
        print('Max, Min, new Min aperture = ',np.nanmax(self.aper[:,3]),np.nanmin(self.aper[:,3]),min_aper)
        print('Number of elements adjusted = ',count)
        return elem, conn


    def adjust_mesh_poro(self, elem, conn, gmap, map, poro):
        """ Adjust mesh based on aperture

        This is version2 of adjust_mesh. It does not directly change mesh
        volumes. Remeshing is performed by adjusting porosity. But the
        fracture-matrix face area needs to the adjusted accordingly
        """
        self.min_aper = 1e-4
        N = self.aper.shape[0]
        # self.min_aper = np.amax([np.nanmax(self.aper[:,3])*0.8, np.nanmin(self.aper[:,3])])
        # loop over all apertures to adjust mesh size
        count = 0
        for ii in range(N):
            if self.aper[ii,3] > 0 and gmap[ii,0] >= 0 and poro[ii][3] > 0:
                aperture = self.aper[ii,3]
                if aperture < self.min_aper:
                    aperture = self.min_aper
                eid = gmap[ii,0]
                # adjust distance and face area
                if map['jP'][eid,1] >= 0:
                    # fracture-fracture area
                    if map['kP'][eid,1] >= 0:
                        conn[map['kP'][eid,1]][5] = aperture * self.dx
                        count += 1
                    if map['iP'][eid,1] >= 0:
                        conn[map['iP'][eid,1]][5] = aperture * self.dz
                        count += 1
        print('Number of elements adjusted = ',count)
        # enforce minimum area between 2 fracture connections
        count = 0
        for ii in range(N):
            if self.aper[ii,3] > 0 and gmap[ii,0] >= 0:
                aperture = self.aper[ii,3]
                if aperture < self.min_aper:
                    aperture = self.min_aper
                eid = gmap[ii,0]
                if map['kM'][eid,1] >= 0:
                    if aperture * self.dx < conn[map['kM'][eid,1]][5]:
                        conn[map['kM'][eid,1]][5] = aperture * self.dx
                        count += 1
                if map['iM'][eid,1] >= 0:
                    if aperture * self.dz < conn[map['iM'][eid,1]][5]:
                        conn[map['iM'][eid,1]][5] = aperture * self.dz
                        count += 1
        print('Number of elements adjusted = ',count)
        return elem, conn


    def update_mesh(self, elem, conn):
        """ Write the new MESH file

        Args:
            elem, conn (list): Obtained from 'adjust_mesh()'
            isbound (nparray): Obtained from 'markBound()'
        """
        fid = open('MESH_NEW', "w")
        # write elements
        fid.write('ELEME\n')
        for ii in range(self.Nele):
            line = elem[ii]
            for jj in range(2,7):
                line[jj] = self.reformat("{:10.4E}".format(line[jj]))
            fid.write(f"{line[0]}          {line[1]}{line[2]}{line[3]}          {line[4]}{line[5]}{line[6]}")
            if line[7] == 1:
                fid.write(f" I")
            fid.write('\n')
        fid.write('     \n')
        # write connections
        fid.write('CONNE\n')
        for ii in range(self.Ncon):
            line = conn[ii]
            line[2] = str(line[2])
            for jj in range(3,6):
                if isinstance(line[jj],tuple):
                    print(line)
                line[jj] = self.reformat("{:.4E}".format(line[jj]))
            fid.write(f"{line[0]}{line[1]}                   {line[2]}{line[3]}{line[4]}{line[5]}")
            if len(line[6]) > 1:
                fid.write(line[6]+'\n')
            else:
                fid.write(line[6])
        fid.write('     \n')
        fid.close()
        sp.call(['rm', 'MESH'])
        sp.call(['mv', 'MESH_NEW', 'MESH'])
        return 0

"""
    ============================================================================
                                End of Class 'Mesh'
    ============================================================================
"""


"""
    ============================================================================
                    Class 'Geos' for manipulating GEOS output
    ============================================================================
"""
class Geos():
    """ A Geos object represents 1 GEOS output file

    Args:
        fname (str)      : Name (with directory) of the GEOS output file
        options (list)   : Type of operations on raw GEOS output

    """
    def __init__(self, fname, offset, xlim, cutoff):
        self.geos = np.genfromtxt(fname, delimiter=',')
        self.offset = offset
        self.xlim = xlim
        # for ii in range(np.shape(self.geos)[0]):
        #     if self.geos[ii,3] <= cutoff:
        #         self.geos[ii,3] = np.nan

    def convert_to_3D_slice(self, geos, yy):
        xx = np.linspace(3.0, 699.0, 117)
        zz = np.linspace(-2.0, -278.0, 70)
        dim = [len(xx),len(zz)]
        geos3d = np.nan * np.ones((dim))
        for ll in range(np.shape(geos)[0]):
            if yy == geos[ll,1]:
                ii = np.where(xx == geos[ll,0])
                kk = np.where(zz == geos[ll,2])
                geos3d[ii,kk] = geos[ll,3]
        return np.transpose(geos3d), xx, zz

    def convert_to_1D(self, geos3d, xx, zz, yy):
        geos = []
        for ii in range(len(xx)):
            for kk in range(len(zz)):
                if not np.isnan(geos3d[kk,ii]):
                    geos.append([xx[ii], yy, zz[kk], geos3d[kk,ii]])
        return np.array(geos)


    def custom_range(self):
        out = []
        geos = self.geos
        for ii in range(np.shape(geos)[0]):
            loc = [geos[ii,0]+self.offset[0], geos[ii,1]+self.offset[1], geos[ii,2]+self.offset[2]]
            if loc[0]>self.xlim[0] and loc[0]<self.xlim[1]:
                if loc[1]>self.xlim[2] and loc[1]<self.xlim[3]:
                    if loc[2]>self.xlim[4] and loc[2]<self.xlim[5]:
                        if not np.isnan(geos[ii,3]):
                            out.append([loc[0], loc[1], loc[2], geos[ii,3]])
        return np.array(out)


    def block_checking(self, data, blk_key, blk_range, coef):
        """ Reduce aperture for the isolated regions to make model converge

        For the time being, the isolated regions have to be manually delineated
        and stored as blk_range.

        Args:
            blk_range (list)    : [xmin, xmax, zmin, zmax, sweep direction]
        """
        aper = copy.deepcopy(data)


        for n_blk in range(len(blk_range)):
            #   extract aperture within the manually-selected isolated region
            blk_aper = []
            for ii in range(np.shape(aper)[0]):
                if aper[ii,0] > blk_range[n_blk][0] and aper[ii,0] < blk_range[n_blk][1]:
                    if aper[ii,2] > blk_range[n_blk][2] and aper[ii,2] < blk_range[n_blk][3]:
                        blk_aper.append([aper[ii,0],aper[ii,1],aper[ii,2],aper[ii,3]])
            blk_aper = np.array(blk_aper)
            xVec = np.unique(blk_aper[:,0])
            zVec = np.unique(blk_aper[:,2])
            blk_aper2D = np.nan * np.ones((len(zVec), len(xVec)))
            #   map aperture to 2D
            for ii in range(np.shape(blk_aper)[0]):
                indx = abs(xVec - blk_aper[ii,0]).argmin()
                indz = abs(zVec - blk_aper[ii,2]).argmin()
                blk_aper2D[indz, indx] = blk_aper[ii,3]
            #   calculate equivalent aperture for the minimum cross-section
            nn = []
            aa = []
            if blk_range[n_blk][4] == 'x+' or blk_range[n_blk][4] == 'x-':
                for ii in range(len(xVec)):
                    slice = blk_aper2D[:,ii]
                    n_elem = 0
                    a_elem = 0
                    for jj in range(len(zVec)):
                        if not np.isnan(slice[jj]):
                            a_elem += blk_aper2D[jj,ii] ** 3.0
                            n_elem += 1
                    if n_elem > 0:
                        nn.append(ii)
                        aa.append((a_elem/n_elem)**(1.0/3.0))
                aa = np.array(aa)
                nn = np.array(nn)
                #   find the x-value with the minimum equivalent aperture
                min_x = xVec[nn[aa.argmin()]]
                min_aper = np.amin(aa)
                print(' >>> Block checking : minimum eq. aperture = ',min_aper,' at x = ',min_x)
                #   set aperture of the isolated zone to the equivalent aperture
                if blk_range[n_blk][4] == 'x-':
                    for ii in range(np.shape(aper)[0]):
                        if aper[ii,0] < min_x and aper[ii,3] >= min_aper * coef:
                            aper[ii,3] = min_aper * coef
                elif blk_range[n_blk][4] == 'x+':
                    for ii in range(np.shape(aper)[0]):
                        if aper[ii,0] > min_x and aper[ii,3] > min_aper * coef:
                            aper[ii,3] = min_aper * coef
            elif blk_range[n_blk][4] == 'z+' or blk_range[n_blk][4] == 'z-':
                for ii in range(len(zVec)):
                    slice = blk_aper2D[ii,:]
                    n_elem = 0
                    a_elem = 0
                    for jj in range(len(xVec)):
                        if not np.isnan(slice[jj]):
                            a_elem += blk_aper2D[ii,jj] ** 3.0
                            n_elem += 1
                    if n_elem > 0:
                        nn.append(ii)
                        aa.append((a_elem/n_elem)**(1.0/3.0))
                aa = np.array(aa)
                nn = np.array(nn)
                #   find the x-value with the minimum equivalent aperture
                min_z = zVec[nn[aa.argmin()]]
                min_aper = np.amin(aa)
                print(' >>> Block checking : minimum eq. aperture = ',min_aper,' at z = ',min_z)
                #   set aperture of the isolated zone to the equivalent aperture
                if blk_range[n_blk][4] == 'z-':
                    for ii in range(np.shape(aper)[0]):
                        if aper[ii,2] < min_z and aper[ii,3] > min_aper * coef:
                            aper[ii,3] = min_aper * coef
                elif blk_range[n_blk][4] == 'z+':
                    coef = 2.5
                    for ii in range(np.shape(aper)[0]):
                        if aper[ii,2] > min_z and aper[ii,3] > min_aper * coef:
                            aper[ii,3] = min_aper * coef
            else:
                raise ValueError("Currently only supports x-sweep!")

        # remove large aperture due to proppants
        geos3d, xx, zz = self.convert_to_3D_slice(aper, 8.0)
        if blk_key == 'upshi3':
            for kk in range(50,80):
                geos3d[63,kk] = 0.66*geos3d[62,kk] + 0.34*geos3d[65,kk]
                geos3d[64,kk] = 0.34*geos3d[62,kk] + 0.66*geos3d[65,kk]
            # for ii in range(40,45):
            #     for kk in range(50,65):
            #         if geos3d[ii,kk] > 0:
            #             geos3d[ii,kk] = geos3d[ii,kk] + 1e-5
            # for kk in range(50,65):
            #     if geos3d[45,kk] > 3e-4:
            #         geos3d[45,kk] = 3e-4
        elif blk_key == 'upshi2':
            for kk in range(40,70):
                geos3d[62,kk] = 0.75*geos3d[61,kk] + 0.25*geos3d[65,kk]
                geos3d[63,kk] = 0.5*geos3d[61,kk] + 0.5*geos3d[65,kk]
                geos3d[64,kk] = 0.25*geos3d[61,kk] + 0.75*geos3d[65,kk]
            for ii in range(40,45):
                for kk in range(50,65):
                    if geos3d[ii,kk] > 0:
                        geos3d[ii,kk] = geos3d[ii,kk] + 1e-5
        elif blk_key == 'upshi1':
            for kk in range(45,75):
                geos3d[61,kk] = 0.75*geos3d[60,kk] + 0.25*geos3d[64,kk]
                geos3d[62,kk] = 0.5*geos3d[60,kk] + 0.5*geos3d[64,kk]
                geos3d[63,kk] = 0.25*geos3d[60,kk] + 0.75*geos3d[64,kk]
            for ii in range(40,45):
                for kk in range(55,65):
                    if geos3d[ii,kk] > 0:
                        geos3d[ii,kk] = geos3d[ii,kk] + 1e-5
        aper = self.convert_to_1D(geos3d, xx, zz, 8.0)

        return aper



    def scale(self, data, scale, min_aper=1e-6, max_aper=5e-2):
        """ Scale geos data to match target mean/std"""
        targ_avg = scale[0]
        targ_std = scale[1]
        out = copy.deepcopy(data)

        #   Adjust aperture to match desired mean and std
        if targ_avg != None:
            actv = out[:,3] > 0.0
            new_mean = np.mean(out[actv,3])
            current_std = np.std(out[actv,3])
            iter = 0
            while abs(new_mean - targ_avg)/targ_avg > 0.05:
                for ii in range(np.shape(out)[0]):
                    if actv[ii] == 1:
                        out[ii,3] = data[ii,3] - (np.mean(data[actv,3]) - targ_avg)
                new_mean = np.mean(out[actv,3])
                if targ_std != None:
                    coeff = targ_std / np.std(data[actv,3])
                    out[actv,3] = (out[actv,3] - new_mean) * coeff + new_mean
                # Iterate to remove small aperture
                for ii in range(np.shape(out)[0]):
                    if actv[ii] == 1 and out[ii,3] <= min_aper:
                        out[ii,3] = min_aper
                new_mean = np.mean(out[actv,3])
                print('Scaling GEOS: iter ',iter,' aperture has mean, std = ',np.mean(out[actv,3]),np.std(out[actv,3]))
                iter += 1
                if iter > 10:
                    break
        elif targ_std != None:
            actv = out[:,3] > 0.0
            new_mean = np.mean(out[actv,3])
            current_std = np.std(out[actv,3])
            coeff = targ_std / np.std(data[actv,3])
            out[actv,3] = (out[actv,3] - new_mean) * coeff + new_mean

            for ii in range(np.shape(out)[0]):
                if actv[ii] == 1 and out[ii,3] <= 0.0:
                    diff = abs(min_aper - out[ii,3])
                    out[ii,3] = min_aper
                    aa = out[:,3] > new_mean
                    incre = diff/sum(aa)
                    print(ii,' of ',np.shape(out)[0],': ',incre)
                    out[aa,3] -= incre
            new_mean = np.mean(out[actv,3])
            print('Scaling GEOS: aperture has mean, std = ',np.mean(out[actv,3]),np.std(out[actv,3]))
        else:
            actv = out[:,3] > 0.0
        #   remove large and small aperture
        # if targ_std == None or targ_std < current_std:
        #     aa = out[:,3] < self.min_aper
        #     out[aa,3] = np.nan
        # aa = out[:,3] > max_aper
        # out[aa,3] = max_aper
        # aa = out[:,3] < min_aper
        # out[aa,3] = min_aper
        for av in range(np.shape(out)[0]):
            if actv[av] == 0:
                out[av,3] = np.nan
        print('Final GEOS: aperture has mean, std = ',np.nanmean(out[actv,3]),np.nanstd(out[actv,3]))
        return out
"""
    ============================================================================
                                End of Class 'Geos'
    ============================================================================
"""



"""
    ============================================================================
                    Class 'Incon' for generating TOUGH INCON
    ============================================================================
"""
class Incon():
    """ INCON for 1 variable

    Args:
        Nele (int)      : Number of elements in the TOUGH MESH file
        const (float)   : Default initial values for the fracture

    Attributes:
        domain (array)  : Nele by 1 array initialized with const
        Nele (int)      : Number of elements in TOUGH MESH

    """
    def __init__(self, const, Nele, import_from_geos):
        self.Nele = Nele
        self.domain = const * np.ones((Nele,), dtype=float)
        self.import_from_geos = import_from_geos


    def insert_geos(self, data, gmap):
        """ Insert GEOS values into TOUGH mesh

        Args:
            data (nparray)  : GEOS results, which is a N by 4 array
            gmap (nparray)  : gmap[row index in 'geos'] = row index in TOUGH MESH

        Returns:
            out1d (nparray) : Initial values for all elements. Values in fractures
                              are read from GEOS. Values in matrix are set to default.
        """
        out1d = self.domain
        if self.import_from_geos:
            N = data.shape[0]
            for ii in range(0,N):
                if not np.isnan(gmap[ii,0]):
                    if not np.isnan(data[ii,3]):
                        out1d[gmap[ii,0]] = data[ii,3]
            print('INCON: (mean,std,max,min = ', np.nanmean(out1d), \
                    np.nanstd(out1d),np.nanmax(out1d),np.nanmin(out1d),')')
        return out1d

"""
    ============================================================================
                                End of Class 'Incon'
    ============================================================================
"""



"""
    ============================================================================
                    Class 'Input' for writing INPUT file
    ============================================================================
"""
class Input():
    def __init__(self, elem, conn, use_minc):
        self.elem = elem
        self.conn = conn
        self.Nele = len(elem)
        self.Ncon = len(conn)
        self.use_minc = use_minc

    def write_incon(self, data, state, zone, isfrac, perm0):
        """ Write the INCON file for inserted geos data """
        fid = open('INCON', 'w')
        fid.write('INCON\n')
        num_blocks = 0
        for ii in range(0,len(self.elem)):
            if self.elem[ii][1] in zone and self.elem[ii][7] == 0 and isfrac[ii] == 1:
                eid = self.elem[ii][0]
                poro = "{:.8E}".format(data[ii,0])
                perm = "{:.8E}".format(data[ii,1])
                pres = "{:.3E}".format(data[ii,2])
                gaor = "{:.3E}".format(data[ii,3])
                sato = "{:.3E}".format(data[ii,4])
                temp = "{:.3E}".format(data[ii,5])
                perm0str = "{:.8E}".format(perm0)
                # line1: element, porosity(3), state(4), permeability(5-7)
                fid.write(f"{eid}           {poro}  {state}                                     {perm} {perm0str} {perm}")
                fid.write('\n')
                # line2: primary variables (2,4,6)
                fid.write(f"           {pres}           {gaor}           {sato}           {temp}")
                fid.write('\n')
                num_blocks += 1
        fid.write("<<<")
        fid.close()
        print("Total number of grid blocks in INCON = ",num_blocks)

    def write_memory(self, finput):
        fin = open(finput, 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'Cartesian' in line:
                fout.write(f"'Cartesian'   {self.Nele}   {self.Ncon}   5   .FALSE.   .FALSE.\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()

    def write_rocks(self, region, boundary):
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'ROCKS' in line:
                fout.write(line)
                i_rock = 1
                for i_reg in range(1, len(region)):
                    for i_subreg in range(1, len(region[i_reg])):
                        reg_name = region[i_reg][0]
                        # if use minc, write minc properties first
                        if self.use_minc:
                            if len(region[i_reg][i_subreg]) <= 7:
                                i_prop = 6
                            else:
                                i_prop = 7
                            param = region[i_reg][i_subreg][i_prop]
                            reg_name2 = reg_name[1:]+'F'
                            reg_name = reg_name[1:]+'M'
                            rho = "{:4.1E}".format(param[1])
                            phi = "{:4.1E}".format(param[2])
                            perm = "{:8.2E}".format(param[3])
                            comp = "{:5.1E}".format(param[4])
                            ind_rp = param[5]
                            ind_cp = param[6]
                            # get rp and cp coefficients
                            rpCoeff = self.relPermParam(ind_rp)
                            cpCoeff = self.capPresParam(ind_cp)
                            fout.write(f"{reg_name2}    2   {rho}   {phi}  {perm}  {perm}  {perm}       4.0    1000.0\n")
                            fout.write(f"   {comp}               1.5e0\n")
                            fout.write(f"    {ind_rp}      ")
                            for ii in range(len(rpCoeff)):
                                coef = "{:7.3E}".format(rpCoeff[ii])
                                fout.write(f"{coef} ")
                            fout.write("\n")
                            fout.write(f"    {ind_cp}      ")
                            for ii in range(len(cpCoeff)):
                                coef = "{:7.3E}".format(cpCoeff[ii])
                                fout.write(f"{coef} ")
                            fout.write("\n")
                        # write regular properties
                        param = region[i_reg][i_subreg][6]
                        rho = "{:4.1E}".format(param[1])
                        phi = "{:4.1E}".format(param[2])
                        perm = "{:8.2E}".format(param[3])
                        comp = "{:5.1E}".format(param[4])
                        ind_rp = param[5]
                        ind_cp = param[6]
                        # get rp and cp coefficients
                        rpCoeff = self.relPermParam(ind_rp)
                        cpCoeff = self.capPresParam(ind_cp)
                        fout.write(f"{reg_name}    2   {rho}   {phi}  {perm}  {perm}  {perm}       4.0    1000.0\n")
                        fout.write(f"   {comp}               1.5e0\n")
                        fout.write(f"    {ind_rp}      ")
                        for ii in range(len(rpCoeff)):
                            coef = "{:7.3E}".format(rpCoeff[ii])
                            fout.write(f"{coef} ")
                        fout.write("\n")
                        fout.write(f"    {ind_cp}      ")
                        for ii in range(len(cpCoeff)):
                            coef = "{:7.3E}".format(cpCoeff[ii])
                            fout.write(f"{coef} ")
                        fout.write("\n")
                # write the well
                for i_reg in range(1, len(boundary)):
                    reg_name = boundary[i_reg][0]
                    param = boundary[i_reg][1][6]
                    rho = "{:4.1E}".format(param[1])
                    phi = "{:4.1E}".format(param[2])
                    perm = "{:8.2E}".format(param[3])
                    comp = "{:5.1E}".format(param[4])
                    ind_rp = param[5]
                    ind_cp = param[6]
                    # get rp and cp coefficients
                    rpCoeff = self.relPermParam(ind_rp)
                    cpCoeff = self.capPresParam(ind_cp)
                    fout.write(f"{reg_name}    2   {rho}   {phi}  {perm}  {perm}  {perm}       4.0    1000.0\n")
                    fout.write(f"   {comp}               1.5e0\n")
                    fout.write(f"    {ind_rp}      ")
                    for ii in range(len(rpCoeff)):
                        coef = "{:7.3E}".format(rpCoeff[ii])
                        fout.write(f"{coef} ")
                    fout.write("\n")
                    fout.write(f"    {ind_cp}      ")
                    for ii in range(len(cpCoeff)):
                        coef = "{:7.3E}".format(cpCoeff[ii])
                        fout.write(f"{coef} ")
                    fout.write("\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()
        sp.call(['rm', 'INPUT0'])

    def write_param(self, time, ic):
        Tend = "{:9.4E}".format(time[0])
        dt = "{:7.2E}".format(time[1])
        dt_max = "{:7.2E}".format(time[2])
        state = ic['Default'][0]
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'PARAM' in line:
                fout.write(line)
                fout.write(f"   3-999    -999100030010000000400803010   0.00E-5                        100\n")
                fout.write(f"          {Tend}  {dt}  {dt_max}              9.8060     1.5e0\n")
                fout.write(f"     1.E-5     1.E01                                  1.0e-8            {state}\n")
                for ii in range(1,len(ic['Default'])):
                    val = "{:10.6E}".format(ic['Default'][ii])
                    fout.write(f"        {val}")
                fout.write(f"\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()
        sp.call(['rm', 'INPUT0'])

    def write_indom(self, ic, use_minc):
        ickeys = ic.keys()
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'INDOM' in line:
                fout.write(line)
                for key in ickeys:
                    if key != 'Default':
                        elem = ic[key]
                        reg_name = key
                        if use_minc:
                            if 'WELL' in reg_name:
                                state = elem[1]
                                fout.write(f"{reg_name}  {state}\n")
                                for ii in range(2,len(elem)):
                                    val = "{:10.6E}".format(elem[ii])
                                    fout.write(f"        {val}")
                                fout.write(f"\n")
                            else:
                                reg_name = reg_name[1:] + 'M'
                                state = elem[1]
                                fout.write(f"{reg_name}  {state}\n")
                                for ii in range(2,len(elem)):
                                    val = "{:10.6E}".format(elem[ii])
                                    fout.write(f"        {val}")
                                fout.write(f"\n")
                                reg_name = reg_name[:-1] + 'F'
                                fout.write(f"{reg_name}  {state}\n")
                                for ii in range(2,len(elem)):
                                    val = "{:10.6E}".format(elem[ii])
                                    fout.write(f"        {val}")
                                fout.write(f"\n")
                        else:
                            state = elem[1]
                            fout.write(f"{reg_name}  {state}\n")
                            for ii in range(2,len(elem)):
                                val = "{:10.6E}".format(elem[ii])
                                fout.write(f"        {val}")
                            fout.write(f"\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()
        sp.call(['rm', 'INPUT0'])

    def write_times(self, tVec):
        Nt = len(tVec)
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'TIMES' in line:
                fout.write(line)
                fout.write(f"{str(Nt)}\n")
                col = 0
                for ii in range(Nt):
                    val = "{:7.3E}".format(tVec[ii])
                    fout.write(f" {val}")
                    col += 1
                    if col == 8 and ii != Nt-1:
                        fout.write(f"\n")
                fout.write(f"\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()
        sp.call(['rm', 'INPUT0'])

    def write_conx(self, conx):
        N = len(conx)
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'Conx_Time_Series' in line:
                fout.write(line)
                for ii in range(N):
                    val = conx[ii]
                    fout.write(f"{val}\n")
                fout.write(f"\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()
        sp.call(['rm', 'INPUT0'])

    def relPermParam(self, num):
        if num == 8:
            param = [0.45, 0.15, 0.01, 2.0, 2.0]
            # param = [0.4, 0.21, 0.01, 2.0, 2.0]
        elif num == 6:
            param = [0.05,0.05,0.0,2.0]
        else:
            raise ValueError('Relative permeability equation number is not available!!!')
        return param

    def capPresParam(self, num):
        if num == 5:
            param = []
        elif num == 8:
            # param = [0.4, 1.85, 10.0, 11.0, 0.21]
            param = [0.45, 1.85, 10.0, 11.0, 0.15]
        else:
            raise ValueError('Capillary pressure equation number is not available!!!')
        return param
"""
    ============================================================================
                                End of Class 'Input'
    ============================================================================
"""
