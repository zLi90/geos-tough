""" A class that read an existing TOUGH mesh, adjust its size to match the
    fracture data (from GEOS), and add boundary tags (I).
    For the time being, fractures are assumed 2D in x-z plane.

    TODO: Validate 'adjustMesh' using real aperture data from GEOS
    TODO: Handle fractures in arbitrary directions
    TODO: Improve search efficiency in 'buildMap'
    TODO: Combine with tough-convert for loading MESH parameters
    TODO: Enable mesh-building when TOUGH and GEOS mesh resolutions are different
"""
import numpy as np
import subprocess as sp
import copy
import os
import math

class Mesh():
    """ The class for building TOUGH mesh with fractures

    Args:
        fmesh (str)     : Name of an existing mesh file, usually 'MESH'
        fmt (str)       : MeshMaker input (used to generate 'MESH' if not provided)
        delta (list)    : Increments in each direction (in meters), [dx, dy, dz]
        dim (list)      : Number of elements in each direction, [Nx, Ny, Nz]
        bound (list)    : List of domain boundaries with fixed values.
        wname (str)     : Media name of the well (assume only 1 well exists for now)

    Attributes:
        dx/dy/dz (float): Element size in each direction
        Nx/Ny/Nz (int)  : Number of elements in each direction
        fmesh (str)     : Name of existing TOUGH mesh file
        fmt (str)       : MeshMaker input
        well (str)      : Media name of the well
        bound (list)    : List of boundaries with fixed boundary values
        Nele (int)      : Number of elements
        Ncon (int)      : Number of connections
        xVec (nparray)  : Coordinate values in x direction
        yVec (nparray)  : Coordinate values in y direction
        zVec (nparray)  : Coordinate values in z direction
        aper (nparray)  : GEOS aperture data read from .csv file

    """
    def __init__(self, fmesh, rocks, delta, dim, bound, bc, wname=['WELL1']):
        self.dx = delta[0]
        self.dy = delta[1]
        self.dz = delta[2]
        self.delta = delta
        self.Nx = dim[0]
        self.Ny = dim[1]
        self.Nz = dim[2]
        self.well = wname
        self.bound = bound
        self.useBC = bc['useBC']
        self.rocks = rocks
        if type(fmesh) is str:
            self.fmesh = fmesh
            self.fmt = None
        else:
            self.fmesh = 'MESH'
            self.fmt = fmesh

    def genMeshInput(self):
        """ Generate input file for MeshMaker

        Returns:
            finput (str): Name of the MeshMaker input file created
        """
        fmt = self.fmt
        Nreg = 1
        for i_reg in range(1,len(fmt['Regions'])):
            reg_range = self.rocks[fmt['Regions'][i_reg]][0]
            Nreg += (len(reg_range) - 1)
        finput = 'meshinput'
        fid = open(finput, 'w')
        fid.write(f"Cartesian grid\n{str(fmt['MaxElem'])} {str(fmt['Longest'])} 5 'Old' 'm' \n")
        fid.write(f"Regions\n{str(Nreg)}\n'{fmt['Regions'][0]}'\n")
        if len(fmt['Regions']) > 1:
            for i_reg in range(1,len(fmt['Regions'])):
                reg_name = fmt['Regions'][i_reg]
                reg_range = self.rocks[reg_name][0]
                if len(reg_range) > 1:
                    for i_subreg in range(1,len(reg_range)):
                        fid.write(f"'{reg_name}'\n'Cartesian'  'm' \n")
                        for jj in range(0,6):
                            fid.write(str(reg_range[i_subreg][jj])+" ")
                        fid.write('\n')
        fid.write(f"XYZ\n       00.    0.0d0    0.0d0    0.0d0\n")
        self.writeIncre(fid, fmt, 'NX')
        self.writeIncre(fid, fmt, 'NY')
        self.writeIncre(fid, fmt, 'NZ')
        fid.write('<<<')
        fid.close()
        return finput

    def writeIncre(self, fid, fmt, dir):
        """ Write element increment section of MeshMaker input

        Args:
            fid         : File handle
            fmt (str)   : MashMaker input
            dir (str)   : Direction of increment ('NX', 'NY' or 'NZ')
        """
        N = "{:5d}".format(fmt[dir][0])
        fid.write(f"{dir}   {N}")
        if type(fmt[dir][1]) is list:
            fid.write('\n')
            count = 0
            for ii in range(len(fmt[dir][1])):
                # if constant grid spacing
                if len(fmt[dir][1][ii]) == 2:
                    for jj in range(fmt[dir][1][ii][0]):
                        fid.write("{:8.4E}".format(fmt[dir][1][ii][1]))
                        count += 1
                        if count % 8 == 0 and count != fmt[dir][0]:
                        # if count % 8 == 0 and ii != len(fmt[dir][1])-1:
                            fid.write('\n')
                # if geometric grid spacing
                elif len(fmt[dir][1][ii]) == 3:
                    a1 = fmt[dir][1][ii][1]
                    a_sum = fmt[dir][1][ii][2]
                    ncell = fmt[dir][1][ii][0]
                    if ncell % 2 == 0:
                        r = self.commonR(a1, a_sum/2, ncell/2)
                        for jj in range(int(ncell/2)):
                            print(fmt[dir][1][ii][1] * r**jj)
                            fid.write("{:8.4E}".format(fmt[dir][1][ii][1] * r**jj))
                            count += 1
                            if count % 8 == 0 and count != fmt[dir][0]:
                                fid.write('\n')
                        for jj in range(int(ncell/2)):
                            print(fmt[dir][1][ii][1] * r**(ncell/2-1) / r**(jj))
                            fid.write("{:8.4E}".format(fmt[dir][1][ii][1] * r**(ncell/2-1) / r**(jj)))
                            count += 1
                            if count % 8 == 0 and count != fmt[dir][0]:
                                fid.write('\n')
                    else:
                        r = self.commonR(a1, a_sum, ncell)
                        for jj in range(ncell):
                            print(fmt[dir][1][ii][1] * r**jj)
                            fid.write("{:8.4E}".format(fmt[dir][1][ii][1] * r**jj))
                            count += 1
                            if count % 8 == 0 and count != fmt[dir][0]:
                                fid.write('\n')
                # if log grid spacing
                elif len(fmt[dir][1][ii]) == 4:
                    a1 = fmt[dir][1][ii][1]
                    a_sum = fmt[dir][1][ii][2]
                    ncell = fmt[dir][1][ii][0]
                    sign = fmt[dir][1][ii][3]
                    x1 = math.exp(a1)
                    if ncell % 2 == 0:
                        dx = self.logRoot(x1, a_sum/2.0, int(ncell/2))
                        x2 = x1 + (ncell/2)*dx
                        area = 0
                        lst = []
                        for jj in range(int(ncell/2)):
                            lst.append(math.log(x1 + jj*dx))
                            fid.write("{:8.4E}".format(math.log(x1 + jj*dx)))
                            count += 1
                            if count % 8 == 0 and count != fmt[dir][0]:
                                fid.write('\n')
                        for jj in range(int(ncell/2)):
                            fid.write("{:8.4E}".format(lst[int(ncell/2)-jj-1]))
                            count += 1
                            if count % 8 == 0 and count != fmt[dir][0]:
                                fid.write('\n')
                    else:
                        dx = self.logRoot(x1, a_sum, ncell)
                        if sign == '+':
                            for jj in range(ncell):
                                fid.write("{:8.4E}".format(math.log(x1 + jj*dx)))
                                count += 1
                                if count % 8 == 0 and count != fmt[dir][0]:
                                    fid.write('\n')
                        elif sign == '-':
                            xn = x1 + ncell*dx
                            for jj in range(ncell):
                                fid.write("{:8.4E}".format(math.log(xn - (jj+1)*dx)))
                                count += 1
                                if count % 8 == 0 and count != fmt[dir][0]:
                                    fid.write('\n')
                else:
                    raise ValueError("Length of a mesh segment should be 2, 3 or 4!")

        else:
            fid.write("{:8.4E}".format(fmt[dir][1]))
        fid.write('\n')

    def logRoot(self, x1, L, n, eps=1, tolerance=1e-8):
        """ Find end point of the log function """
        x_left = (x1-1)
        x_right = 1.0
        f_left = self.logSum(x1, x_left, n) - L
        f_right = self.logSum(x1, x_right, n) - L
        while eps > tolerance:
            if f_left * f_right < 0:
                x_right -= 0.5 * (x_right - x_left)
            else:
                x_right += (x_right - x_left)
                x_left += 0.5 * (x_right - x_left)
            f_left = self.logSum(x1, x_left, n) - L
            f_right = self.logSum(x1, x_right, n) - L
            eps = x_right - x_left
        return 0.5 * (x_right+x_left)

    def logSum(self, x1, dx, n):
        xsum = 0
        for ii in range(n):
            xsum += math.log(x1+ii*dx)
        return xsum

    def commonR(self, a1, a_sum, n, eps=1, tolerance=1e-8):
        """ Find common ratio of a grometric series using bisection """
        x_left = 1
        x_right = a_sum/a1
        f_left = self.geometricSum(a1, n, x_left) - a_sum
        f_right = self.geometricSum(a1, n, x_right) - a_sum
        while eps > tolerance:
            if f_left * f_right < 0:
                x_right -= 0.5 * (x_right - x_left)
            else:
                x_right += (x_right - x_left)
                x_left += 0.5 * (x_right - x_left)
            f_left = self.geometricSum(a1, n, x_left) - a_sum
            f_right = self.geometricSum(a1, n, x_right) - a_sum
            eps = x_right - x_left
        return 0.5 * (x_right+x_left)

    def geometricSum(self, a1, n, r):
        """ Returns the sum of a geometric series """
        if r != 1:
            return a1 * (1 - r**n) / (1 - r)
        else:
            return a1 * n

    def genMesh(self, fname):
        """ Calls MeshMake to generate 'MESH'

        Args:
            fname (str) : Name of MeshMaker input file created by 'genMeshInput'
        """
        with open(fname,'r') as fid:
            fmtstr = fid.read()
        proc = sp.Popen(['MeshMaker'], stdin=sp.PIPE)
        proc.communicate(str.encode(fmtstr))
        return 0

    def readMesh(self):
        """ Read an existing MESH file

        Returns:
            elem (list): List of all elements in ELEME block
                         Each element is a list with format:
                         [elem id, media, volume, area, x, y, z]
            conn (list): List of all connections in CONNE block
                         Each connection is a list with format:
                         [conn id1, id2, axis(x,y,z), l1, l2, area]
        """
        fm = open(self.fmesh, 'r')
        elem = []
        conn = []
        zone = []
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
                # elem has the structure (element id, media, V, A, x, y, z)
                elem.append([line[0:5], line[15:20], float(line[20:30]), float(line[30:40]),
                    float(line[50:60]), float(line[60:70]), float(line[70:80])])
        fm.close()
        self.Nele = len(elem)
        self.Ncon = len(conn)

        # read (x,y,z) from elem
        self.xVec = np.sort(np.unique(np.array([item[4] for item in elem])))
        self.yVec = np.sort(np.unique(np.array([item[5] for item in elem])))
        self.zVec = np.sort(np.unique(np.array([item[6] for item in elem])))
        self.zVec[::-1].sort()
        return elem, conn


    def coordNearFrac(self, elem, yfrac):
        """ Identify repeated coordinates near fracture
            Note that these coordinates do not affect mesh size, but they
            are used as reference when inserting GEOS values.

            Assumptions:
                - Fracture in x-z plane
                - Only 1 layer of repeated coordinates on each side of fracture

            Returns:
                elem    : elem with updated coordinates
                remesh  : Boolean to indicate if remesh is needed
        """

        elem_y = []
        elem_ind = []
        for ll in range(len(elem)):
            if elem[ll][5] in yfrac:
                elem_y.append(elem[ll][5])
                elem_ind.append(ll)

        elem_y_np = np.array(elem_y)

        remesh = 0
        for ff in range(len(yfrac)):
            yy = yfrac[ff]
            aa = elem_y_np == yy
            # if NY dimension is the largest
            if self.Ny > self.Nx and self.Ny > self.Nz:
                if np.nansum(aa) == 3 * self.Nx*self.Nz:
                    count = 0
                    remesh = 1
                    for ll in range(len(elem)):
                        if elem[ll][5] == yy:
                            count += 1
                            if count <= self.Nx*self.Nz:
                                elem[ll][5] -= 0.01
                            elif count > 2*self.Nx*self.Nz:
                                elem[ll][5] += 0.01
                elif np.nansum(aa) == 2 * self.Nx*self.Nz:
                    count = 0
                    remesh = 1
                    for ll in range(len(elem)):
                        if elem[ll][5] == yy:
                            count += 1
                            if count <= self.Nx*self.Nz:
                                elem[ll][5] -= 0.01
            # if NY dimension is the smallest
            elif self.Ny < self.Nx and self.Ny < self.Nz:
                i_start = 0
                while i_start < len(elem)-self.Ny:
                    count = 0
                    for ii in range(self.Ny-1):
                        if elem[i_start+ii][5] == yy and count == 0:
                            elem[i_start+ii][5] -= 0.01
                            count += 1
                        elif elem[i_start+ii][5] == yy and count == 1:
                            remesh = 1
                            count += 1
                        elif elem[i_start+ii][5] == yy and count == 2:
                            elem[i_start+ii][5] += 0.01
                            break;
                    i_start += self.Ny
            # for 2D case
            elif self.Nz == 1 and self.Ny < self.Nx:
                i_start = 0
                while i_start < len(elem)-self.Ny:
                    count = 0
                    for ii in range(self.Ny-1):
                        if elem[i_start+ii][5] == yy and count == 0:
                            elem[i_start+ii][5] -= 0.01
                            count += 1
                        elif elem[i_start+ii][5] == yy and count == 1:
                            remesh = 1
                            count += 1
                        elif elem[i_start+ii][5] == yy and count == 2:
                            elem[i_start+ii][5] += 0.01
                            break;
                    i_start += self.Ny
            else:
                raise ValueError("Dimension NY is in between NX and NZ!")
        return elem, remesh



    def buildMap(self, elem, conn, geos, zone):
        """ Build map between GEOS and TOUGH mesh

        Args:
            elem, conn (list): Obtained from readMesh()
            geos (nparray)   : Aperture data from GEOS
                               Has size (N by 4) with rows [x, y, z, aperture]
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
        print(self.xVec)
        print(self.yVec)
        print(self.zVec)
        for ii in range(dim[0]):
            # print(np.shape(self.aper))
            rowx = self.findInd(self.aper[ii,0], self.xVec, self.dx)
            rowy = self.findInd(self.aper[ii,1], self.yVec, self.dy)
            rowz = self.findInd(self.aper[ii,2], self.zVec, self.dz)
            # Note that the order of elements in TOUGH seems depend on the
            # dimension of the domain. So here we sort Nx, Ny and Nz to get
            # the correct element index.
            # if any([np.isnan(rowx[0][0]), np.isnan(rowy[0][0]), np.isnan(rowz[0][0])]):
            if any([isinstance(rowx,float),isinstance(rowy,float),isinstance(rowz,float)]):
                gmap[ii] = -1
            else:
                rowarray = np.array([[self.Nz, rowz[0]],[self.Ny, rowy[0]],[self.Nx, rowx[0]]])
                rows = rowarray[rowarray[:,0].argsort()]
                row = rows[2,1] * rows[0,0] * rows[1,0] + rows[1,1] * rows[0,0] + rows[0,1]
                # print(row)
                if len(row) > 1:
                    print(self.aper[ii,:])
                    print(rowx,rowy,rowz)
                row = int(row)
                if elem[row][1] in zone:
                    if np.isnan(self.aper[ii,3]):
                        gmap[ii] = -1
                        isfrac[row] = False
                    else:
                        gmap[ii] = row
                        isfrac[row] = True
                else:
                    gmap[ii] = -1
                    isfrac[row] = False
        return gmap, isfrac

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

    def getNeighbor(self, elem, conn):
        """ Build map among TOUGH elements

        Args:
            elem, conn (list): Obtained from readMesh()

        Returns:
            map (dict): Dict with keys 'iP'(x-plus), 'iM'(x-minus), 'jP', 'jM', 'kP', 'kM'
                        e.g., map['iP'] is an (N by 2) np array with rows:
                        [ind_elem, ind_conn]
                        ind_elem: index of x-plus element of an element
                        ind_conn: index of connection between an element and its x-plus element
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


    def adjustMesh(self, elem, conn, gmap, map, min_ratio=20.0):
        """ Adjust mesh based on aperture

        Change the element volume, interface area and distance to interface
        based on the aperture data. After this adjustment, aperture is
        represented by exactly one element.

        Args:
            elem, conn (list): Obtained from readMesh()
            gmap (np array)  : Obatined from buildMap()
            map (np array)   : Obtained from getNeighbor()

        Returns:
            elem, conn (list): 'elem' and 'conn' with updated values
        """
        N = self.aper.shape[0]
        # set the minimum aperture
        min_aper = np.nanmax(self.aper[:,3]) / min_ratio
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
                        # print('z-area corrected from ',conn[map['kM'][eid,1]][5],' to ',self.aper[ii,3] * self.dx)
                        conn[map['kM'][eid,1]][5] = aperture * self.dx
                if map['iM'][eid,1] >= 0:
                    if aperture * self.dz < conn[map['iM'][eid,1]][5]:
                        # print('x-area corrected from ',conn[map['iM'][eid,1]][5],' to ',self.aper[ii,3] * self.dz)
                        conn[map['iM'][eid,1]][5] = aperture * self.dz

            # check results
        print('Max, Min, new Min aperture = ',np.nanmax(self.aper[:,3]),np.nanmin(self.aper[:,3]),min_aper)
        print('Number of elements adjusted = ',count)
        # for ii in range(len(conn)):
        #     dist1 = conn[ii][3]
        #     dist2 = conn[ii][4]
        #     area = float(conn[ii][5])
        #
        #     if dist1 <= 1e-5 or dist1 >= 0.0125:
        #         print('dist1 = ',dist1,' conn = ',conn[ii][0],conn[ii][1])
        #     if dist2 <= 1e-5 or dist2 >= 0.0125:
        #         print('dist2 = ',dist2,' conn = ',conn[ii][0],conn[ii][1])
        #     if area <= 1e-4:
        #         print('area = ',area, ' conn = ',conn[ii][0],conn[ii][1])

        return elem, conn

    def markBound(self, elem):
        ''' Identify boundary elements

        Add label 'I' to the end of elements to indicate boundary condition
        Add label 'Vxx' to indicate time-variable boundary conditions

        Args:
            elem (list): Obtained from 'adjustMesh'

        Returns:
            isbound (nparray): Equals 1 if an element is tagged with 'I', 0 otherwise
            iswell (nparray) : Equals 1 if an element belongs to a well, 0 otherwise
        '''
        isbound = np.zeros((self.Nele,1), dtype=int)
        iswell = np.zeros((self.Nele,1), dtype=int)
        xrange = [np.amin(self.xVec), np.amax(self.xVec)]
        yrange = [np.amin(self.yVec), np.amax(self.yVec)]
        zrange = [np.amin(self.zVec), np.amax(self.zVec)]
        for ii in range(self.Nele):
            media = elem[ii][1]
            # mark domain boundaries
            if 'x+' in self.bound and elem[ii][4] >= xrange[1]:
                isbound[ii] = 1
            if 'x-' in self.bound and elem[ii][4] <= xrange[0]:
                isbound[ii] = 1
            if 'y+' in self.bound and elem[ii][5] >= yrange[1]:
                isbound[ii] = 1
            if 'y-' in self.bound and elem[ii][5] <= yrange[0]:
                isbound[ii] = 1
            if 'z+' in self.bound and elem[ii][6] >= zrange[1]:
                isbound[ii] = 1
            if 'z-' in self.bound and elem[ii][6] <= zrange[0]:
                isbound[ii] = 1
            # well from zones
            if media in self.well:
                iswell[ii] = 1
                ind = self.well.index(media)
                isbound[ii] = ind+2
        return isbound, iswell

    def getWellInterface(self, elem, conn, map, iswell):
        """ Print connection index for well connections

        This function will print connection indexes if they connect to the well.
        These indexes can be copied into the TOUGH input file for calculating
        well flow rate using the INTERFACE block.

        Args:
            elem (list) : Obtained from 'adjustMesh'
            conn (list) : Obtained from 'adjustMesh'

        Returns:
            allid (dict): All CONNE id for each well
        """
        count = 0
        allid = {}
        dir = ['iP','iM','kP','kM']
        for ii in range(1,len(self.fmt['Regions'])):
            allid[self.fmt['Regions'][ii]] = []
        for ii in range(self.Nele):
            if iswell[ii][0] == 1:
                wellid = elem[ii][0]
                for jj in range(self.Ncon):
                    # Only consider well-matrix connections (not well-well connections)
                    # Also, this assumes well is in y direction, so it only accounts for
                    # flow entering the well from matrix, not between well elements
                    if wellid == conn[jj][0] and conn[jj][2] != 2:
                        for kk in range(len(dir)):
                            if conn[jj][1] == elem[map[dir[kk]][ii,0]][0] and elem[map[dir[kk]][ii,0]][1] != elem[ii][1]:
                                allid[elem[ii][1]].append(conn[jj][0]+conn[jj][1])
                                count += 1
                    elif wellid == conn[jj][1] and conn[jj][2] != 2:
                        for kk in range(len(dir)):
                            if conn[jj][0] == elem[map[dir[kk]][ii,0]][0] and elem[map[dir[kk]][ii,0]][1] != elem[ii][1]:
                                allid[elem[ii][1]].append(conn[jj][0]+conn[jj][1])
                                count += 1
        return allid

    def getMatInterface(self, elem, conn, faces, allid):
        """ Get other user-defined interfaces

        Args:
            faces (list)    : Coordinates of each interface in the form of
                            [[name, [x1, x2, y1, y2, axis]],...]

        """
        axlst = [1, 2, 3]
        # Loop over all interfaces
        for ii in range(len(faces)):
            allid[faces[ii][0]] = []
            # Loop over all segments of one interface
            for jj in range(1,len(faces[ii])):
                plane = copy.deepcopy(axlst)
                plane.remove(faces[ii][jj][4])
                # Loop over all elements to find the interface connections
                for kk in range(self.Nele):
                    if elem[kk][plane[0]+3] >= faces[ii][jj][0] and elem[kk][plane[0]+3] < faces[ii][jj][1]:
                        if elem[kk][plane[1]+3] >= faces[ii][jj][2] and elem[kk][plane[1]+3] < faces[ii][jj][3]:
                            if abs(elem[kk][faces[ii][jj][4]+3] - faces[ii][jj][5]) < 0.5 * self.delta[faces[ii][jj][4]-1]:
                                # print(elem[kk], faces[ii][jj], self.delta[faces[ii][jj][4]-1])
                                for cc in range(self.Ncon):
                                    # Each element is only count once (unlike for wells)
                                    if elem[kk][0] == conn[cc][0] and conn[cc][2] == faces[ii][jj][4]:
                                        allid[faces[ii][0]].append(conn[cc][0]+conn[cc][1])
        return allid

    def updateINPUT(self, finput, allid, wellid=None):
        """ Update the TOUGH INPUT

        This function writes a new TOUGH INPUT file. It adds/replaces the following
        information:
            1 - Adds index of the INTERFACE block
            2 - Write the actual number of elements and connections
        Other info needs to be updated in INPUT should be added here in the future.

        Args:
            finput (str)    : Name of an existing TOUGH INPUT template
            allid (dict)    : All CONNE id for all wells, obtained from 'getInterface'
        """
        fin = open(finput, 'r')
        fout = open('INPUT', 'w')
        iwell = 1
        fkeys = allid.keys()
        Nface = len(fkeys)
        for line in fin:
            if 'INTERFACE' in line:
                fout.write(line)
                fout.write(f"&Interface_General_Info number_of_interfaces = {Nface} /\n")
                for item in fkeys:
                    Nc = len(allid[item])
                    fout.write(f"        &Individual_Interface_Specifics interface_name = '{item}',\n")
                    fout.write(f"                                        number_of_surfaces = 1,\n")
                    fout.write(f"                                        sign_of_flow_direction = 'DIR' /\n")
                    fout.write(f"                &Surface_Specifics      definition_mode = 'NameList', \n")
                    fout.write(f"                                        number_of_connections = {Nc}, \n")
                    fout.write(f"                                        format_to_read_data = '(A10)' /\n")
                    for jj in range(len(allid[item])):
                        fout.write(allid[item][jj]+'\n')
                    fout.write('\n')
            elif 'Conx_Time_Series' in line:
                fout.write(line)
                if wellid is None:
                    raise ValueError('If use Conx_Time_Series, user must provide well id!')
                item = wellid
                for jj in range(len(allid[item])):
                    fout.write(allid[item][jj]+'\n')
                fout.write('\n')
            elif 'Cartesian' in line:
                fout.write(f"'Cartesian'   {self.Nele}   {self.Ncon}   5   .FALSE.   .FALSE.\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()

    def wellMask(self, elem, maskzone):
        """ Mask the well elements

        This function is used to remove well elements from INCON, such that
        a different initial condition can be applied to well in INDOM block.

        Args:
            elem (list) : Obtained from 'adjustMesh'

        Returns:
            mask (nparray)  : Indicates if an element is well or not
        """
        mask = np.zeros((self.Nele, 1), dtype=bool)
        for ii in range(self.Nele):
            if elem[ii][1] in self.well:
                mask[ii,0] = True
            elif elem[ii][1] in maskzone:
                mask[ii,0] = True
        return mask

    def fracMask(self, elem, isfrac, map, frac_zone, mat_zone):
        """ Mask the frac zones where no geos data is provided

        This function is useful when the fracture from geos has irregular
        shape. Since MeshMaker only creates rectangular fracture, this function
        helps to convert blocks back to the matrix zone when it does not belong
        to a fracture

        Args:
            elem (list) : The ELEME
            isfrac (bool) : Whether or not an elem belongs to fracture
            zone (list) : Name of fracture zones
            map (list) : Neighbor map for elem

        Returns:
            elem (list) : elem with new zonal distribution

        """
        for ii in range(self.Nele):
            if elem[ii][1] in frac_zone and not isfrac[ii]:
                elem[ii][1] = mat_zone
                jj = ii
                while map['jP'][jj,0] != -1 and elem[map['jP'][jj,0]][1] not in frac_zone:
                    if elem[map['jP'][jj,0]][1] != mat_zone:
                        elem[map['jP'][jj,0]][1] = mat_zone
                    jj = map['jP'][jj,0]
                jj = ii
                while map['jM'][jj,0] != -1 and elem[map['jM'][jj,0]][1] not in frac_zone:
                    if elem[map['jM'][jj,0]][1] != mat_zone:
                        elem[map['jM'][jj,0]][1] = mat_zone
                    jj = map['jM'][jj,0]
        return elem



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

    def writeNewMesh(self, elem, conn, isbound, bcname='I'):
        """ Write the new mesh file

        Args:
            elem, conn (list): Obtained from 'adjustMesh()'
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
            if isbound[ii] == 1:
                fid.write(' I\n')
            elif isbound[ii] > 1:
                if self.useBC == 1:
                    # wname = bcname[isbound[ii][0]-2]
#                    wname = bcname
#                    fid.write(f" {wname}\n")
                    # Oil code does not support V flag, so use I for all wells, ZhiLi20200421
                    fid.write(' I\n')
                else:
                    fid.write(' I\n')
            else:
                fid.write('     \n')
        fid.write('     \n')
        # write connections
        fid.write('CONNE\n')
        for ii in range(self.Ncon):
            line = conn[ii]
            line[2] = str(line[2])
            for jj in range(3,6):
                if isinstance(line[jj],tuple):
                    print(line)
                # print(line[jj])
                line[jj] = self.reformat("{:.4E}".format(line[jj]))
            fid.write(f"{line[0]}{line[1]}                   {line[2]}{line[3]}{line[4]}{line[5]}")
            if len(line[6]) > 1:
                fid.write(line[6]+'\n')
            else:
                fid.write(line[6])
        fid.write('     \n')
        fid.close()
        sp.call(['rm', 'MESH'])
        # sp.call(['rm', 'meshinput'])
        sp.call(['mv', 'MESH_NEW', 'MESH'])
        return 0
