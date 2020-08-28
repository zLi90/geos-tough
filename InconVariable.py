""" Create data array (matches TOUGH MESH elements) for one TOUGH INCON variable.
    Data could be read from GEOS or could be constant.
    ZhiLi 20191101

    Todo:
        * Enable reading INCON variable from existing TOUGH outputs
"""
import numpy as np


class InconVariable():
    """ The class for reading one INCON variable

    Args:
        source (list)   : Contains user-defined settings for this variable with format:
                        [source of variable, file name, default value, unit conversion]
        Nele (int)      : Number of elements in the TOUGH MESH file
        isfrac (nparray): Whether the elements represent fractures or not
        out_domain (Bool): Whether or not output GEOS data as a separate file
        nswarm (int)    : Number of small fractures in a fracture swarm

    Attributes:
        source (str)    : Source of variable, could be 'geos' if it comes from GEOS
                          simulation, or 'const' if it is a constant, or 'binary'
                          if it has 2 constant values for matrix and fracture
                          respectively.
        fname (str)     : File name if data is read from an external file (source='geos')
        const (float)   : Default initial values for the matrix
        constF(float)   : Default initial values for the fractures
        convert (str)   : Unit conversion option
        Nele (int)      : Number of elements in TOUGH MESH
        domain (nparray): 1d array with size=Nele and values=const
        isfrac (nparray): Array indicates if an element represents a fracture or not

    """
    def __init__(self, source, Nele, isfrac, out_domain=False, nswarm=1):
        self.source = source[0]
        self.fname = source[1]
        self.const = source[2][0]
        self.constF = source[2][1]
        self.convert = source[3]
        self.Nele = Nele
        self.domain = self.const * np.ones((Nele, 1), dtype=float)
        self.isfrac = isfrac
        self.output_domain = out_domain
        self.nswarm = nswarm

    def loadGEOS(self):
        """ Load GEOS output file (in .csv format)

        Returns:
            Numpy array with row format [x, y, z, geos output value]
        """
        return np.genfromtxt(self.fname, delimiter=',')

    def convertGEOS(self, x, n):
        """ Perform unit conversion on GEOS outputs

        Args:
            data (nparray)  : Output from 'loadGEOS()'

        Returns:
            data (nparray)  : Same format with input but values are converted to SI units
                              If self.convert == 'aper2perm', it converts aperture to
                              permeability following cubic law.
        """
        if self.convert == 'none':
            cfun = lambda x: x
        elif self.convert == 'psi2Pa':
            cfun = lambda x: x * 6894.76
        elif self.convert == 'aper2perm':
            cfun = lambda x: (x**2.0)/12.0/(n*n)
        elif self.convert == 'F2C':
            cfun = lambda x: (x - 32.0) * 5.0 / 9.0
        else:
            raise ValueError('Undefined unit conversion!')
        # data[:,3] = cfun(data[:,3])
        y = cfun(x)
        return y

    def insertGEOS(self, data, gmap):
        """ Insert GEOS values into TOUGH mesh

        Args:
            data (nparray)  : Output from 'convertGEOS()'
            gmap (nparray)  : gmap[row index in 'geos'] = row index in TOUGH MESH

        Returns:
            out1d (nparray) : Initial values for all elements. Values in fractures
                              are read from GEOS. Values in matrix are set to default.
        """
        out1d = self.domain
        N = data.shape[0]
        out_data = []
        count = 0
        for ii in range(0,N):
            if not np.isnan(gmap[ii,0]):
                if not np.isnan(data[ii,3]):
                    out1d[gmap[ii,0]] = self.convertGEOS(data[ii,3], self.nswarm)
                    # out1d[gmap[ii,0]] = data[ii,3]
                    out_data.append([data[ii,:]])
                    count += 1
        out_data = np.reshape(out_data, (count,4))
        if self.output_domain is True and self.source == 'geos':
            self.writeInsertedData(out_data)
        return out1d

    def writeInsertedData(self, out):
        """ Write an output file for the GEOS data within the domain range

        """
        fname = self.fname+"_domain"
        print('Writing GEOS output for the domain as ', fname)
        fid = open(fname,"w")
        for ii in range(np.shape(out)[0]):
            fid.write(f"{out[ii,0]} {out[ii,1]} {out[ii,2]} {out[ii,3]}\n")
        fid.close()

    def targMean(self, out1d, targ):
        """ Adjust GEOS permeability to get a targeted mean value """
        if targ > 0:
            # k_targ = (targ ** 2.0) / 12.0
            k_targ = targ
            a_targ = (12.0 * k_targ) ** 0.5

            outfrac = []
            for ii in range(len(out1d)):
                if out1d[ii] != self.const and not np.isnan(out1d[ii]):
                    aper = (12.0 * out1d[ii]) ** 0.5
                    outfrac.append([ii,aper])
            outfrac_np = np.array(outfrac)
            avg = np.nanmean(outfrac_np[:,1])
            print('GEOS avg aperture = ',avg,', targ = ',a_targ)
            ratio = a_targ / avg
            outfrac_np[:,1] = outfrac_np[:,1] * ratio
            for ii in range(np.shape(outfrac_np)[0]):
                out1d[outfrac_np[ii,0]] = (outfrac_np[ii,1] ** 2.0) / 12.0
                # out1d[outfrac_np[ii,0]] = 1e-10

        return out1d


    def processField(self, gmap, targ=-1):
        """ Perform complete processing of one INCON variable

        Args:
            gmap (nparray)  : gmap[row index in 'geos'] = row index in TOUGH MESH

        Returns:
            data (nparray)  : 1d array containing final INCON values for all elements
        """
        if self.source == 'geos':
            geos = self.loadGEOS()
            # geos = self.convertGEOS(geos)
            data = self.insertGEOS(geos, gmap)
            # data = self.targMean(data, targ)
        elif self.source == 'binary':
            data = self.domain
            for ii in range(self.Nele):
                if self.isfrac[ii] == True:
                    data[ii] = self.constF
        elif self.source == 'const':
            data = self.domain
        else:
            raise ValueError('Undefined source of the INCON variable!')
        return data
