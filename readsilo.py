""" Read from GEOS outputs (.silo) and save them into .csv file (ZhiLi20191101)

Inputs:
    fname     : Name of .silo file
    ftype     : Field type in .silo file that will be parsed
    fields    : All field names in 'ftype'
    subfolder : Directory contains GEOS output from all processors
    nproc     : Number of processors used in GEOS simulation
    saveFields: Field name and file name of the variable to be saved
                    (field name must match names listed in 'fields')
"""
import numpy as np
import silo_parser as sp
from gplot import *

#
#   User settings
#
# fname = 'FiveFrac_base_737154'
# subfolder = 'GEOS_upscaling/A_Upscaling/sub'

# fname = 'FiveFrac_base_258004'
# subfolder = 'GEOS_upscaling/B_Non_Upscaling/sub'

# fname = 'plot_703036'
# subfolder = 'sub'

fname = 'FiveFrac_base_687492'
subfolder = 'GEOS_upscaling202008/Case_2/Case_2/HighLandingPoint'
# subfolder = 'GEOS_upscaling202008/Case_1/Case_1'

ftype = 'FaceFields'
fields = ['Aperture','FaceArea','Pressure','ProppantVolumeFraction',
          'Volume','birthTime',
          'permeability','stressNOnFace','totalLeakedThickness']

nproc = 16
save = True
saveFields = {
                'Pressure'  :   'GEOS_upshi_P.csv',
                'Aperture'  :   'GEOS_upshi_A.csv',
                'permeability': 'GEOS_upshi_Perm.csv',
                'totalLeakedThickness': 'GEOS_upshi_Leak.csv'
             }
# saveFields = {}

plotFields = ['Aperture']
# plotFields = fields

#
#   Functions
#

def main():
    """Main function converts .silo to .csv

    It reads the .silo file using 'silo_parser.parse_file', extract
    the variables defined in 'saveFields', and save them in .csv files.

    """
    data = sp.parse_file(fname, ftype, fields, subfolder, nproc)
    # from IPython import embed; embed()
    print('Available variable names to be extracted are: ', data.keys())
    for key in plotFields:
        makePlot(data, key)
    for key in saveFields:
        extractField(data, key, save, saveFields[key])

def extractField(data, field, savefile, savename):
    """Function to extract one field from .silo file

    It creates a .csv file with columns [x, y, z, value]

    Args:
        data (dict): GEOS outputs from 'silo_parser.parse_file()'
        field (str): Name of field to be extracted, which matches those in 'fields'
        savename (str): Name of .csv file that will be generated
    """
    print('Extracting '+field)
    out = np.r_['1,2,0',data['FaceCenter[0]'],data['FaceCenter[1]']
        ,data['FaceCenter[2]'],data[field]]
    if savefile == True:
        np.savetxt(savename, out, delimiter=",")

# Call main function
if __name__=="__main__":
    main()
