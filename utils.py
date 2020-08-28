""" Utility functions (mainly reading/writing functions) for INCON generator
    ZhiLi20191104
"""
import numpy as np
import xlrd
import time
import datetime

def writeIncon(elem, state, data, mask, perm0):
    """ Write data into INCON file

    Args:
        state (str)   : Model state of TOUGH simulation
        data (nparray): Values of INCON variables
        mask (nparray): Indicates if an element belongs to the well
    """
    fi = open("INCON","w")
    fi.write("INCON\n")
    for ii in range(0,len(elem)):
        if mask[ii,0] == False:
            eid = elem[ii][0]
            poro = "{:.8E}".format(data[ii,0])
            perm = "{:.8E}".format(data[ii,1])
            pres = "{:.3E}".format(data[ii,2])
            satu = "{:.3E}".format(data[ii,3])
            temp = "{:.3E}".format(data[ii,4])
            perm0str = "{:.8E}".format(perm0)
            # line1: element, porosity(3), state(4), permeability(5-7)
            fi.write(f"{eid}           {poro}  {state}                                     {perm} {perm0str} {perm}")
            fi.write('\n')
            # line2: primary variables (2,4,6)
            fi.write(f"           {pres}           {satu}           {temp}")
            fi.write('\n')
    fi.write("<<<")
    fi.close()

def writeInconOil(elem, state, data, mask, perm0, map, wellid, ic):
    """ Write data into INCON file (for the Hyd+Oil code) """
    num_blocks = 0
    fi = open("INCON","w")
    fi.write("INCON\n")
    for ii in range(0,len(elem)):
        if mask[ii,0] == False:
            eid = elem[ii][0]
            poro = "{:.8E}".format(data[ii,0])
            perm = "{:.8E}".format(data[ii,1])
            pres = "{:.3E}".format(data[ii,2])
            satu = "{:.3E}".format(data[ii,3])
            sato = "{:.3E}".format(data[ii,4])
            temp = "{:.3E}".format(data[ii,5])
            perm0str = "{:.8E}".format(perm0)
            # line1: element, porosity(3), state(4), permeability(5-7)
            fi.write(f"{eid}           {poro}  {state}                                     {perm} {perm0str} {perm}")
            fi.write('\n')
            # line2: primary variables (2,4,6)
            fi.write(f"           {pres}           {satu}           {sato}           {temp}")
            fi.write('\n')
            num_blocks += 1
        else:
            if elem[ii][1] in wellid:
                if map['iP'][ii,0] > 0 and mask[map['iP'][ii,0],0] == False:
                    eid = elem[ii][0]
                    fid = map['iP'][ii,0]
                    wstate = ic[elem[ii][1]][0]
                    poro = "{:.8E}".format(1.0)
                    perm = "{:.8E}".format(data[fid,1])
                    pres = "{:.3E}".format(ic[elem[ii][1]][1])
                    satu = "{:.3E}".format(ic[elem[ii][1]][2])
                    sato = "{:.3E}".format(ic[elem[ii][1]][3])
                    temp = "{:.3E}".format(ic[elem[ii][1]][4])
                    perm0str = "{:.8E}".format(perm0)
                    fi.write(f"{eid}           {poro}  {wstate}                                     {perm} {perm0str} {perm}")
                    fi.write('\n')
                    fi.write(f"           {pres}           {satu}           {sato}           {temp}")
                    fi.write('\n')
                    num_blocks += 1
                elif map['iM'][ii,0] > 0 and mask[map['iM'][ii,0],0] == False:
                    eid = elem[ii][0]
                    fid = map['iM'][ii,0]
                    wstate = ic[elem[ii][1]][0]
                    poro = "{:.8E}".format(1.0)
                    perm = "{:.8E}".format(data[fid,1])
                    pres = "{:.3E}".format(ic[elem[ii][1]][1])
                    satu = "{:.3E}".format(ic[elem[ii][1]][2])
                    sato = "{:.3E}".format(ic[elem[ii][1]][3])
                    temp = "{:.3E}".format(ic[elem[ii][1]][4])
                    perm0str = "{:.8E}".format(perm0)
                    fi.write(f"{eid}           {poro}  {wstate}                                     {perm} {perm0str} {perm}")
                    fi.write('\n')
                    fi.write(f"           {pres}           {satu}           {sato}           {temp}")
                    fi.write('\n')
                    num_blocks += 1
    fi.write("<<<")
    fi.close()
    print("Total number of grid blocks in INCON = ",num_blocks)


def excel_date(date1):
    """ Convert date string to excel date number

    Args:
        date1 (str): Date string in the format of yyyymmdd (later than 19991230)

    Returns:
        Excel date number
    """
    temp = time.mktime(datetime.datetime.strptime('19991230','%Y%m%d').timetuple())
    base = 36525
    delta = date1 - temp
    return base + delta/86400

def readBC(fname, bc):
    """ Read time-variable boundary condition from excel

    This function is designed specifically for HFTS Production Data (in excel).
    For other types of input data, should define other functions.

    Args:
        fname (str)     : Excel file name to read
        bc (dict)       : Dict containing bc settings, which includes:
            isheet (int)    : Index of sheet to read
            period (list)   : [start end] date in the format of yyyymmdd
            col (list)      : Columns in excel to be read
            units (list)    : Unit conversion function for each column

    Returns:
        data (nparray)  : Data extracted from excel file within the range of 'period'
    """
    # top 4 lines are not date numbers, so it is removed
    offset = 13
    # load excel file
    datafile = xlrd.open_workbook(fname)
    for ss in range(len(bc['isheet'])):
        sheet = datafile.sheet_by_index(bc['isheet'][ss])
        times = sheet.col_values(0)
        kk = 0
        while kk < len(times):
            if type(times[kk]) != float:
                del times[kk]
            else:
                kk += 1
        times = np.array(times, dtype=float)
        # find row index corresponds to the desired time period
        t0 = time.mktime(datetime.datetime.strptime(bc['period'][0],'%Y%m%d').timetuple())
        t1 = time.mktime(datetime.datetime.strptime(bc['period'][1],'%Y%m%d').timetuple())
        ind0 = np.where(abs(times - excel_date(t0)) < 0.1)[0][0] + offset - 1
        ind1 = np.where(abs(times - excel_date(t1)) < 0.1)[0][0] + offset - 1
        # extract data into numpy array
        data = np.zeros((ind1-ind0,len(bc['col'])))
        counter = 0
        for kk in range(ind0, ind1):
            for ii in range(len(bc['col'])):
                if isinstance(bc['col'][ii], list):
                    # calculate water saturation
                    if type(sheet.cell_value(kk,bc['col'][ii][0])) is float and type(sheet.cell_value(kk,bc['col'][ii][1])) is float:
                        comp1 = unitConv(sheet.cell_value(kk,bc['col'][ii][0]), bc['units'][ii][0])
                        comp2 = unitConv(sheet.cell_value(kk,bc['col'][ii][1]), bc['units'][ii][1])
                        data[counter,ii] = comp2 / (comp1 + comp2)
                    else:
                        data[counter,ii] = np.nan
                else:
                    # extract pressure and temperature
                    if type(sheet.cell_value(kk,bc['col'][ii])) is float:
                        data[counter,ii] = sheet.cell_value(kk,bc['col'][ii])
                    else:
                        data[counter,ii] = np.nan
                    data[counter,ii] = unitConv(data[counter,ii], bc['units'][ii])
            counter += 1
        data[:,0] -= data[0,0]
        # print boundary conditions
        print('\nPrint boundary conditions for :', bc['name'][ss])
        for kk in range(data.shape[0]):
            print("{:.4E}".format(data[kk,0]), "{:.4E}".format(data[kk,1]), "{:.4E}".format(data[kk,2]), "{:.4E}".format(data[kk,3]))
        print('Total number of BC rows: ',data.shape[0])
    return data


def readBCandAvg(fname, bc):
    """ Read time-variable boundary condition from excel

    Read BC from well production data and average to get piecewise-constant
    boundary conditions. This is used for simulation without BOUNDARIES block.
    Note that structure of input variable 'bc' is different from previous
    'readBC' function.

    Args:
        fname (str)     : Excel file name to read
        bc (dict)       : Dict containing bc settings, which includes:
            useBC (int)     : Use BC or not
            name (list)     : Name of BC, must match media name in INDOM
            state (list)    : State index of each BC
            isheet (list)   : Index of sheet to read (variable index), can be None
            period (list)   : [start end] date in the format of yyyymmdd
            col (list)      : Columns in excel to be read (well index), can be None
            units (list)    : Unit conversion function for each column

    Returns:
        data_out (list)     : List of BC for all boundaries averaged over 'period'
    """
    #   top 4 lines are not date numbers, so it is removed
    offset = 13
    #   load excel file
    datafile = xlrd.open_workbook(fname)
    #   Get time vector
    sheet = datafile.sheet_by_index(bc['isheet'][0])
    times = sheet.col_values(0)
    kk = 0
    while kk < len(times):
        if type(times[kk]) != float:
            del times[kk]
        else:
            kk += 1
    times = np.array(times, dtype=float)
    t0 = time.mktime(datetime.datetime.strptime(bc['period'][0],'%Y%m%d').timetuple())
    t1 = time.mktime(datetime.datetime.strptime(bc['period'][1],'%Y%m%d').timetuple())
    ind0 = np.where(abs(times - excel_date(t0)) < 0.1)[0][0] + offset - 1
    ind1 = np.where(abs(times - excel_date(t1)) < 0.1)[0][0] + offset - 1
    #   Extract BC data from excel file
    data_out = []
    print(bc['name'])

    for ibc in range(len(bc['name'])):
        name = bc['name'][ibc]
        state = bc['state'][ibc]
        data = np.zeros((ind1-ind0,len(bc['col'])))
        for ss in range(len(bc['isheet'])):
            if bc['isheet'][ss] is not None:
                counter = 0
                sheet = datafile.sheet_by_index(bc['isheet'][ss])
                values = sheet.col_values(bc['col'][ss][ibc])
                for kk in range(ind0, ind1):
                    if type(values[kk]) is float:
                        data[counter,ss] = values[kk]
                    else:
                        data[counter,ss] = np.nan
                    data[counter,ss] = unitConv(data[counter,ss], bc['units'][ss])
                    counter += 1
            else:
                data[:,ss] = np.nan
        data[:,0] -= data[0,0]
        #   Calculate averages
        data_avg = [name, state]
        for ii in range(data.shape[1]-1):
            data_avg.append(np.nanmean(data[:,ii+1]))
        data_out.append(data_avg)

    return data_out


def genConstIC(fname, dim, delta, yfrac, value, margin):
    """ Generate constant initial condition for fractures

    This function is used to create fake GEOS output, i.e., in the same format
    with GEOS output, but with constant values.

    Args:
        fname (str)     : Name of the output file saved
        dim (list)      : Dimension of domain
        delta (list)    : Element increments of domain
        yfrac (nparray) : y-coordinates of fractures
        value (float)   : Constant field value to be assigned to the fractures
        margin (list)   : Distance (in number of elements) from fracture edge to
                            domain boundary
    """
    aper = np.zeros((len(yfrac)*(dim[0]-margin[0])*(dim[2]-margin[2]),4))
    ll = 0
    for ii in range(dim[0]-margin[0]):
        for jj in range(dim[2]-margin[2]):
            for kk in range(len(yfrac)):
                aper[ll,0] = ii * delta[0] + 0.5*delta[0]
                aper[ll,1] = yfrac[kk]
                aper[ll,2] = -(dim[2] * delta[2] - jj * delta[2] - 0.5*delta[2])
                aper[ll,3] = value
                ll += 1
    np.savetxt(fname, aper, delimiter=",")
    return 0

def unitConv(value, func):
    """ Unit conversion functions """
    if func == 'none':
        cfun = lambda x: x
    elif func == 'psi2Pa':
        cfun = lambda x: x * 6894.76
    elif func == 'd2sec':
        cfun = lambda x: x * 86400.0
    elif func == 'F2C':
        cfun = lambda x: (x - 32.0) * 5.0 / 9.0
    elif func == 'bbl2kg':
        cfun = lambda x: x / 0.00839
    elif func == 'mcf2kg':
        cfun = lambda x: x / 0.023
    elif func == 'bbl2m3':
        cfun = lambda x: x * 0.11924
    elif func == 'mcf2m3':
        cfun = lambda x: x * 28.3168
    else:
        raise ValueError('Undefined unit conversion!')
    return cfun(value)


def readTOUGH(time, vid):
    """ Read initial values from TOUGH output

    This function needs further developing and debugging! Ignore it for now!
    """
    fid = open('Plot_Data_Elem', 'r')
    title = str.split(fid.readline())
    counter = 0
    while True:
        counter += 1
        line = str.split(fid.readline())
        if counter > 5e6:
            raise ValueError("'time' is not found in the first 5e6 lines of datafile!")
        if len(line) > 0:
            if line[0] == "ZONE":
                t = line[3]
                t = t[:-2]
                # search for output data at 'time'
                if float(t) == time:
                    while len(line) > 0:
                        line = str.split(fid.readline())
                        data.append(float(line[vid]))
                    break
    fid.close()
    return data
