""" Create mapping from TH output (Plot_Data_Elem) to 3D Cartesian np array """
import numpy as np

def readConne(fname, dim):
    """ Read CONNE from TH MESH """
    conn = []
    flag = 0
    N = (dim[0]-1) * dim[1] * dim[2] + \
        (dim[1]-1) * dim[2] * dim[0] + \
        (dim[2]-1) * dim[0] * dim[1]
    fid = open(fname,'r')
    while True:
        if flag == 1:
            break
        line = fid.readline()
        if len(line) >= 5 and line[0:5] == 'CONNE':
            # Get connections
            for ii in range(N):
                line = fid.readline()
                # each element in conn has format [conn1, conn2, direction]
                conn.append([line[0:5], line[5:10], line[29]])
                if ii == N-1:
                    flag = 1
    fid.close()
    return conn

def readEleme(fname, dim):
    """ Read ELEME from TH MESH """
    elem = []
    coord = []
    flag = 0
    N = dim[0]*dim[1]*dim[2]
    fid = open(fname,'r')
    while True:
        if flag == 1:
            break
        line = fid.readline()
        if len(line) >= 5 and line[0:5] == 'ELEME':
            # Get mesh ID
            ii = 0
            while True:
                line = fid.readline()
                if line[0].isupper() == True:
                    elem.append(line[0:5])
                    coord.append([float(line[50:60]), float(line[60:70]), float(line[70:80])])
                    ii += 1
                if ii == N:
                    flag = 1
                    break
    fid.close()
    return elem, coord



def getLstValue(ind, trigger, i_current, nmax=100):
    lst = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
        'U','V','W','X','Y','Z','end']
    if i_current == 3:
        ind = int(ind) + 1
        if ind >= nmax:
            ind = 0
            trigger[i_current-1] = 1
        ind = "{:02d}".format(ind)
    else:
        if trigger[i_current] == 1:
            lst_id = lst.index(ind)
            ind = lst[lst_id+1]
            trigger[i_current] = 0
            if ind == 'end':
                ind = '0'
                trigger[i_current-1] = 1
    return trigger, ind



def buildMap(dim, order=[2,0,1]):
    """ Build map from MESH to 3D np array """
    map = ['A0000']
    ind0 = 'A'
    ind1 = '0'
    ind2 = '0'
    ind3 = '00'
    trigger = [0,0,0,0]
    # for 3D problem
    if dim[order[0]] > 1:
        for ii in range(dim[order[2]]):
            for jj in range(dim[order[1]]):
                for kk in range(dim[order[0]]):
                    # 4th and 5th index of the element id
                    trigger, ind3 = getLstValue(ind3, trigger, 3, dim[order[0]])
                    # 3rd index
                    trigger, ind2 = getLstValue(ind2, trigger, 2)
                    # 2rd index
                    trigger, ind1 = getLstValue(ind1, trigger, 1)
                    # 1st index
                    trigger, ind0 = getLstValue(ind0, trigger, 0)
                    # combine all index to form element id
                    eid = ind0 + ind1 + ind2 + ind3
                    map.append(eid)
    # for 2D problem
    else:
        for ii in range(dim[order[2]]):
            for jj in range(dim[order[1]]):
                # 4th and 5th index of the element id
                trigger, ind3 = getLstValue(ind3, trigger, 3, dim[order[1]])
                # 3rd index
                trigger, ind2 = getLstValue(ind2, trigger, 2)
                # 2rd index
                trigger, ind1 = getLstValue(ind1, trigger, 1)
                # 1st index
                trigger, ind0 = getLstValue(ind0, trigger, 0)
                # combine all index to form element id
                eid = ind0 + ind1 + ind2 + ind3
                map.append(eid)

    return map[:-1]


def backMap(eid, dim, order=[2,0,1]):
    """ From element id to element index """
    lst = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
        'U','V','W','X','Y','Z']
    ielem = 0
    # 3D problems
    if dim[order[0]] > 1:
        for ii in range(4):
            ind = lst.index(eid[ii])
            if ii == 0:
                ielem += (ind - 10) * dim[order[0]] * len(lst) * len(lst)
            elif ii == 1:
                ielem += ind * dim[order[0]] * len(lst)
            elif ii == 2:
                ielem += ind * dim[order[0]]
            else:
                ielem += int(eid[ii:])
    # 2D problems
    else:
        for ii in range(4):
            ind = lst.index(eid[ii])
            if ii == 0:
                ielem += (ind - 10) * dim[order[1]] * len(lst) * len(lst)
            elif ii == 1:
                ielem += ind * dim[order[1]] * len(lst)
            elif ii == 2:
                ielem += ind * dim[order[1]]
            else:
                ielem += int(eid[ii:])
    return ielem


def get3Ddata(folders, dim, dim_order, fields, out_time):
    all_fields = ['x','y','z','P','T','S_hyd','S_aqu','S_gas','S_org','X_inh','k_rg','k_rw','k_ro','k_adj_F','perm_abs','porosity','Pcap(aq)','Pcap(go)']
    data3D = {}
    eid = []
    N = dim[0]*dim[1]*dim[2]
    col = []
    data = {}
    for ii in range(len(all_fields)):
        col.append(20 + ii*14)
    for folder in folders:
        elem, coord = readEleme(folder+'MESH', dim)
        fname = folder+'Plot_Data_Elem'
        times = []
        xVec = []
        yVec = []
        zVec = []
        data[folder] = {}
        ind = {}
        for field in fields:
            data[folder][field] = []
        # Load Plot_Data_Elem
        tt = -1

        with open(fname) as fid:
            for line in fid:
                if len(line) > 6 and line[0:6] == 'ZONE T':
                    tt += 1
                    times.append(float(line[10:24]))
                    # Get cell coordinates
                    ii = 0
                    while True:
                        line = fid.readline()
                        linelst = line.split()
                        if linelst[1][0].isupper() == True:
                            if tt == 0:
                                eid.append(linelst[1])
                            for field in fields:
                                ind = all_fields.index(field)
                                try:
                                    data[folder][field].append(float(line[col[ind]:col[ind]+14]))
                                except:
                                    data[folder][field].append(-100)
                            ii += 1
                        if ii == N:
                            break

        print('PDE data available at T = ', times)
        if len(out_time) == 1:
            t_ind = times.index(out_time[0])
        # Allocate data into 3D matrix
        data3D[folder] = {}
        for field in fields:
            data3D[folder][field] = np.zeros((dim[dim_order[2]],dim[dim_order[1]],dim[dim_order[0]],len(out_time)))
            # loop over all time slices
            elem_ind = []
            if len(out_time) == 1:
                data1D = np.zeros((N,1))
                for ll in range(N):
                    elem_ind.append(backMap(eid[ll], dim, dim_order))
                    data1D[elem_ind[ll]] = data[folder][field][N*t_ind+ll]
                data3D[folder][field][:,:,:,0] = np.reshape(data1D,(dim[dim_order[2]],dim[dim_order[1]],dim[dim_order[0]]),order='C')
            else:
                for tt in range(len(out_time)):
                    # loop over all grid cells
                    data1D = np.zeros((N,1))
                    for ll in range(N):
                        if tt == 0:
                            elem_ind.append(backMap(eid[ll], dim, dim_order))
                        data1D[elem_ind[ll]] = data[folder][field][N*tt+ll]
                    data3D[folder][field][:,:,:,tt] = np.reshape(data1D,(dim[dim_order[2]],dim[dim_order[1]],dim[dim_order[0]]),order='C')
    return data3D

def getSliceID(elem, coord, xlim):
    eid = {}
    cid = []
    kk = 0
    # extract element id and coordinates for a sub-region of the domain
    for ii in range(len(coord)):
        if coord[ii][0] >= xlim[0] and coord[ii][0] < xlim[1]:
            if coord[ii][1] >= xlim[2] and coord[ii][1] < xlim[3]:
                if coord[ii][2] >= xlim[4] and coord[ii][2] < xlim[5]:
                    cid.append([kk, coord[ii][0], coord[ii][1], coord[ii][2]])
                    eid[kk] = elem[ii]
                    kk += 1
    # reorder based on coordinates
    cid = np.array(cid)

    cid = cid[cid[:,1].argsort()]
    cid = cid[cid[:,2].argsort(kind='mergesort')]
    cid = cid[cid[:,3].argsort(kind='mergesort')]
    out = []
    for ii in range(np.shape(cid)[0]):
        out.append(eid[cid[ii,0]])
    return out, cid

def getSlice(folders, eid, fields, dim, out_time):
    all_fields = ['x','y','z','P','T','S_hyd','S_aqu','S_gas','S_org','X_inh','k_rg','k_rw','k_ro','k_adj_F','perm_abs','porosity','Pcap(aq)','Pcap(go)']
    data3D = {}
    N = dim[0]*dim[1]*dim[2]
    col = []
    data = {}
    for ii in range(len(all_fields)):
        col.append(20 + ii*14)
    for folder in folders:
        fname = folder+'Plot_Data_Elem'
        times = []
        data[folder] = {}
        ind = {}
        for field in fields:
            data[folder][field] = np.zeros((len(eid)))
        # Load Plot_Data_Elem
        tt = -1

        with open(fname) as fid:
            for line in fid:
                if len(line) > 6 and line[0:6] == 'ZONE T':
                    tt += 1
                    times.append(float(line[10:24]))
                    # Get cell coordinates
                    ii = 0
                    if times[-1] == out_time[0]:
                        while True:
                            line = fid.readline()
                            linelst = line.split()
                            if linelst[1][0].isupper() == True:
                                if linelst[1] in eid:
                                    kk = eid.index(linelst[1])
                                    for field in fields:
                                        ind = all_fields.index(field)
                                        try:
                                            data[folder][field][kk] = float(line[col[ind]:col[ind]+14])
                                        except:
                                            data[folder][field][kk] = -100.0
                                ii += 1
                            if ii == N:
                                break
        print('PDE data available at T = ', times)
    return data

def from1Dto2D(data1D, dim, dx, offset, y_frac):
    for ll in range(np.shape(data1D)[0]):
        for ii in range(3):
            data1D[ll,ii] += offset[ii]
    data2D= np.nan * np.ones((dim[0], dim[2]))
    xVec = np.linspace(0.5*dx[0], dim[0]*dx[0]-0.5*dx[0], int(dim[0]))
    zVec = np.linspace(-dim[2]*dx[2]+0.5*dx[2], -0.5*dx[2], int(dim[2]))
    for ll in range(np.shape(data1D)[0]):
        if abs(y_frac - data1D[ll,1]) < 0.5*dx[1]:
            ii = abs(xVec - data1D[ll,0]).argmin()
            kk = abs(zVec - data1D[ll,2]).argmin()
            data2D[ii,kk] = data1D[ll,3]
    return data2D
