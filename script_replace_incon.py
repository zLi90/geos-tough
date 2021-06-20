"""
    Replace HF elements in INCON by closed permeability
"""
import numpy as np

"""
    Background permeability for the matrix
"""
state = [0.07, 1e-20]

"""
    Load INCON and SAVE files
"""
def load_file(fname):
    lines = []
    with open(fname) as fid:
        for line in fid:
            lines.append(line)
    return lines

lst_incon = load_file('INCON')
lst_save = load_file('SAVE')

"""
    Create dict of fracture elements
"""
def get_hf_index(incon):
    out = {}
    for elem in incon:
        if elem[0:5] != '     ' and elem[0:5] != 'INCON' and elem[0:5] != '<<<':
            out[elem[0:5]] = float(elem[71:86])
            # if elem[0].isupper() and elem[3:5] == '17':
                # out[elem[0:5]] = float(elem[71:86])
    return out
perm = get_hf_index(lst_incon)

"""
    Insert INCON into SAVE
"""
def insert_close_perm(perm, incon):
    out = []
    for elem in incon:
        elem = elem.split(' ')
        if elem[0] != '':
            if 'AqO:' in elem:
                id = elem[0]
                kx = float(elem[34])
                ky = float(elem[35])
                kz = float(elem[36])
                if id in perm:
                    elem[34] = "{:.8E}".format(perm[id])
                    elem[35] = "{:.8E}".format(ky/2.0)
                    elem[36] = "{:.8E}".format(perm[id]) + '\n'
                else:
                    elem[34] = "{:.8E}".format(kx/2.0)
                    elem[35] = "{:.8E}".format(ky/2.0)
                    elem[36] = "{:.8E}".format(kz/2.0) + '\n'
            elif 'AOG:' in elem:
                id = elem[0]
                kx = float(elem[30])
                ky = float(elem[31])
                kz = float(elem[32])
                if id in perm:
                    elem[30] = "{:.8E}".format(perm[id])
                    elem[31] = "{:.8E}".format(ky/2.0)
                    elem[32] = "{:.8E}".format(perm[id]) + '\n'
                else:
                    elem[30] = "{:.8E}".format(kx/2.0)
                    elem[31] = "{:.8E}".format(ky/2.0)
                    elem[32] = "{:.8E}".format(kz/2.0) + '\n'
        for ii in range(len(elem)):
            if elem[ii] == '':
                elem[ii] = ' '
            elif ii != 0:
                elem[ii] = ' '+elem[ii]
        if elem[0] == ' ':
            elem[0] = ''
        out.append(elem)
    return out
new_incon = insert_close_perm(perm, lst_save)

"""
    Save output
"""
fid = open('INCON_NEW','w')
for elem in new_incon:
    joined = ''.join([item for item in elem])
    # print(elem)
    for ii in range(len(elem)):
        # print(elem[ii])
        fid.write(elem[ii])
        # fid.write(joined)
fid.close()
