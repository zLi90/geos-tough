""" Remove MINC fracture from Plot_Data_Elem """
import numpy as np

# fdir = 'sens_ntwk3/upslo2_gor2800_150/'
fdir = 'sens_ntwk3/upslo2_movie_hf1500/'


with open(fdir+'Plot_Data_Elem') as fid:
    lines = []
    for line in fid:
        if len(line) > 11:
            id = line[11]
            if not id.islower():
                lines.append(line)
        else:
            lines.append(line)


fid = open('PDE_NEW','w')
for ii in range(len(lines)):
    # print(elem[ii])
    fid.write(lines[ii])
fid.close()


with open(fdir+'MESH') as fid:
    lines = []
    for line in fid:
        if len(line) > 11:
            id = line[0]
            id2 = line[5]
            if not id.islower():
                if not id2.islower():
                    lines.append(line)
        else:
            lines.append(line)


fid = open('MESH_NEW','w')
for ii in range(len(lines)):
    # print(elem[ii])
    fid.write(lines[ii])
fid.close()
