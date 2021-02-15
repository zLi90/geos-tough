""" Write INPUT file for TOUGH HYD+OIL simulations """
import os
import subprocess as sp
import numpy as np
import copy
from InconVariable import *
from BuildMesh import *
from utils import *


class Input():
    def __init__(self, elem, conn):
        self.Nele = len(elem)
        self.Ncon = len(conn)

    def writeMEMORY(self, finput):
        fin = open(finput, 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'Cartesian' in line:
                fout.write(f"'Cartesian'   {self.Nele}   {self.Ncon}   5   .FALSE.   .FALSE.\n")
            else:
                fout.write(line)
        fin.close()
        fout.close()

    def writeROCKS(self, rocks):
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'ROCKS' in line:
                fout.write(line)
                for reg_name in rocks.keys():
                    param = rocks[reg_name][1][1]
                    rho = "{:4.1E}".format(param[0])
                    phi = "{:4.1E}".format(param[1])
                    perm = "{:8.2E}".format(param[2])
                    comp = "{:5.1E}".format(param[3])
                    ind_rp = param[4]
                    ind_cp = param[5]
                    # get rp and cp coefficients
                    rpCoeff = self.relPermParam(ind_rp)
                    cpCoeff = self.capPresParam(ind_cp)
                    # write lines
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

    def writePARAM(self, time, ic):
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

    def writeINDOM(self, databc):
        sp.call(['mv', 'INPUT', 'INPUT0'])
        fin = open('INPUT0', 'r')
        fout = open('INPUT', 'w')
        for line in fin:
            if 'INDOM' in line:
                fout.write(line)
                for elem in databc:
                    reg_name = elem[0]
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

    def writeTIMES(self, tVec):
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

    def writeCONX(self, conx):
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
        elif num == 6:
            param = [0.05,0.05,0.0,2.0]
        else:
            raise ValueError('Relative permeability equation number is not available!!!')
        return param

    def capPresParam(self, num):
        if num == 5:
            param = []
        elif num == 8:
            param = [0.4, 1.85, 10.0, 11.0, 0.21]
        else:
            raise ValueError('Capillary pressure equation number is not available!!!')
        return param
