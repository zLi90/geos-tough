""" Plot relative permeability functions used in TOUGH """
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# residual saturation
Srw = 0.25
Srg = 0.01
Sro = 0.15

Sfw = 0.05
Sfg = 0.0
Sfo = 0.05

# exponents
a1 = 4.0
a2 = 2.0
a1f = 1.0
a2f = 1.0
b2 = 1.85
b3 = 10.0
b4 = 11.0

# line color
co = ['b','r','k']
fs = 8

# relative permeability functions

def rpw(Sw, Srw, a1):
    if Sw > Srw:
        return ((Sw-Srw)/(1.0-Srw))**a1
    return 0.0

def rpg(Sg, Srg, Srw, a2):
    if Sg > Srg:
        return ((Sg-Srg)/(1.0-Srw))**a2
    return 0.0

def rpo(So, Sro, Srw, a1):
    if So > Sro:
        temp = 1.0 - Srw - Sro
        return ((So - Sro)/temp)**a1
    return 0.0

# capillary pressure functions
def pwo(Sw, So, Srw, Sro, n, b):
    rho = 1000.0
    g = 9.81
    invn = 1.0/n
    invm = 1.0 / (1.0 - 1.0 / n)
    Swbar = (Sw - Srw) / (1.0 - Srw)
    coef = rho * g / b
    if Swbar > 1.0:
        Swbar = 1.0
    elif Swbar < 0.0:
        Swbar = 0.0
    if Swbar > 0.1:
        return -coef * (Swbar**(-invm) - 1.0)**invn
    else:
        part1 = -coef * (0.1**(-invm) - 1.0)**invn
        part2 = coef * ((0.10**(-invm)-1.0)**invn - (0.1**(-invm)-1.0)**invn)
        part3 = 0.1 - Swbar
        return part1 - 500.0 * part2 * part3

def pgo(Sw, So, Srw, Sro, n, b):
    rho = 1000.0
    g = 9.81
    invn = 1.0/n
    invm = 1.0 / (1.0 - 1.0 / n)
    Stbar = (Sw + So - Srw - Sro) / (1.0 - Srw - Sro)
    coef = rho * g / b
    if Stbar > 1.0:
        Stbar = 1.0
    elif Stbar < 0.0:
        Stbar = 0.0
    if Stbar > 0.1:
        return -coef * (Stbar**(-invm) - 1.0)**invn
    else:
        part1 = -coef * (0.1**(-invm) - 1.0)**invn
        part2 = coef * ((0.1**(-invm)-1.0)**invn - (0.1**(-invm)-1.0)**invn)
        part3 = 0.1 - Stbar
        return part1 - 500.0 * part2 * part3



# calculate relative permeabilities
s = np.linspace(0,1,1000)
kr = []
kf = []
cp = []
for ii in range(len(s)):
    #
    #   relative permeability
    #
    krw = 0.0
    krg = 0.0
    kro = 0.0
    kfw = 0.0
    kfg = 0.0
    kfo = 0.0
    # water
    if s[ii] > Srw:
        krw = np.minimum(((s[ii]-Srw)/(1.0-Srw))**a1, 1.0)
    if s[ii] > Sfw:
        kfw = np.minimum(((s[ii]-Sfw)/(1.0-Sfw))**a1f, 1.0)
    # gas
    if s[ii] > Srg:
        krg = np.minimum(((s[ii]-Srg)/(1.0-Srw))**a2, 1.0)
    if s[ii] > Sfg:
        kfg = np.minimum(((s[ii]-Sfg)/(1.0-Sfw))**a2f, 1.0)
    # oil
    if s[ii] > Sro:
        temp = 1.0 - Srw - Sro
        kro = np.minimum(((s[ii] - Sro)/temp)**a1, 1.0)
        if abs(s[ii]-0.6) < 1e-3:
            print('krg, kro = ',krg,kro)
    if s[ii] > Sfo:
        temp = 1.0 - Sfw - Sfo
        kfo = np.minimum(((s[ii] - Sfo)/temp)**a1f, 1.0)
    kr.append([krw, krg, kro])
    kf.append([kfw, kfg, kfo])
    #
    #   capillary pressure
    #
    cpow = np.nan
    cpog = np.nan
    if s[ii] > Srw:
        cpow = -pwo(s[ii], 1-s[ii], Srw, Sro, b2, b4)
        if cpow < 100:
            cpow = np.nan
    if s[ii] > Sro:
        cpog = -pgo(Srw, s[ii], Srw, Sro, b2, b3)
        if cpog < 100:
            cpog = np.nan
    cp.append([cpow, cpog])

kr = np.array(kr)
kf = np.array(kf)
cp = np.array(cp)

# make plot
cm = 1.0 / 2.54
font = {'family' : 'Arial',
        'size'   : fs}
matplotlib.rc('font', **font)

#
#   Plot relative permeability
#

fig = plt.figure(1, figsize=[18*cm,6*cm])

plt.subplot(1,4,1)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 - 0.05, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(s, kr[:,0], color='b', linewidth=1.0)
plt.plot(1-s, kr[:,2], color='r', linewidth=1.0)
plt.xlabel('Sw',fontsize=fs)
plt.ylabel('Relative permeability',fontsize=fs)
plt.legend(['$k_{rw}$','$k_{rg}$'], fontsize=fs)
# plt.annotate('(a)', (0.5,0.9*12), color='k', fontsize=fs)
plt.title('(a) Matrix and SRV',fontsize=fs)


plt.subplot(1,4,2)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 - 0.02, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(1-s, kr[:,2], color='k', linewidth=1.0)
plt.plot(s, kr[:,1], color='r', linewidth=1.0)
plt.xlabel('Sg',fontsize=fs)
plt.legend(['$k_{ro}$'], loc='lower center',fontsize=fs)
plt.title('(b) Matrix and SRV',fontsize=fs)


plt.subplot(1,4,3)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 + 0.01, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(s, kf[:,0], color='b', linewidth=1.0)
plt.plot(1-s, kf[:,2], color='r', linewidth=1.0)
# plt.plot(s, cp[:,0]/1000.0, color='b', linewidth=1.0)
plt.xlabel('Sw',fontsize=fs)
plt.title('(c) Hydraulic fracture',fontsize=fs)

plt.subplot(1,4,4)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 + 0.05, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(1-s, kf[:,2], color='k', linewidth=1.0)
plt.plot(s, kf[:,1], color='r', linewidth=1.0)
# plt.plot(1-s, cp[:,1]/1000.0, color='r', linewidth=1.0)
plt.xlabel('Sg',fontsize=fs)
plt.title('(d) Hydraulic fracture',fontsize=fs)

# plt.savefig('Figure_6.eps',format='eps')



#
#   Plot capillary pressure
#

fig = plt.figure(2, figsize=[9*cm,6*cm])
plt.subplot(1,2,1)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 + 0.02, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(s, cp[:,0]/1000.0, color='b', linewidth=1.0)
plt.xlabel('Sw',fontsize=fs)
plt.ylabel('Capillary pressure [kPa]',fontsize=fs)
plt.title('(a) $P_{cow}$',fontsize=fs)


plt.subplot(1,2,2)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0 + 0.05, pos1.y0 + 0.04, pos1.width, pos1.height]
ax.set_position(pos2)
plt.plot(1-s, cp[:,1]/1000.0, color='r', linewidth=1.0)
plt.xlabel('Sg',fontsize=fs)
plt.title('(b) $P_{cog}$',fontsize=fs)

# plt.savefig('Figure_7.eps',format='eps')


plt.show()
