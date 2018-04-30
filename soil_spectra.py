from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

import sys
sys.path.append('../general')
from quick_xmu import * 

#Experiment specific stuff
datafile = '1238_Luke'

scans = {}
scans['I0'] = (range(22,25)+[49], [20,21])


scans['Fe'] = (range(14,17), [17,19])
scans['FeS'] = (range(42,46), [47,48])
scans['Fe2O3'] = (range(35,40), [40,41])

#scans['Soil#1'] = (range(28, 33), [26,27])
#scans['Soil#2'] = (range(51, 56), [57,58])
#scans['Soil#3'] = (range(59, 66), [67,68])


#scans['245-16 Py13 #1'] = (range(81,124), [125,126])
#scans['245-14 Py10 #1'] = (range(127,136), [136,137])


scans['I0_s'] = ([168,175,182,189,196,203,210,217,224,231,240], [160,161,232,233])
#scans['I0_s'] = ([168,182,196,210,224,240], [160,161,232,233])
scans['245-16 Py13 #1'] = (range(164,168),(162,163))
scans['245-16 Py13 #2'] = (range(171,175),(169,170))
scans['245-14 Py10 #1'] = (range(178,182),(176,177))
scans['245-14 Py10 #2'] = (range(185,189),(183,184))
scans['101-18 Ku11 #1'] = (range(192,196),(190,191))
scans['101-18 Ku11 #2'] = (range(199,203),(197,198))
scans['101-41 Py07 #1'] = (range(206,210),(204,205))
scans['101-41 Py07 #2'] = (range(213,217),(211,212))
scans['Ym 5-10 cm #1'] = (range(220,224),(218,219))
scans['Ym 5-10 cm #2'] = (range(227,231),(225,226))

scans['"FeCl2"'] = (range(236,240),(234,235))
scans['Fe3Co2OHx'] = (range(243,246),(241,242))

hkl = [5,3,1]

#Helpful functions

def xmu(sample_str,plot=False,scans=scans):
    return quick_xmu(scans['I0'][0],scans[sample_str][0],scans['I0'][1],scans[sample_str][1],plot,datafile)

def absorption(sample_str,plot=False,remove_bg = False,bg_ind = [0],normalize=False,scans=scans):
    th_s, xmu_s = xmu(sample_str,False)
    th_corrected = th_s + th_calibration
    e_corrected = energy(th_corrected, hkl)

    xmu_corrected = xmu_s
    if remove_bg:
        p = np.polyfit(e_corrected[bg_ind],xmu_s[bg_ind],1)
        xmu_corrected = xmu_s - np.polyval(p,e_corrected)
        #xmu_corrected = xmu_s - np.min(xmu_s)

    if normalize:
        #p = np.polyfit(e_corrected[3*np.max(bg_ind):],xmu_corrected[3*np.max(bg_ind):],1)       
        #xmu_corrected = xmu_corrected/np.polyval(p,e_corrected)
        #xmu_corrected = xmu_corrected/np.max(xmu_corrected)
        xmu_corrected = xmu_corrected/np.mean(xmu_corrected[:])

    if plot:
        plt.plot(e_corrected, xmu_corrected,label=sample_str)
    return e_corrected, xmu_corrected


#Theta calibration

'''
th_foil, xmu_foil = xmu('Fe')

#1st derivative max
dth = 0.5*(th_foil[:-1]+th_foil[1:])
dxmu = -np.diff(xmu_foil)/np.diff(th_foil)

p = np.polyfit(dth[94:99], dxmu[94:99],2)

plt.figure()
plt.plot(xmu_foil)
plt.plot(dxmu)

plt.figure()
plt.plot(dth,dxmu)
plt.plot(dth[94:99],np.polyval(p,dth[94:99]))

dmax_theta = -p[1]/(2*p[0])
nominal_theta = theta(7112,hkl)

th_calibration = nominal_theta-dmax_theta

print th_calibration
'''
th_calibration = 0.162700072135

#e,mu = absorption('Co',False,True,range(10))
#A=np.loadtxt('/home/xasadmin/references/Co.dat')
#plt.plot(e,mu/np.max(mu))
#plt.plot(A[:,0],(A[:,1]-np.min(A[:,1]))/(np.max(A[:,1])-np.min(A[:,1])))

e={}
mu = {}

def savespectrum(sample_str):
    ener,xmu = absorption(sample_str,False,False,range(100),False)
    A = np.array([ener,xmu]).T
    np.savetxt(sample_str+'.dat',A)


#e['Soil#1'],mu['Soil#1'] = absorption('Soil#1',True,True,range(50),True)
#e['Soil#2'],mu['Soil#2'] = absorption('Soil#2', True,True,range(50),True)
#e['Soil#3'],mu['Soil#3'] = absorption('Soil#3', True,True,range(50),True)

#e['245-16 Py13 #1'],mu['245-16 Py13 #1'] = absorption('245-16 Py13 #1',True,True,range(50),True)
#e['245-14 Py10 #1'],mu['245-14 Py10 #1'] = absorption('245-14 Py10 #1',True,True,range(50),True)

nbgr = True

#e['Fe'],mu['Fe'] = absorption('Fe',True,nbgr,range(50),nbgr)
#e['FeS'],mu['FeS'] = absorption('FeS',True,nbgr,range(50),nbgr)
#e['Fe2O3'],mu['Fe2O3'] = absorption('Fe2O3',True,nbgr,range(50),nbgr)

scans['I0']=scans['I0_s']
e['245-16 Py13 #1'],mu['245-16 Py13 #1'] = absorption('245-16 Py13 #1', True,nbgr,range(50),nbgr)
e['245-16 Py13 #2'],mu['245-16 Py13 #2'] = absorption('245-16 Py13 #2', True,nbgr,range(50),nbgr)
e['245-14 Py10 #1'],mu['245-14 Py10 #1'] = absorption('245-14 Py10 #1', True,nbgr,range(50),nbgr)
e['245-14 Py10 #2'],mu['245-14 Py10 #2'] = absorption('245-14 Py10 #2', True,nbgr,range(50),nbgr)
e['101-18 Ku11 #1'],mu['101-18 Ku11 #1'] = absorption('101-18 Ku11 #1', True,nbgr,range(50),nbgr)
e['101-18 Ku11 #2'],mu['101-18 Ku11 #2'] = absorption('101-18 Ku11 #2', True,nbgr,range(50),nbgr)
e['101-41 Py07 #1'],mu['101-41 Py07 #1'] = absorption('101-41 Py07 #1', True,nbgr,range(50),nbgr)
e['101-41 Py07 #2'],mu['101-41 Py07 #2'] = absorption('101-41 Py07 #2', True,nbgr,range(50),nbgr)
e['Ym 5-10 cm #1'],mu['Ym 5-10 cm #1'] = absorption('Ym 5-10 cm #1', True,nbgr,range(50),nbgr)
e['Ym 5-10 cm #2'],mu['Ym 5-10 cm #2'] = absorption('Ym 5-10 cm #2', True,nbgr,range(50),nbgr)

e['"FeCl2"'],mu['"FeCl2"'] = absorption('"FeCl2"', True,nbgr,range(50),nbgr)
e['Fe3Co2OHx'],mu['Fe3Co2OHx'] = absorption('Fe3Co2OHx', True,nbgr,range(50),nbgr)


plt.xlabel('Energy (eV)')
plt.ylabel('Absorption (arb. units)')
plt.legend()


#savespectrum('I0')


#A=np.array([e['Se(0)'],mu['Se(0)'],mu['Se(+4)'],mu['Sample#1'],mu['Sample#2'],mu['Null']]).T

#np.savetxt('Se_spectra_first_test.dat',A,header='#Energy Se(0) Se(+2) Sample#1 Sample#2 Starch')



plt.show()
