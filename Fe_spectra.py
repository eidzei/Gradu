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

scans['245-16 Py13 #1'] = (range(81,124), [125,126])
scans['245-14 Py10 #1'] = (range(127,136), [136,137])

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
        xmu_corrected = xmu_corrected/np.mean(xmu_corrected[150:])

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

e['Fe'],mu['Fe'] = absorption('Fe',True,True,range(50),True)
e['FeS'],mu['FeS'] = absorption('FeS',True,True,range(50),True)
e['Fe2O3'],mu['Fe2O3'] = absorption('Fe2O3',True,True,range(50),True)

#e['Soil#1'],mu['Soil#1'] = absorption('Soil#1',True,True,range(50),True)
#e['Soil#2'],mu['Soil#2'] = absorption('Soil#2', True,True,range(50),True)
#e['Soil#3'],mu['Soil#3'] = absorption('Soil#3', True,True,range(50),True)

e['245-16 Py13 #1'],mu['245-16 Py13 #1'] = absorption('245-16 Py13 #1',True,True,range(50),True)
e['245-14 Py10 #1'],mu['245-14 Py10 #1'] = absorption('245-14 Py10 #1',True,True,range(50),True)

plt.xlabel('Energy (eV)')
plt.ylabel('Absorption (arb. units)')
plt.legend()


#savespectrum('I0')


#A=np.array([e['Se(0)'],mu['Se(0)'],mu['Se(+4)'],mu['Sample#1'],mu['Sample#2'],mu['Null']]).T

#np.savetxt('Se_spectra_first_test.dat',A,header='#Energy Se(0) Se(+2) Sample#1 Sample#2 Starch')



plt.show()
