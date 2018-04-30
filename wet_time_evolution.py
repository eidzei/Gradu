from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from pyspec.spec import SpecDataFile

datapath = '/home/xasadmin/data/'

datafile = '1238_Luke'

sd = SpecDataFile(datapath + datafile)

data = []

for i in range(369, 475):
   data.append(sd[i].Detector)

data = np.array(data).T

time = []
for i in range(0,10):
   time.append(i*2*1.0)
   
time = np.array(time)
   

print time
  
    
ax = plt.subplot()

ax.plot(data.T)               
ax.xaxis.set_ticks(np.arange(0, 110, 2*5.474325))
ax.set_xticklabels(time)
plt.xlabel('time(h)')
#plt.ylabel('absortion (arb. units)')
#plt.legend()

plt.show()
