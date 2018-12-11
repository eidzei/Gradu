from __future__ import division, print_function

from helxas import HelXAS

from bayesfit.fit import least_squares, outlier_fit

import numpy as np
import matplotlib.pyplot as plt

helxas = HelXAS('1251_Luke','/home/xasadmin/data')
helxas.set_analyser('Si',531)
helxas.theta_calibration = 0

helxas_ref = HelXAS('1238_Luke','/home/xasadmin/data')
helxas_ref.set_analyser('Si',531)
helxas.theta_calibration = 0

SI = (range(22,25)+[49],[20,21])

SFe = (range(14,17), [17,19])
SFeS = (range(42,46), [47,48])
SFe2O3 = (range(35,40),[40,41])

I0_ini = ([196],[160,161])
I_m3_ini = (range(192,196),[190,191])


I0_ind = ([167,218],[221,222])
m1_ind = ([155,158,161,164],[153,154])
metal1_ind = ([206,209,212,215],[204,205])
I0s_ind = ([275,326],[329,330])
m1s_ind = ([263,266,269,272],[261,262])
metal1s_ind = ([314,317,320,323],[312,313])

I0_ind2 = ([168,219],[221,222])
m1_ind2 = ([156,159,162,165],[153,154])
metal2_ind = ([207,210,213,216],[204,205])
I0s_ind2 = ([276,327],[329,330])
m1s_ind2 = ([264,267,270,273],[261,262])
metal2s_ind = ([315,318,321,324],[312,313])

I0_ind3 = ([169,220],[221,222])
m1_ind3 = ([157,160,163,166],[153,154])
metal3_ind = ([208,211,214,217],[204,205])
I0s_ind3 = ([277,328],[329,330])
m1s_ind3 = ([265,268,271,274],[261,262])
metal3s_ind = ([316,319,322,325],[312,313])

# Part 1
helxas.read_I0(I0_ind[0],I0_ind[1],20)
helxas.read_I('m11',m1_ind[0],m1_ind[1],20)
helxas.read_I('metal11', metal1_ind[0], metal1_ind[1], 20)

helxas.set_background_fit_order(1,False)

e,mux,mux_err = helxas.get_spectrum('m11')

em, muxm, muxm_err = helxas.get_spectrum('metal11')

helxas.read_I0(I0s_ind[0],I0s_ind[1],20)
helxas.read_I('m11s',m1s_ind[0],m1s_ind[1],20)
helxas.read_I('metal21', metal1s_ind[0], metal1s_ind[1], 20)

helxas.set_background_fit_order(1,False)

es,muxs,muxs_err = helxas.get_spectrum('m11s')
ems, muxms, muxms_err = helxas.get_spectrum('metal21')

np.savetxt('test1.dat',np.array([es,muxs,muxs_err]).T)                               

# Part 2
helxas.read_I0(I0_ind2[0],I0_ind2[1],20)
helxas.read_I('m12',m1_ind2[0],m1_ind2[1],20)
helxas.read_I('metal12', metal2_ind[0], metal2_ind[1], 20)

helxas.set_background_fit_order(1,False)

e2,mux2,mux2_err = helxas.get_spectrum('m12')
em2, muxm2, muxm2_err = helxas.get_spectrum('metal12')

helxas.read_I0(I0s_ind2[0],I0s_ind2[1],20)
helxas.read_I('m12s',m1s_ind2[0],m1s_ind2[1],20)
helxas.read_I('metal22', metal2s_ind[0], metal2s_ind[1], 20)

helxas.set_background_fit_order(1,False)

es2,muxs2,muxs2_err = helxas.get_spectrum('m12s')
ems2, muxms2, muxms2_err = helxas.get_spectrum('metal22')

plt.figure()
#plt.errorbar(e,mux,mux_err)
#plt.errorbar(e2,mux2,mux2_err)

np.savetxt('test2.dat',np.array([es2,muxs2,muxs2_err]).T)

# Part 3
helxas.read_I0(I0_ind3[0],I0_ind3[1],20)
helxas.read_I('m13',m1_ind3[0],m1_ind3[1],20)
helxas.read_I('metal13', metal3_ind[0], metal3_ind[1], 20)

helxas.set_background_fit_order(1,False)

e3,mux3,mux3_err = helxas.get_spectrum('m13')
em3, muxm3, muxm3_err = helxas.get_spectrum('metal13')

helxas.read_I0(I0s_ind3[0],I0s_ind3[1],20)
helxas.read_I('m13s',m1s_ind3[0],m1s_ind3[1],20)
helxas.read_I('metal23', metal3s_ind[0], metal3s_ind[1], 20)

helxas.set_background_fit_order(1,False)

es3,muxs3,muxs3_err = helxas.get_spectrum('m13s')
ems3, muxms3, muxms3_err = helxas.get_spectrum('metal23')

np.savetxt('test3.dat',np.array([es3,muxs3,muxs3_err]).T)
#plt.errorbar(e3,mux3,mux3_err)

# Stitch I

mux_s1 = np.concatenate((mux[e < 7.100423329842326893e+03],mux2[e2 >= 7.100423329842326893e+03]))
mux_s1_err = np.concatenate((mux_err[e < 7.100423329842326893e+03],mux2_err[e2 >= 7.100423329842326893e+03]))
e_s1 = np.concatenate((e[e < 7.100423329842326893e+03],e2[e2 >= 7.100423329842326893e+03]))

mux_s1s = np.concatenate((muxs[es < 7.100423329842326893e+03],muxs2[es2 >= 7.100423329842326893e+03]))
mux_s1s_err = np.concatenate((muxs_err[es < 7.100423329842326893e+03],muxs2_err[es2 >= 7.100423329842326893e+03]))
e_s1s = np.concatenate((es[es < 7.100423329842326893e+03],es2[es2 >= 7.100423329842326893e+03]))

mux_s2 = np.concatenate((mux_s1[e_s1 <= 7.177215327513202283e+03],mux3[e3 > 7.177215327513202283e+03]))
mux_s2_err = np.concatenate((mux_s1_err[e_s1 <= 7.177215327513202283e+03],mux3_err[e3 > 7.177215327513202283e+03]))
e_s2 = np.concatenate((e_s1[e_s1 <= 7.177215327513202283e+03],e3[e3 > 7.177215327513202283e+03]))

mux_s2s = np.concatenate((mux_s1s[e_s1s <= 7.177215327513202283e+03],muxs3[es3 > 7.177215327513202283e+03]))
mux_s2s_err = np.concatenate((mux_s1s_err[e_s1s <= 7.177215327513202283e+03],muxs3_err[es3 > 7.177215327513202283e+03]))
e_s2s = np.concatenate((e_s1s[e_s1s <= 7.177215327513202283e+03],es3[es3 > 7.177215327513202283e+03]))

# Stitch metal foil
mux_m1 = np.concatenate((muxm[em < 7.100423329842326893e+03],muxm2[em2 >= 7.100423329842326893e+03]))
mux_m1_err = np.concatenate((muxm_err[em < 7.100423329842326893e+03],muxm2_err[em2 >= 7.100423329842326893e+03]))
e_m1 = np.concatenate((em[em < 7.100423329842326893e+03],em2[em2 >= 7.100423329842326893e+03]))

mux_m1s = np.concatenate((muxms[ems < 7.100423329842326893e+03],muxms2[ems2 >= 7.100423329842326893e+03]))
mux_m1s_err = np.concatenate((muxms_err[ems < 7.100423329842326893e+03],muxms2_err[ems2 >= 7.100423329842326893e+03]))
e_m1s = np.concatenate((ems[ems < 7.100423329842326893e+03],ems2[ems2 >= 7.100423329842326893e+03]))

mux_m2 = np.concatenate((mux_m1[e_m1 <= 7.177215327513202283e+03],muxm3[em3 > 7.177215327513202283e+03]))
mux_m2_err = np.concatenate((mux_m1_err[e_m1 <= 7.177215327513202283e+03],muxm3_err[em3 > 7.177215327513202283e+03]))
e_m2 = np.concatenate((e_s1[e_s1 <= 7.177215327513202283e+03],e3[e3 > 7.177215327513202283e+03]))

mux_m2s = np.concatenate((mux_m1s[e_m1s <= 7.177215327513202283e+03],muxms3[ems3 > 7.177215327513202283e+03]))
mux_m2s_err = np.concatenate((mux_m1s_err[e_m1s <= 7.177215327513202283e+03],muxms3_err[ems3 > 7.177215327513202283e+03]))
e_m2s = np.concatenate((e_m1s[e_m1s <= 7.177215327513202283e+03],ems3[ems3 > 7.177215327513202283e+03]))

#Reference samples
helxas_ref.read_I0(SI[0],SI[1],8)
helxas_ref.read_I('Fe',SFe[0],SFe[1],8)
helxas_ref.read_I('FeS', SFeS[0],SFeS[1],8)
helxas_ref.read_I('Fe2O3', SFe2O3[0],SFe2O3[1],8)

helxas_ref.set_background_fit_order(3,False)
e_Fe, mux_Fe, mux_err_Fe = helxas_ref.get_spectrum('Fe')
e_FeS, mux_FeS, mux_err_FeS = helxas_ref.get_spectrum('FeS')
e_Fe2O3, mux_Fe2O3, mux_err_Fe2O3 = helxas_ref.get_spectrum('Fe2O3')


helxas_ref.read_I0(I0_ini[0],I0_ini[1],40)
helxas_ref.read_I('M3_dry',I_m3_ini[0],I_m3_ini[1],40)

helxas_ref.set_background_fit_order(3,False)
e_dry, mux_dry, mux_err_dry = helxas_ref.get_spectrum('M3_dry')

e_Fe = e_Fe - 6.5
e_FeS  = e_FeS - 8.4
e_dry = e_dry - 8.4
e_Fe2O3 = e_Fe2O3 -8.4

#Nollatason leikkaus

n_s2 = np.polyfit(e_s2[:86],mux_s2[:86],1)
n_s2s = np.polyfit(e_s2s[:86],mux_s2s[:86],1)
n_m2 = np.polyfit(em,muxm,1)
n_m2s = np.polyfit(ems,muxms,1)
fit_e = np.linspace(7.061455915797724629e+03,7.233339179611556574e+03,301)
mux_FeS = np.interp(fit_e, e_FeS, mux_FeS)
mux_Fe2O3 = np.interp(fit_e, e_Fe2O3, mux_Fe2O3)
mux_dry = np.interp(fit_e, e_dry, mux_dry)
n_FeS = np.polyfit(fit_e[:86],mux_FeS[:86],1)
n_Fe2O3 = np.polyfit(fit_e[:86],mux_Fe2O3[:86],1)
n_dry = np.polyfit(fit_e[:86], mux_dry[:86],1)


mux_s2 = mux_s2 - np.polyval(n_s2,fit_e)-0.005
mux_s2s = mux_s2s - np.polyval(n_s2s,fit_e)
mux_m2 = mux_m2 - np.polyval(n_m2,fit_e)
mux_m2s = mux_m2s - np.polyval(n_m2s,fit_e)
mux_Fe = mux_Fe - np.min(mux_Fe)
mux_FeS = mux_FeS - np.polyval(n_FeS,fit_e)
mux_Fe2O3 = mux_Fe2O3 - np.polyval(n_Fe2O3,fit_e)
mux_dry  = mux_dry -np.polyval(n_dry,fit_e)

#Normalisointi
mux_s2s = mux_s2s*(1+0.0015*(e_s2s-7220))
mux_s2s_err = mux_s2s_err/np.mean(mux_s2s[np.logical_and(e_s2s<7230,e_s2s>7220)])


mux_s2 = mux_s2/np.mean(mux_s2[np.logical_and(e_s2<7230,e_s2>7220)])
mux_s2s = mux_s2s/np.mean(mux_s2s[np.logical_and(e_s2s<7230,e_s2s>7220)])
mux_FeS = mux_FeS/np.mean(mux_FeS[np.logical_and(fit_e<7230,fit_e>7220)])
mux_Fe2O3 = mux_Fe2O3/np.mean(mux_Fe2O3[np.logical_and(fit_e<7230,fit_e>7220)])
mux_dry = mux_dry/np.mean(mux_dry[np.logical_and(fit_e<7230,fit_e>7220)])

#mux_s2 = mux_s2/np.trapz(mux_s2,dx=1)
#mux_s2s = mux_s2s/np.trapz(mux_s2s,dx=1)
#mux_FeS = mux_FeS/np.trapz(mux_FeS,dx=1)

mux_s2s = mux_s2s*(1+0.0018*(e_s2s-7175))
mux_FeS = mux_FeS*(1-0.001*(e_s2s-7200))

#Prosentin selvitys

ind=np.logical_and(fit_e < 7170, fit_e > 7100)


def model(x, A, B):
    y_Fe2O3 = np.interp(x, fit_e, mux_Fe2O3)
    y_dry = np.interp(x, fit_e, mux_dry)

    return A*y_Fe2O3 + B*y_dry 

#fit = outlier_fit(model, [0.2,0.8], e_s2s[ind],mux_s2s[ind], mux_s2s_err[ind], method='cauchy')
#fit = outlier_fit(model, [0.2,0.8], e_s2s[ind],mux_s2s[ind], np.min(mux_s2s_err),method='conservative')
fit = least_squares(model, [0.2,0.8], e_s2s,mux_s2s, mux_s2s_err)
res_p = fit.get_result()[0]
res_cov = fit.get_result()[1]


#Plotting
np.savetxt('stich_test.dat',np.array([e_s2,mux_s2,mux_s2_err]).T)
plt.plot(e_s2,mux_s2, label='M3 C0 S1 a')
plt.plot(fit_e,mux_dry, label='M3 dry') 
#plt.plot(e_s2s,mux_s2s, label='M3 C1 S1 a')
#plt.errorbar(e_m2,mux_m2/np.max(mux_m2), label='metal') 
#plt.errorbar(e_m2s,mux_m2s/np.max(mux_m2s), label='metal s')
#plt.errorbar(e_Fe,mux_Fe/np.max(mux_Fe), label='ref metal')
#plt.plot(fit_e,mux_FeS, label='ref FeS')
plt.plot(fit_e,mux_Fe2O3, label='Fe2O3')
#plt.plot(fit_e, model(fit_e,fit.p[0],fit.p[1]), label='fit FeS = %.4f $\pm$ %.4f, M3 = %.4f $\pm$ %.4f' %(res_p[0], res_cov[0], res_p[1], res_cov[1]))
plt.legend()
plt.show()

