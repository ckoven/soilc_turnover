## script to generate the graphics in figure 3 of koven et al. paper
## written by c. koven
## dependent on several libraries: PyNgl, python-netCDF4, numpy, rpy2, and my personal plotting library, which is here: https://github.com/ckoven/ckplotlib

import numpy as np
import map_funcs
import netCDF4 as nc
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
stats = importr('stats')
base = importr('base')


def q10func(temp):
    k_ref = 0.06
    q10 = 1.5
    tref = 15.
    output = k_ref * q10 ** ((temp[:] - tref)/10.)
    return output

def q10func_frozen(temp):
    k_ref = 0.06
    q10 = 1.5
    tref = 15.
    output =  k_ref * q10 ** ((temp[:] - tref)/10.)
    output[temp < 0.] = 0.
    return output

def q10func_frozenq10(temp):
    k_ref = 0.06
    q10 = 1.5
    q10_frozen=5.
    tref = 15.
    output =  k_ref * q10 ** ((temp[:] - tref)/10.)
    output[temp < 0.] =  k_ref * (q10 ** ((0. - tref)/10.)) * q10_frozen ** ((temp[temp < 0.] - 0)/10.)
    return output

def lloyd_taylor_eq11(temp):
    k_ref = .03820684  ### at 10C
    output = k_ref * np.exp(308.56 * (1./56.02 - 1/(temp + 273.15 - 227.13)))
    return output

coloreddotfig_dotsize = 0.005

textfile = open('quadratic_parameters_strawmodels.txt', 'w')


# ### the quadratic fit to filtered observational data from other script
# quadratic_fit = np.array([  7.70383464e-04,  -4.00038392e-02,   1.57171229e+00])
# xd = np.linspace(-25, 29, 60)
# yd_quad2_kspace = 1./(10. ** (quadratic_fit[2] + quadratic_fit[1]*xd + quadratic_fit[0]*xd**2 ))

### read quadratic fit and prediction intervals from other script
x_confint = np.loadtxt('x_confint.txt')
y_confint = np.loadtxt('y_confint.txt')
y_confint_q10 = np.row_stack([y_confint,1./q10func(x_confint)])
y_confint_froz_q10 = np.ma.masked_array(np.row_stack([y_confint,1./q10func(x_confint)]))
y_confint_froz_q10[3,x_confint[:] < 0.] = np.ma.masked_all((x_confint[:] < 0.).sum())
y_confint_lt = np.row_stack([y_confint,1./lloyd_taylor_eq11(x_confint)])



f = nc.MFDataset("datasets/clm5_cdk_r162_2degGSWP3_1850spin_allN_modpsimin_v01.clm2.h1.0*.nc")
tsa = f.variables['TSA']
tsoi = f.variables['TSOI']
lat = f.variables['lat']
lon = f.variables['lon']
levgrnd = f.variables['levgrnd']

f2 = nc.Dataset("datasets/clm5_cdk_r162_2degGSWP3_1850spin_allN_modpsimin_v01__dzsoi.nc")
dzsoi = f2.variables['DZSOI'][:]
f2.close()


f3 = nc.Dataset("datasets/surfdata_1.9x2.5_16pftsmidarctic_simyr1850_c160112.nc")
icefrac = f3.variables["PCT_GLACIER"][:] / 100.
f3.close()

dz = dzsoi[:,0,0]
lev_edges = np.zeros(len(dz)+1)
lev_edges[1:] = dz.cumsum()


depth_1m = 1.
nlev_lt_1m = (lev_edges < depth_1m).sum()
dz_to1m = dz[0:nlev_lt_1m].data
dz_to1m[nlev_lt_1m-1] = dz_to1m[nlev_lt_1m-1] - (lev_edges[nlev_lt_1m] - depth_1m)

frozen = (tsoi[:,0:nlev_lt_1m, :,:] < 273.15).mean(axis=0)

frozen_1m = (np.rollaxis(np.rollaxis(frozen,1), 2, 1) * dz_to1m).sum(axis=2)

tsa_mean  = tsa[:].mean(axis=0)

temprange = [-22,30]


index_10cm = 2
q10func_10cm = q10func(tsoi[:,index_10cm,:,:]-273.15).mean(axis=0)
froz_q10func_10cm = q10func_frozen(tsoi[:,index_10cm,:,:]-273.15).mean(axis=0)

froz_q10func_0to1mm = (np.rollaxis(np.rollaxis(q10func_frozen(tsoi[:,0:nlev_lt_1m,:,:]-273.15).mean(axis=0),1), 2, 1) * dz_to1m).sum(axis=2)

lloydtayl_10cm = lloyd_taylor_eq11(tsoi[:,index_10cm,:,:]-273.15).mean(axis=0)


### mask values below certain level
mintemp = 273.15 -30.
icefracmax = 0.2
glaciercoldmask = np.logical_and(tsa_mean > mintemp, icefrac < icefracmax)
lloydtayl_10cm = lloydtayl_10cm[glaciercoldmask]
q10func_10cm = q10func_10cm[glaciercoldmask]
froz_q10func_10cm = froz_q10func_10cm[glaciercoldmask]
froz_q10func_0to1mm = froz_q10func_0to1mm[glaciercoldmask]
tsa_mean = tsa_mean[glaciercoldmask]

yrange=[1e0, 5e3]


map_funcs.xyplot(tsa_mean-273.15, 1./q10func_10cm,  ylog=True, yrange=yrange, xrange=temprange, dots=True, overlay_x = x_confint, overlay_y = y_confint_q10,overlay_color=['red','red','red','blue'],overlay_linethickness=[2.5,1.,1.,2.5], xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Turnover time of respiration function (yr)', inset_title='Q~B~10~N~=1.5 at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3, file='diagnosed_straight_q10_10cm')

map_funcs.xyplot(tsa_mean-273.15, 1./np.maximum(froz_q10func_10cm, 1e-4),  ylog=True, yrange=yrange, xrange=temprange, dots=True, overlay_x = x_confint, overlay_y = y_confint_froz_q10,overlay_color=['red','red','red','blue'],overlay_linethickness=[2.5,1.,1.,2.5], xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Turnover time of respiration function (yr)', inset_title='Thawed-only Q~B~10~N~=1.5 at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3, file='diagnosed_frozen_q10_10cm')

map_funcs.xyplot(tsa_mean-273.15, 1./np.maximum(lloydtayl_10cm, 1e-4),  ylog=True, yrange=yrange, xrange=temprange, dots=True, overlay_x = x_confint, overlay_y = y_confint_lt,overlay_color=['red','red','red','blue'],overlay_linethickness=[2.5,1.,1.,2.5], xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Turnover time of respiration function (yr)', inset_title='Lloyd-Taylor at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3, file='diagnosed_lloyd-taylor_10cm')

map_funcs.xyplot(tsa_mean-273.15, 1./np.maximum(froz_q10func_0to1mm, 1e-4),  ylog=True, yrange=yrange, xrange=temprange, dots=True, overlay_x = x_confint, overlay_y = y_confint_froz_q10,overlay_color=['red','red','red','blue'],overlay_linethickness=[2.5,1.,1.,2.5], xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Turnover time of respiration function (yr)', inset_title='Thawed-only Q~B~10~N~=1.5 over 0-1m interval', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3, file='diagnosed_frozen_q10_0to1m_integral')

### for each of these, calculate the quadratic regression as for the obs and ESMs


### first the simple q10 case
xdata_nomask = tsa_mean-273.15
ydata_nomask = np.log10(1./q10func_10cm)

xdata = xdata_nomask[np.logical_not(ydata_nomask.mask)]
ydata = ydata_nomask[np.logical_not(ydata_nomask.mask)]

xdata_rvect_simpleq10 = robjects.FloatVector(xdata)
ydata_rvect_simpleq10 = robjects.FloatVector(ydata)
#
robjects.globalenv["xdata_rvect_simpleq10"] = xdata_rvect_simpleq10
robjects.globalenv["ydata_rvect_simpleq10"] = ydata_rvect_simpleq10
#
quadreg_r_simpleq10 = stats.lm("ydata_rvect_simpleq10 ~ poly(xdata_rvect_simpleq10,2)")
robjects.globalenv["quadreg_r_simpleq10"] = quadreg_r_simpleq10
rconfint_simpleq10 = robjects.r['confint']
rpredict_simpleq10 = robjects.r['predict']
rsummary_simpleq10 = robjects.r['summary']
thepredictint_simpleq10 = rpredict_simpleq10(quadreg_r_simpleq10, interval='prediction', level=0.50)
print(rsummary_simpleq10(quadreg_r_simpleq10))
predictint_simpleq10 = np.array(thepredictint_simpleq10)
n_toshow_predictlines_simpleq10 = 50
n_toskip_predictlines_simpleq10 = xdata.shape[0]/n_toshow_predictlines_simpleq10
indices_toshow_simpleq10 = xdata.argsort()[::n_toskip_predictlines_simpleq10]
x_confint_simpleq10 = xdata[indices_toshow_simpleq10]
y_confint_simpleq10 = 10.**(predictint_simpleq10[indices_toshow_simpleq10,:].transpose())
#
# do in numpy for coefficients
quadratic_fit_simpleq10 = np.polyfit(xdata,ydata,2)
textfile.write('quadratic_fit_simpleq10 '+str(quadratic_fit_simpleq10)+'\n')
#
map_funcs.xyplot(tsa_mean-273.15, 1./q10func_10cm, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='simpleq10_MRT_soilc_temp_quadraticregression_r50pctpredint_', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_simpleq10, overlay_y = y_confint_simpleq10,overlay_color='red',overlay_linethickness=[2.5,1.,1.], inset_title='Q~B~10~N~=1.5 at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3 )


### next the lloyd-taylor
xdata_nomask = tsa_mean-273.15
ydata_nomask = np.log10(1./lloydtayl_10cm)

xdata = xdata_nomask[np.logical_not(ydata_nomask.mask)]
ydata = ydata_nomask[np.logical_not(ydata_nomask.mask)]

xdata_rvect_lloydtaylor = robjects.FloatVector(xdata)
ydata_rvect_lloydtaylor = robjects.FloatVector(ydata)
#
robjects.globalenv["xdata_rvect_lloydtaylor"] = xdata_rvect_lloydtaylor
robjects.globalenv["ydata_rvect_lloydtaylor"] = ydata_rvect_lloydtaylor
#
quadreg_r_lloydtaylor = stats.lm("ydata_rvect_lloydtaylor ~ poly(xdata_rvect_lloydtaylor,2)")
robjects.globalenv["quadreg_r_lloydtaylor"] = quadreg_r_lloydtaylor
rconfint_lloydtaylor = robjects.r['confint']
rpredict_lloydtaylor = robjects.r['predict']
rsummary_lloydtaylor = robjects.r['summary']
thepredictint_lloydtaylor = rpredict_lloydtaylor(quadreg_r_lloydtaylor, interval='prediction', level=0.50)
print(rsummary_lloydtaylor(quadreg_r_lloydtaylor))
predictint_lloydtaylor = np.array(thepredictint_lloydtaylor)
n_toshow_predictlines_lloydtaylor = 50
n_toskip_predictlines_lloydtaylor = xdata.shape[0]/n_toshow_predictlines_lloydtaylor
indices_toshow_lloydtaylor = xdata.argsort()[::n_toskip_predictlines_lloydtaylor]
x_confint_lloydtaylor = xdata[indices_toshow_lloydtaylor]
y_confint_lloydtaylor = 10.**(predictint_lloydtaylor[indices_toshow_lloydtaylor,:].transpose())
#
# do in numpy for coefficients
quadratic_fit_lloydtaylor = np.polyfit(xdata,ydata,2)
textfile.write('quadratic_fit_lloydtaylor '+str(quadratic_fit_lloydtaylor)+'\n')

#
#
map_funcs.xyplot(tsa_mean-273.15, 1./lloydtayl_10cm, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='lloydtaylor_MRT_soilc_temp_quadraticregression_r50pctpredint_', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_lloydtaylor, overlay_y = y_confint_lloydtaylor,overlay_color='red',overlay_linethickness=[2.5,1.,1.], inset_title='Lloyd-Taylor at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3,  )


### next frozen q10 case
xdata_nomask = tsa_mean-273.15
ydata_nomask = np.log10(1./froz_q10func_10cm)

xdata = xdata_nomask[np.logical_not(ydata_nomask.mask)]
ydata = ydata_nomask[np.logical_not(ydata_nomask.mask)]

xdata_rvect_froz_q10func_10cm = robjects.FloatVector(xdata)
ydata_rvect_froz_q10func_10cm = robjects.FloatVector(ydata)
#
robjects.globalenv["xdata_rvect_froz_q10func_10cm"] = xdata_rvect_froz_q10func_10cm
robjects.globalenv["ydata_rvect_froz_q10func_10cm"] = ydata_rvect_froz_q10func_10cm
#
quadreg_r_froz_q10func_10cm = stats.lm("ydata_rvect_froz_q10func_10cm ~ poly(xdata_rvect_froz_q10func_10cm,2)")
robjects.globalenv["quadreg_r_froz_q10func_10cm"] = quadreg_r_froz_q10func_10cm
rconfint_froz_q10func_10cm = robjects.r['confint']
rpredict_froz_q10func_10cm = robjects.r['predict']
rsummary_froz_q10func_10cm = robjects.r['summary']
thepredictint_froz_q10func_10cm = rpredict_froz_q10func_10cm(quadreg_r_froz_q10func_10cm, interval='prediction', level=0.50)
print(rsummary_froz_q10func_10cm(quadreg_r_froz_q10func_10cm))
predictint_froz_q10func_10cm = np.array(thepredictint_froz_q10func_10cm)
n_toshow_predictlines_froz_q10func_10cm = 50
n_toskip_predictlines_froz_q10func_10cm = xdata.shape[0]/n_toshow_predictlines_froz_q10func_10cm
indices_toshow_froz_q10func_10cm = xdata.argsort()[::n_toskip_predictlines_froz_q10func_10cm]
x_confint_froz_q10func_10cm = xdata[indices_toshow_froz_q10func_10cm]
y_confint_froz_q10func_10cm = 10.**(predictint_froz_q10func_10cm[indices_toshow_froz_q10func_10cm,:].transpose())
#
# do in numpy for coefficients
quadratic_fit_froz_q10func_10cm = np.polyfit(xdata,ydata,2)
textfile.write('quadratic_fit_froz_q10func_10cm '+str(quadratic_fit_froz_q10func_10cm)+'\n')
#
#
map_funcs.xyplot(tsa_mean-273.15, 1./froz_q10func_10cm, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='froz_q10func_10cm_MRT_soilc_temp_quadraticregression_r50pctpredint_', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_froz_q10func_10cm, overlay_y = y_confint_froz_q10func_10cm,overlay_color='red',overlay_linethickness=[2.5,1.,1.], inset_title='Thawed-only Q~B~10~N~=1.5 at 10cm', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3 )


### and the depth-averaged frozen q10 case
xdata_nomask = tsa_mean-273.15
ydata_nomask = np.log10(1./froz_q10func_0to1mm)

xdata = xdata_nomask[np.logical_not(ydata_nomask.mask)]
ydata = ydata_nomask[np.logical_not(ydata_nomask.mask)]

xdata_rvect_froz_q10func_0to1mm = robjects.FloatVector(xdata)
ydata_rvect_froz_q10func_0to1mm = robjects.FloatVector(ydata)
#
robjects.globalenv["xdata_rvect_froz_q10func_0to1mm"] = xdata_rvect_froz_q10func_0to1mm
robjects.globalenv["ydata_rvect_froz_q10func_0to1mm"] = ydata_rvect_froz_q10func_0to1mm
#
quadreg_r_froz_q10func_0to1mm = stats.lm("ydata_rvect_froz_q10func_0to1mm ~ poly(xdata_rvect_froz_q10func_0to1mm,2)")
robjects.globalenv["quadreg_r_froz_q10func_0to1mm"] = quadreg_r_froz_q10func_0to1mm
rconfint_froz_q10func_0to1mm = robjects.r['confint']
rpredict_froz_q10func_0to1mm = robjects.r['predict']
rsummary_froz_q10func_0to1mm = robjects.r['summary']
thepredictint_froz_q10func_0to1mm = rpredict_froz_q10func_0to1mm(quadreg_r_froz_q10func_0to1mm, interval='prediction', level=0.50)
print(rsummary_froz_q10func_0to1mm(quadreg_r_froz_q10func_0to1mm))
predictint_froz_q10func_0to1mm = np.array(thepredictint_froz_q10func_0to1mm)
n_toshow_predictlines_froz_q10func_0to1mm = 50
n_toskip_predictlines_froz_q10func_0to1mm = xdata.shape[0]/n_toshow_predictlines_froz_q10func_0to1mm
indices_toshow_froz_q10func_0to1mm = xdata.argsort()[::n_toskip_predictlines_froz_q10func_0to1mm]
x_confint_froz_q10func_0to1mm = xdata[indices_toshow_froz_q10func_0to1mm]
y_confint_froz_q10func_0to1mm = 10.**(predictint_froz_q10func_0to1mm[indices_toshow_froz_q10func_0to1mm,:].transpose())
#
# do in numpy for coefficients
quadratic_fit_froz_q10func_0to1mm = np.polyfit(xdata,ydata,2)
textfile.write('quadratic_fit_froz_q10func_0to1mm ' +str(quadratic_fit_froz_q10func_0to1mm)+'\n')#
#
#
map_funcs.xyplot(tsa_mean-273.15, 1./froz_q10func_0to1mm, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='froz_q10func_0to1mm_MRT_soilc_temp_quadraticregression_r50pctpredint_', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_froz_q10func_0to1mm, overlay_y = y_confint_froz_q10func_0to1mm,overlay_color='red',overlay_linethickness=[2.5,1.,1.], inset_title='Thawed-only Q~B~10~N~=1.5 over 0-1m interval', inset_title_x=27.5, inset_textjust="CenterRight", inset_title_y=2.9e3 )

temps = np.arange(-30,30,0.1)
map_funcs.xyplot(temps, 10**(1.54871640e+00 + temps * (-3.68422434e-02) + temps **2 * (5.79263319e-04) ), ylog=True, yrange=[1., 5e3], xrange=temprange, file='sanitycheck')

textfile.close()
