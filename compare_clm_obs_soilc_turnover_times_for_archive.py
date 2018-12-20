## script to generate the graphics in koven et al. paper on the temperatuire sensitivity of soil carbon turnover times
## written by c. koven
## dependent on several libraries: PyNgl, python-netCDF4, numpy, rpy2, scipy, and my personal plotting library, which is here: https://github.com/ckoven/ckplotlib


import Nio
import map_funcs
import numpy as np
import geog_funcs
import linreg
from scipy import optimize
import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
stats = importr('stats')
base = importr('base')
import netCDF4 as nc
from interplib import conservative_regrid

make_data_plots = True
make_cmip5_plots  = True
make_clm5_plots = False
make_clm4_5_plots = True
make_mimics_plots = True
make_multimodel_cmip5_plots = make_cmip5_plots and make_clm5_plots
load_thematic_maps = True

read_clay_map = True

min_pmpet = -1000.
    
hist_yedges = np.arange(0, 4, 0.1)
hist_xedges = np.arange(-30,34., 2.)

log_mrt_levels = 10.**np.arange(.6666666666, 3.33333333, 0.33333333)

rain_colormap = 'GMT_drywet'
rain_colorlevels=np.arange(0,2000., 100)

coloreddotfig_dotsize = 0.005

temprange = [-22,30]

min_T = -15.
max_T = 28.
Tbin_width = 1.
n_T_bins = int((max_T - min_T) / Tbin_width)
T_bin_centers = np.arange(n_T_bins) + min_T + Tbin_width/2.

textfile = open('quadratic_parameters_obs_esms.txt', 'w')

soils_colormap = 'wh-bl-gr-ye-re'
npp_colormap = 'WhiteGreen'
            
if make_data_plots:
    ### first open the observational datasets
    ncscd_datafilename = 'datasets/NCSCDV22_soilc_0.5x0.5.nc'
    hwsd_datafilename = 'datasets/HWSD_soilc_0.5x0.5.nc'
    cru_datafilename = 'datasets/cru_ts_3.1_climo_1961-1990.nc'
    modisnpp_datafilename = 'datasets/MOD17A3_Science_NPP_mean_00_14_regridhalfdegree.nc'
    #
    print( ' opening file '+ncscd_datafilename)
    ncscd_datafile = Nio.open_file(ncscd_datafilename)
    ncscd_data = ncscd_datafile.variables['soilc'][:]
    ncscd_lats = ncscd_datafile.variables['lat'][:]
    ncscd_lons = ncscd_datafile.variables['lon'][:]
    #
    print( ' opening file '+hwsd_datafilename)
    hwsd_datafile = Nio.open_file(hwsd_datafilename)
    hwsd_data = hwsd_datafile.variables['soilc'][:]
    hwsd_lats = hwsd_datafile.variables['lat'][:]
    hwsd_lons = hwsd_datafile.variables['lon'][:]
    #
    print( ' opening file '+cru_datafilename)
    cru_datafile = Nio.open_file(cru_datafilename)
    cru_data = cru_datafile.variables['tmp'][:].mean(axis=0)
    cru_lats = cru_datafile.variables['lat'][:]
    cru_lons = cru_datafile.variables['lon'][:]
    #
    print( ' opening file '+modisnpp_datafilename)
    modisnpp_datafile = Nio.open_file(modisnpp_datafilename)
    modisnpp_data = modisnpp_datafile.variables['npp'][:]
    modisnpp_lats = modisnpp_datafile.variables['lat'][:]
    modisnpp_lons = modisnpp_datafile.variables['lon'][:]
    #
    lats_common = ncscd_lats.copy()
    lons_common = ncscd_lons.copy()
    #
    #
    merged_hwsd_ncscd = hwsd_data.copy()
    merged_hwsd_ncscd[ncscd_data[:] > 0.] = ncscd_data[ncscd_data[:] > 0.] 
    #
    #
    mrt = np.ma.masked_array(merged_hwsd_ncscd/(modisnpp_data * 1e-3), mask= modisnpp_data<0.1)
    numerator = np.ma.masked_array(merged_hwsd_ncscd, mask= modisnpp_data<0.1)
    denominatorator = np.ma.masked_array(modisnpp_data * 1e-3, mask= modisnpp_data<0.1)
    #
    map_funcs.xyplot(cru_data.flatten(), mrt.flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)')
    #
    map_funcs.xyplot(cru_data.flatten(), numerator.flatten(), dots=True, ylog=True, yrange=[.1, 5e2], xrange=temprange, file='obs_soilc_temp', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Soil Carbon to 1m (kg C m~S~-2~N~)')
    #
    map_funcs.xyplot(cru_data.flatten(), denominatorator.flatten(), dots=True, ylog=True, yrange=[1.e-3, 5.], xrange=temprange, file='obs_npp_temp', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='MODIS NPP (kg C m~S~-2~N~ y~S~-1~N~)')
    #
    map_funcs.xyplot(cru_data.flatten(), 1./denominatorator.flatten(), dots=True, ylog=True, yrange=[3.e-1, 1.e3], xrange=temprange, file='obs_oneover_npp_temp', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='1 / MODIS NPP (kg C~S~-1~N~ m~S~2~N~ y~S~~N~)')
    #
    common_mask = np.logical_or(mrt.mask[:], cru_data.mask[:])
    #
    map_funcs.fill(cru_data, lats_common, lons_common, levels=np.arange(-30,35,5), file='map_cru_obs', title='CRU Mean Annual Air Temp (~S~o~N~C)', specialprojection='global_noant')
    #
    cru_data = np.ma.masked_array(cru_data, mask=common_mask)
    mrt = np.ma.masked_array(mrt, mask=common_mask)
    #
    map_funcs.fill(mrt, lats_common, lons_common, levels=log_mrt_levels, file='map_mrt_obs', title='Inferred Turnover Time (years)', specialprojection='global_noant')
    map_funcs.fill(merged_hwsd_ncscd, lats_common, lons_common, levels=np.arange(0,55,5), file='map_soilc_obs', title='HWSD & NCSCD Soil C to 1m (kg m~S~-2~N~)', specialprojection='global_noant', colormap=soils_colormap)
    map_funcs.fill(modisnpp_data, lats_common, lons_common, levels=np.arange(0,1100,100), file='map_npp_obs', title='MODIS NPP (g C m~S~-2~N~ y~S~-1~N~)', specialprojection='global_noant', colormap=npp_colormap)
    #
    #
    ## now plot as a 2d histogram:
    hist_x = cru_data.flatten()
    hist_y = np.log10(mrt.flatten())
    areas = geog_funcs.gridcell_areas(lats_common, lons_common, mask=np.logical_or(cru_data.mask[:], mrt.mask[:])).flatten()
    #
    ### make common mask for data
    commonmask = np.logical_or(hist_x.mask[:], hist_y.mask[:])
    hist_x = hist_x[np.logical_not(commonmask)]
    hist_y = hist_y[np.logical_not(commonmask)]
    areas = areas[np.logical_not(commonmask)]
    #
    #
    common_mask = np.logical_or(cru_data.mask[:], mrt.mask[:])
    common_mask = np.logical_or(common_mask[:], mrt < 1e-1)
    select = np.logical_not(common_mask)
    toskip=2
    #
    map_funcs.xyplot(cru_data[select].flatten()[::toskip].data, mrt[select].flatten()[::toskip].data, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_coloredbsoilC', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', shaded_dot_data=merged_hwsd_ncscd[select].flatten()[::toskip].data, shaded_dot_levels=np.arange(0,55,5), colormap=soils_colormap, shuffle_shaded_dots=True)
    #
    #
    #
    #### now get organic soil maps:
    histel_map = 'datasets/NCSCD_Circumarctic_histel_pct_05deg.nc'
    print(' opening file: '+histel_map)
    histel_file = nc.Dataset(histel_map)
    histel_pct = histel_file.variables['NCSCD_Circumarctic_histel_pct_05deg.tif'][:]
    histel_lats = histel_file.variables['lat'][:]
    histel_lons = histel_file.variables['lon'][:]
    histel_file.close()
    #
    histosol_map = 'datasets/NCSCD_Circumarctic_histosol_pct_05deg.nc'
    print(' opening file: '+histosol_map)
    histosol_file = nc.Dataset(histosol_map)
    histosol_pct = histosol_file.variables['NCSCD_Circumarctic_histosol_pct_05deg.tif'][:]
    histosol_lats = histosol_file.variables['lat'][:]
    histosol_lons = histosol_file.variables['lon'][:]
    histosol_file.close()
    #
    max_hist_frac = 50.
    #
    organics = histel_pct.astype('float') + histosol_pct.astype('float')
    #
    ## map_funcs.fill(organics[::-1,:].astype('float'), histosol_lats[::-1], histosol_lons[:])
    organics_lat_offset = int((histosol_lats.min() - lats_common.min() )*2)
    IM_common = len(lons_common)
    JM_common = len(lats_common)
    organics_commonmap = np.zeros([JM_common, IM_common], dtype=np.bool)
    organics_commonmap[organics_lat_offset:organics_lat_offset+len(histosol_lats),:] = organics[::-1,:] > max_hist_frac
    # map_funcs.fill(organics_commonmap, lats_common, lons_common)
    map_funcs.xyplot(np.ma.masked_array(cru_data, mask=organics_commonmap).flatten(), np.ma.masked_array(mrt, mask=organics_commonmap).flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_peatlands_separate', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title='mineral soils black; peatlands red', overlay_x=np.ma.masked_array(cru_data, mask=np.logical_not(organics_commonmap)).flatten(), overlay_y=np.ma.masked_array(mrt, mask=np.logical_not(organics_commonmap)).flatten() , overlay_dots=True, overlay_color='red')
    #
    #
    ################################################################################
    if load_thematic_maps:
        print('loading wise30sec_output file')
        print(' opening file: '+'datasets/gustafs_outputvars.nc')
        wise30sec_output_hwsd_file = nc.Dataset('datasets/gustafs_outputvars.nc')
        gf_lats = wise30sec_output_hwsd_file.variables['lat'][:]
        gf_lons = wise30sec_output_hwsd_file.variables['lon'][:]
        Ctot_0_20  = wise30sec_output_hwsd_file.variables['Ctot_0_20_wgs84_05.tif'][:]
        Ctot_100_150  = wise30sec_output_hwsd_file.variables['Ctot_100_150_wgs84_05.tif'][:]
        Ctot_100  = wise30sec_output_hwsd_file.variables['Ctot_100_wgs84_05.tif'][:]
        Ctot_150_200  = wise30sec_output_hwsd_file.variables['Ctot_150_200_wgs84_05.tif'][:]
        Ctot_200  = wise30sec_output_hwsd_file.variables['Ctot_200_wgs84_05.tif'][:]
        Ctot_20_40  = wise30sec_output_hwsd_file.variables['Ctot_20_40_wgs84_05.tif'][:]
        Ctot_40_60  = wise30sec_output_hwsd_file.variables['Ctot_40_60_wgs84_05.tif'][:]
        Ctot_60_80  = wise30sec_output_hwsd_file.variables['Ctot_60_80_wgs84_05.tif'][:]
        Ctot_80_100  = wise30sec_output_hwsd_file.variables['Ctot_80_100_wgs84_05.tif'][:]
        Fract_Peat  = wise30sec_output_hwsd_file.variables['Fract_Peat_wgs84_05.tif'][:]
        Fract_Pfmin  = wise30sec_output_hwsd_file.variables['Fract_Pfmin_wgs84_05.tif'][:]
        Fract_Shlw  = wise30sec_output_hwsd_file.variables['Fract_Shlw_wgs84_05.tif'][:]
        Fract_Water  = wise30sec_output_hwsd_file.variables['Fract_Water_wgs84_05.tif'][:]
        Fract_Weak  = wise30sec_output_hwsd_file.variables['Fract_Weak_wgs84_05.tif'][:]
        Fract_dune  = wise30sec_output_hwsd_file.variables['Fract_dune_wgs84_05.tif'][:]
        Fract_ice  = wise30sec_output_hwsd_file.variables['Fract_ice_wgs84_05.tif'][:]
        Fract_misc  = wise30sec_output_hwsd_file.variables['Fract_misc_wgs84_05.tif'][:]
        Fract_rock  = wise30sec_output_hwsd_file.variables['Fract_rock_wgs84_05.tif'][:]
        Fract_salt  = wise30sec_output_hwsd_file.variables['Fract_salt_wgs84_05.tif'][:]
        Fract_soil  = wise30sec_output_hwsd_file.variables['Fract_soil_wgs84_05.tif'][:]
        Fract_wet  = wise30sec_output_hwsd_file.variables['Fract_wet_wgs84_05.tif'][:]    
        wise30sec_output_hwsd_file.close()
        lat_offset_min = abs(lats_common - gf_lats.min()).argmin()
        Ctot_0_20_commongrid = np.zeros([JM_common,IM_common])
        Ctot_0_20_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_0_20[::-1,:]
        Ctot_100_150_commongrid = np.zeros([JM_common,IM_common])
        Ctot_100_150_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_100_150[::-1,:]
        Ctot_100_commongrid = np.zeros([JM_common,IM_common])
        Ctot_100_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_100[::-1,:]
        Ctot_150_200_commongrid = np.zeros([JM_common,IM_common])
        Ctot_150_200_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_150_200[::-1,:]
        Ctot_200_commongrid = np.zeros([JM_common,IM_common])
        Ctot_200_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_200[::-1,:]
        Ctot_20_40_commongrid = np.zeros([JM_common,IM_common])
        Ctot_20_40_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_20_40[::-1,:]
        Ctot_40_60_commongrid = np.zeros([JM_common,IM_common])
        Ctot_40_60_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_40_60[::-1,:]
        Ctot_60_80_commongrid = np.zeros([JM_common,IM_common])
        Ctot_60_80_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_60_80[::-1,:]
        Ctot_80_100_commongrid = np.zeros([JM_common,IM_common])
        Ctot_80_100_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Ctot_80_100[::-1,:]
        Fract_Peat_commongrid = np.zeros([JM_common,IM_common])
        Fract_Peat_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_Peat[::-1,:]
        Fract_Pfmin_commongrid = np.zeros([JM_common,IM_common])
        Fract_Pfmin_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_Pfmin[::-1,:]
        Fract_Shlw_commongrid = np.zeros([JM_common,IM_common])
        Fract_Shlw_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_Shlw[::-1,:]
        Fract_Water_commongrid = np.zeros([JM_common,IM_common])
        Fract_Water_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_Water[::-1,:]
        Fract_Weak_commongrid = np.zeros([JM_common,IM_common])
        Fract_Weak_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_Weak[::-1,:]
        Fract_dune_commongrid = np.zeros([JM_common,IM_common])
        Fract_dune_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_dune[::-1,:]
        Fract_ice_commongrid = np.zeros([JM_common,IM_common])
        Fract_ice_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_ice[::-1,:]
        Fract_misc_commongrid = np.zeros([JM_common,IM_common])
        Fract_misc_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_misc[::-1,:]
        Fract_rock_commongrid = np.zeros([JM_common,IM_common])
        Fract_rock_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_rock[::-1,:]
        Fract_salt_commongrid = np.zeros([JM_common,IM_common])
        Fract_salt_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_salt[::-1,:]
        Fract_soil_commongrid = np.zeros([JM_common,IM_common])
        Fract_soil_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_soil[::-1,:]
        Fract_wet_commongrid = np.zeros([JM_common,IM_common])
        Fract_wet_commongrid[lat_offset_min:lat_offset_min+len(gf_lats),:] = Fract_wet[::-1,:]
        ## try some plots of this now
        select_general = np.logical_not(common_mask)
        #
        ### also make peatlands mask using either hwsd or ncscd peatland fraction
        organics_commonmap_combined = np.logical_or(organics_commonmap, Fract_Peat_commongrid > max_hist_frac)
        map_funcs.xyplot(np.ma.masked_array(cru_data, mask=organics_commonmap_combined).flatten(), np.ma.masked_array(mrt, mask=organics_commonmap_combined).flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_peatlands_HWSD_NCSCD_separate', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x=np.ma.masked_array(cru_data, mask=np.logical_not(organics_commonmap_combined)).flatten(), overlay_y=np.ma.masked_array(mrt, mask=np.logical_not(organics_commonmap_combined)).flatten() , overlay_dots=True, overlay_color='red', labels=['Peat soils','Non-peat soils'], labelcolors=['red','black'], label_xstart=-15, label_ystart=5., label_yspace=2., labelfontsize=.03)
    #
    ###  and now the precip data from gpcc
    gpcc_filename = 'datasets/precip.mon.total.v6.nc'
    print( ' opening file '+gpcc_filename)
    gpcc_file = Nio.open_file(gpcc_filename)
    gpcc_lats = gpcc_file.variables['lat'][:]
    gpcc_lons = gpcc_file.variables['lon'][:]
    start_date_index = 960
    gpcc_data = gpcc_file.variables['precip'][start_date_index:,::-1,:].mean(axis=0) * 12
    gpcc_data_temp = gpcc_data.copy()
    gpcc_data[:,0:360] = gpcc_data_temp[:,360:]
    gpcc_data[:,360:] = gpcc_data_temp[:,0:360]
    #
    map_funcs.fill(gpcc_data, lats_common, lons_common, levels=np.arange(100,2500,100), file='map_gpcc_obs', title='GPCC mean precipitation (mm/yr)', specialprojection='global_noant')
    #
    common_mask = np.logical_or(cru_data.mask[:], mrt.mask[:])
    common_mask = np.logical_or(common_mask[:], gpcc_data.mask[:])
    common_mask = np.logical_or(common_mask[:], mrt < 1e-1)
    select = np.logical_not(common_mask)
    toskip=2
    #
    map_funcs.xyplot(cru_data[select].flatten()[::toskip].data, mrt[select].flatten()[::toskip].data, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_coloredbyprecip', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', shaded_dot_data=gpcc_data[select].flatten()[::toskip].data, shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, shuffle_shaded_dots=True)
    #
    map_funcs.fill(gpcc_data, lats_common, lons_common, levels=rain_colorlevels, colormap=rain_colormap, file='obs_gpcc_precip_map', specialprojection='global_noant', title='Precipitation (mm/yr)')
    #
    #
    ##############################
    # now load the PET data:  datasets/MOD16A3_Science_PET_mean_00_13_regridhalfdegree.nc
    modispet_datafilename = 'datasets/MOD16A3_Science_PET_mean_00_13_regridhalfdegree.nc'
    print( ' opening file '+modispet_datafilename)
    modispet_datafile = Nio.open_file(modispet_datafilename)
    modispet_data = modispet_datafile.variables['pet'][:]
    modispet_lats = modispet_datafile.variables['lat'][:]
    modispet_lons = modispet_datafile.variables['lon'][:]
    modispet_datafile.close()
    #
    common_mask_pet = np.logical_or(common_mask[:], modispet_data.mask[:])
    select = np.logical_not(common_mask_pet)
    toskip=2
    #
    pmpet_colorlevels = np.arange(-2000,2200,200)
    pmpet_colormap='MPL_RdYlBu'
    #
    map_funcs.xyplot(cru_data[select].flatten()[::toskip].data, mrt[select].flatten()[::toskip].data, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_coloredby_pmpet', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', shaded_dot_data=(gpcc_data - modispet_data)[select].flatten()[::toskip].data, shaded_dot_levels=pmpet_colorlevels, colormap=pmpet_colormap, shuffle_shaded_dots=True)
    #
    map_funcs.fill((gpcc_data - modispet_data), lats_common, lons_common, levels=pmpet_colorlevels, colormap=pmpet_colormap, file='obs_gpcc_pmpet_map', overlay_contour_data=(gpcc_data - modispet_data), overlay_contour_levels=[min_pmpet], overlay_contour_thickness=1.0, specialprojection='global_noant', title='Precip minus PET  (mm/yr)')
    #
    #
    ##############################
    # and now mask by minimum value of P minus PET
    common_mask_minpmpet = np.logical_or(common_mask_pet[:], (gpcc_data - modispet_data)[:] < min_pmpet)
    common_mask_minpmpet_notpeat = np.logical_or(common_mask_minpmpet[:], organics_commonmap_combined[:])
    select = np.logical_not(common_mask_minpmpet_notpeat)
    #
    #
    xdata_linreg_filtered = cru_data.flatten()[select.flatten()]
    ydata_linreg_filtered = np.log10(mrt.flatten()[select.flatten()])
    #
    ### try a polynomial regression
    quadratic_fit = np.polyfit(xdata_linreg_filtered,ydata_linreg_filtered,2)
    quadraticfunc = np.poly1d(quadratic_fit)
    residual_variance = np.var(ydata_linreg_filtered - quadraticfunc(xdata_linreg_filtered))
    textfile.write('quadratic_fit_obsquadratic_fit '+str(quadratic_fit)+', residual variance: '+str(residual_variance)+'\n')
    xd_logquadfit = np.linspace(xdata_linreg_filtered.min(), xdata_linreg_filtered.max(), 60)
    #
    yd_quad = 10. ** quadraticfunc(xd_logquadfit)
    yd_quad2 = 10. ** (quadratic_fit[2] + quadratic_fit[1]*xd_logquadfit + quadratic_fit[0]*xd_logquadfit**2 )
    slope_logspace = quadratic_fit[1] + 2 * quadratic_fit[0]*xd_logquadfit
    map_funcs.xyplot(cru_data[select].flatten(), mrt[select].flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_quadraticregression_filtered_pmpet', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title='filtered observations', overlay_x = xd_logquadfit, overlay_y = yd_quad2,overlay_color='red' )
    #
    #
    inferred_q10 = 10**(10.*(-slope_logspace))
    map_funcs.xyplot(xd_logquadfit, inferred_q10, file='inferred_Q10', xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Climatological Q~B~10~N~', aspect_ratio=1.5, xrange=temprange)
    #
    ### I want to try to do this in R to get confidence intervals
    xdata_rvect = robjects.FloatVector(xdata_linreg_filtered)
    ydata_rvect = robjects.FloatVector(ydata_linreg_filtered)
    #
    robjects.globalenv["xdata_rvect"] = xdata_rvect
    robjects.globalenv["ydata_rvect"] = ydata_rvect
    #
    quadreg_r = stats.lm("ydata_rvect ~ poly(xdata_rvect,2)")
    robjects.globalenv["quadreg_r"] = quadreg_r
    rconfint = robjects.r['confint']
    rpredict = robjects.r['predict']
    rsummary = robjects.r['summary']
    thepredictint = rpredict(quadreg_r, interval='prediction', level=0.50)
    print(rsummary(quadreg_r))
    predictint = np.array(thepredictint)
    n_toshow_predictlines = 50
    n_toskip_predictlines = xdata_linreg_filtered.shape[0]/n_toshow_predictlines
    indices_toshow = xdata_linreg_filtered.argsort()[::n_toskip_predictlines]
    x_confint = xdata_linreg_filtered[indices_toshow]
    y_confint = 10.**(predictint[indices_toshow,:].transpose())
    #
    map_funcs.xyplot(cru_data[select].flatten(), mrt[select].flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_quadraticregression_filtered_pmpet_r50pctpredint', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint, overlay_y = y_confint,overlay_color='red',overlay_linethickness=[2.5,1.,1.] )
    #
    np.savetxt('x_confint.txt', x_confint)
    np.savetxt('y_confint.txt', y_confint)
    #
    #
    #
    if read_clay_map:
        print('reading clay map')
        clayfilename = "datasets/PercClay_HWSD.nc"
        print(' opening file: '+clayfilename)
        clayfile = nc.Dataset(clayfilename)
        claylats = clayfile.variables['lat'][:]
        claylons = clayfile.variables['lon'][:]
        claymap_in = clayfile.variables['T_CLAY.tif'][:]
        claymap_halfdegree = conservative_regrid(claymap_in[::-1,:], claylats[::-1], claylons, lats_common, lons_common, centers_out=True, dilute_w_masked_data=False)
        #
        claylevels=np.arange(0,60,5)
        map_funcs.fill(claymap_halfdegree, lats_common, lons_common, levels=claylevels, file='clayperc_map', colormap='MPL_YlOrRd', specialprojection='global_noant', title='Mean clay percentage')
        #
        residual = np.log10(mrt) - quadraticfunc(cru_data)
        select_general = select.copy()
        map_funcs.xyplot(claymap_halfdegree, residual, dots=True, file='residual_logtau_vs_claypercent_commongrid', dotsize=0.002, xtitle='Clay Percentage', ytitle='residual log(tau)')
        toskip = 2
        map_funcs.xyplot(cru_data[select_general].flatten()[::toskip].data, mrt[select_general].flatten()[::toskip].data, dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='obs_MRT_soilc_temp_coloredby_clayperc', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', shaded_dot_data=claymap_halfdegree[select_general].flatten()[::toskip], shaded_dot_levels=claylevels, shuffle_shaded_dots=True, colormap='MPL_YlOrRd')
        map_funcs.xyplot(cru_data, claymap_halfdegree, dots=True, dotsize=coloreddotfig_dotsize/2., xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Clay percentage', regress=True, xrange=temprange, file='clay_vs_maat')
    #
    #
    # as a sanity check, make a plot of the central relationship and then a bunch of Q10 curves at points along the curve to make sure that they are tangent.
    # map_funcs.xyplot(xd_logquadfit, yd_quad2, ylog=True, yrange=[1., 5e3], xrange=temprange)
    ntangents = 15
    tangent_xcenters = np.linspace(-18.,28., ntangents)
    tangent_q10s = 10**(10.*(-(quadratic_fit[1] + 2 * quadratic_fit[0]*tangent_xcenters)))
    npoints_tangents = 30
    length_tangents = 100.
    tangentcurves_x = np.zeros([ntangents, npoints_tangents])
    tangentcurves_y = np.zeros([ntangents, npoints_tangents])
    for ii in range(ntangents):
        tangentcurves_x[ii,:] = np.linspace(tangent_xcenters[ii]-length_tangents/2.,tangent_xcenters[ii]+length_tangents/2.,npoints_tangents)
        tangentcurves_y[ii,:] = 10. ** quadraticfunc(tangent_xcenters[ii]) * (1./(tangent_q10s[ii]**((tangentcurves_x[ii,:]-tangent_xcenters[ii])/10.)))
    map_funcs.xyplot(xd_logquadfit, yd_quad2, ylog=True, yrange=[1., 5e3], xrange=temprange, overlay_x=tangentcurves_x, overlay_y=tangentcurves_y, overlay_color='red', file='Q10_tangent_calculation')
        
    
################################################################################
### now load some of this same info from CMIP5 models
### lets use the standard years of 1980-2005 where possible.
rmse_list = []
quadratic_regression_coefficients_list = []
if make_cmip5_plots:
    cmip5_data_in = [
        ['datasets/npp_Lmon_CCSM4_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/cSoil_Lmon_CCSM4_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/tas_Amon_CCSM4_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
        'CCSM4',
        'datasets/pr_Amon_CCSM4_historical_r1i1p1_185001-200512.nc', 130*12, 156*12],
        ['datasets/npp_Lmon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/cSoil_Lmon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/tas_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'MPI-ESM-LR',
        'datasets/pr_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12],
        ['datasets/npp_Lmon_GFDL-ESM2G_historical_r1i1p1_186101-200512.nc', 119*12, 145*12,
         'datasets/cSoil_Lmon_GFDL-ESM2G_historical_r1i1p1_199601-200512.nc', 0, 10*12,
         'datasets/tas_Amon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc', 0, 5*12,
         'GFDL-ESM2G',
        'datasets/pr_Amon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc', 0, 5*12],
        ['datasets/npp_Lmon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc', 0, 21*12,
         'datasets/cSoil_Lmon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc',  0, 21*12,
         'datasets/tas_Amon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc',  0, 21*12,
         'HadGEM2',
        'datasets/pr_Amon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc', 0, 21*12],
        ['datasets/npp_Lmon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/cSoil_Lmon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/tas_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'IPSL-CM5A-LR',
        'datasets/pr_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc', 130*12, 156*12],
        # ['datasets/npp_Lmon_CanESM2_historical_r1i1p1_185001-200512.nc',
        #  '',  soil C has zero size
        #  '',
        #  '']
        ['datasets/npp_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/cSoil_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'datasets/tas_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc', 130*12, 156*12,
         'MIROC-ESM',
        'datasets/pr_Amon_MIROC-ESM_historical_r1i1p1_185001-200512.nc', 130*12, 156*12]
         ]
    #
    n_models = len(cmip5_data_in)
    for model_i in range(n_models):
        print( ' opening file '+cmip5_data_in[model_i][0])
        nppfile_cmip5 = nc.Dataset(cmip5_data_in[model_i][0])
        lats_cmip5 = nppfile_cmip5.variables['lat'][:]
        lons_cmip5 = nppfile_cmip5.variables['lon'][:]
        npp_cmip5 = np.ma.masked_array(nppfile_cmip5.variables['npp'][cmip5_data_in[model_i][1]:cmip5_data_in[model_i][2],:,:].mean(axis=0))
        nppfile_cmip5.close()
        #
        print( ' opening file '+cmip5_data_in[model_i][3])
        cSoilfile_cmip5 = nc.Dataset(cmip5_data_in[model_i][3])
        cSoil_cmip5 = np.ma.masked_array(cSoilfile_cmip5.variables['cSoil'][cmip5_data_in[model_i][4]:cmip5_data_in[model_i][5],:,:].mean(axis=0))
        cSoilfile_cmip5.close()
        #
        print( ' opening file '+cmip5_data_in[model_i][6])
        tasfile_cmip5 = nc.Dataset(cmip5_data_in[model_i][6])
        tas_cmip5 = np.ma.masked_array(tasfile_cmip5.variables['tas'][cmip5_data_in[model_i][7]:cmip5_data_in[model_i][8],:,:].mean(axis=0))
        tasfile_cmip5.close()
        #
        print( ' opening file '+cmip5_data_in[model_i][10])
        prfile_cmip5 = nc.Dataset(cmip5_data_in[model_i][10])
        pr_cmip5 = np.ma.masked_array(prfile_cmip5.variables['pr'][cmip5_data_in[model_i][11]:cmip5_data_in[model_i][12],:,:].mean(axis=0) * 86400. * 365.)
        prfile_cmip5.close()
        #
        #
        yvar =  cSoil_cmip5/(npp_cmip5 * 86400.*365.)
        yvar_masked = np.ma.masked_array(yvar, mask=np.logical_or(yvar < 1, yvar > 1e4))
        xvar = tas_cmip5-273.15
        map_funcs.xyplot(xvar.flatten(), yvar_masked.flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_CMIP5_'+cmip5_data_in[model_i][9], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title=cmip5_data_in[model_i][9])
        #
        #
        hist_x = xvar.flatten()
        hist_y = np.log10(yvar_masked.flatten())
        #areas_cmip5 = geog_funcs.gridcell_areas(lats_cmip5, lons_cmip5, mask=np.logical_or(np.ma.masked_array(tas_cmip5[:]).mask[:], np.ma.masked_array(cSoil_cmip5[:]).mask[:])).flatten()
        #
        ### make common mask for data
        commonmask = np.logical_or(hist_x[:].mask, hist_y[:].mask)
        hist_x = hist_x[np.logical_not(commonmask)]
        hist_y = hist_y[np.logical_not(commonmask)]
        #areas_cmip5 = areas_cmip5[np.logical_not(commonmask)]
        #
        #
        map_funcs.fill(yvar_masked, lats_cmip5, lons_cmip5, levels=log_mrt_levels, file='map_mrt_'+cmip5_data_in[model_i][9], title='MRT (yr)')
        #
        map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_CMIP5_coloredbyprecip_'+cmip5_data_in[model_i][9], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title=cmip5_data_in[model_i][9]+'; Colored by Precip', shaded_dot_data=pr_cmip5[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, shuffle_shaded_dots=True)
        #
        ### plot with overlay of the observed data trend
        if make_data_plots:
            map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_CMIP5_coloredbyprecip_'+cmip5_data_in[model_i][9]+'_withobscurve', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', inset_title=cmip5_data_in[model_i][9], inset_title_x=25., inset_textjust="CenterRight", inset_title_y=2.5e3, shaded_dot_data=pr_cmip5[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, overlay_x = x_confint, overlay_y = y_confint,overlay_color='black',overlay_linethickness=[2.5,1.,1.], shuffle_shaded_dots=True)
            #
            # regrid the modispet_data
            modispet_data_regrid_cmip5 = conservative_regrid(modispet_data[:], modispet_lats[:], modispet_lons, lats_cmip5, lons_cmip5, centers_out=True, dilute_w_masked_data=False)
            p_minus_pet_cmip5 = pr_cmip5 - modispet_data_regrid_cmip5
            #
            p_minus_pet_mask = np.logical_or(np.logical_or(yvar_masked.mask, (p_minus_pet_cmip5 < min_pmpet)), p_minus_pet_cmip5.mask)
            xvar_petmask = np.ma.masked_array(xvar, mask=p_minus_pet_mask)
            yvar_petmask = np.ma.masked_array(yvar_masked, mask=p_minus_pet_mask)
            #
            xdata_linreg_filtered = xvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:]
            ydata_linreg_filtered = np.log10(yvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:])
            #
            #
            # do the quadtratic fit using numpy directly to get coefficients
            quadratic_fit_cmip5 = np.polyfit(xdata_linreg_filtered,ydata_linreg_filtered,2)
            print(quadratic_fit_cmip5)
            quadratic_regression_coefficients_list.append([cmip5_data_in[model_i][9],quadratic_fit_cmip5])
            quadraticfunc_cmip5 = np.poly1d(quadratic_fit_cmip5)
            residual_variance_cmip5 =np.var(ydata_linreg_filtered - quadraticfunc_cmip5(xdata_linreg_filtered))
            textfile.write(cmip5_data_in[model_i][9]+' '+str(quadratic_fit_cmip5)+', residual variance: '+str(residual_variance_cmip5)+'\n')
            #
            ### next do the same regression as the data and above, but also in R to get prediction intervals
            xdata_rvect_cmip5 = robjects.FloatVector(xdata_linreg_filtered)
            ydata_rvect_cmip5 = robjects.FloatVector(ydata_linreg_filtered)
            #
            robjects.globalenv["xdata_rvect_cmip5"] = xdata_rvect_cmip5
            robjects.globalenv["ydata_rvect_cmip5"] = ydata_rvect_cmip5
            #
            quadreg_r_cmip5 = stats.lm("ydata_rvect_cmip5 ~ poly(xdata_rvect_cmip5,2)")
            robjects.globalenv["quadreg_r_cmip5"] = quadreg_r_cmip5
            rconfint_cmip5 = robjects.r['confint']
            rpredict_cmip5 = robjects.r['predict']
            rsummary_cmip5 = robjects.r['summary']
            thepredictint_cmip5 = rpredict_cmip5(quadreg_r_cmip5, interval='prediction', level=0.50)
            print(rsummary_cmip5(quadreg_r_cmip5))
            predictint_cmip5 = np.array(thepredictint_cmip5)
            n_toshow_predictlines_cmip5 = 50
            n_toskip_predictlines_cmip5 = xdata_linreg_filtered.shape[0]/n_toshow_predictlines_cmip5
            indices_toshow_cmip5 = xdata_linreg_filtered.argsort()[::n_toskip_predictlines_cmip5]
            x_confint_cmip5 = xdata_linreg_filtered[indices_toshow_cmip5]
            y_confint_cmip5 = 10.**(predictint_cmip5[indices_toshow_cmip5,:].transpose())
            #
            select = np.logical_not(p_minus_pet_mask)
            #
            map_funcs.xyplot(xvar[select].flatten(), yvar_masked[select].flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='cmip5_MRT_soilc_temp_quadraticregression_filtered_pmpet_r50pctpredint_'+cmip5_data_in[model_i][9], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_cmip5, overlay_y = y_confint_cmip5,overlay_color='red',overlay_linethickness=[2.5,1.,1.] )
            #
            np.savetxt('x_confint_'+cmip5_data_in[model_i][9]+'.txt', x_confint)
            np.savetxt('y_confint_'+cmip5_data_in[model_i][9]+'.txt', y_confint)
            #
            #
            #
            ### calculate the rmse of the mean model-predicted value per temperature increment
            mean_y_cmip = np.zeros(n_T_bins)
            mean_x_cmip = np.zeros(n_T_bins)
            count_Tbin_cmip = np.zeros(n_T_bins)
            for i_T in range(n_T_bins):
                T_bin_mask = np.logical_not(np.logical_and(xvar_petmask > (T_bin_centers[i_T] - Tbin_width/2.), xvar_petmask <= (T_bin_centers[i_T] + Tbin_width/2.)))
                mean_y_cmip[i_T] = np.log10(np.ma.masked_array(yvar_petmask, mask=T_bin_mask)).mean()
                mean_x_cmip[i_T] = np.ma.masked_array(xvar_petmask, mask=T_bin_mask).mean()
                count_Tbin_cmip[i_T] = np.ma.count(np.ma.masked_array(xvar_petmask, mask=T_bin_mask))
            #
            rmse = np.mean((mean_y_cmip - quadraticfunc(T_bin_centers))**2.)**(0.5)
            print(cmip5_data_in[model_i][9] + ' RMSE: ' + str(rmse))
            rmse_list.append([cmip5_data_in[model_i][9], rmse])
            #
            #
        #
        #




################################################################################################################################################################
### now open up clm4.5 files
################################################################################################################################################################
if make_clm4_5_plots:
    #
    clm45_input_filenames = [
        ['clm45bgc_hrv_1deg4508_hist.NPP.monthly.1850-2005.nc',
         'clm45bgc_hrv_1deg4508_hist.TOTSOMC_1m.annual.1850-2005.nc',
         'clm45bgc_hrv_1deg4508_hist.TSA.monthly.1850-2005.nc',
         'clm45bgc_hrv_1deg4508_hist.RAIN.monthly.1850-2005.nc',
         'clm45bgc_hrv_1deg4508_hist.SNOW.monthly.1850-2005.nc'],
        ['clm45bgc_hrv_cdk_10DDD_1deg4508_hist.NPP.monthly.1850-2005.nc',
         'clm45bgc_hrv_cdk_10DDD_1deg4508_hist.TOTSOMC_1m.annual.1850-2005.nc',
         'clm45bgc_hrv_cdk_10DDD_1deg4508_hist.TSA.monthly.1850-2005.nc',
         'clm45bgc_hrv_cdk_10DDD_1deg4508_hist.RAIN.monthly.1850-2005.nc',
         'clm45bgc_hrv_cdk_10DDD_1deg4508_hist.SNOW.monthly.1850-2005.nc'],
         ]
    titles = ['CLM4.5; Z~B~~723~~N~ = 0.5m','CLM4.5; Z~B~~723~~N~ = 10m']
    case_names = ['clm45_zt05m','clm45_zt10m']
    clm45_datadir = 'datasets/'
    for model_i in range(len(clm45_input_filenames)):
        #
        clm45_time_strt_monthly = 125 * 12
        clm45_time_end_monthly = 155 * 12
        #
        clm45_time_strt_annual = 125
        clm45_time_end_annual = 155
        #
        print( ' opening file '+clm45_datadir+clm45_input_filenames[model_i][0])
        nppfile_clm45 = Nio.open_file(clm45_datadir+clm45_input_filenames[model_i][0])
        lats_clm45 = nppfile_clm45.variables['lat'][:]
        lons_clm45 = nppfile_clm45.variables['lon'][:]
        npp_clm45 = nppfile_clm45.variables['NPP'][clm45_time_strt_monthly:clm45_time_end_monthly,:,:].mean(axis=0) * 86400.*365.
        nppfile_clm45.close()
        #
        print( ' opening file '+clm45_datadir+clm45_input_filenames[model_i][1])
        cSoilfile_clm45 = Nio.open_file(clm45_datadir+clm45_input_filenames[model_i][1])
        cSoil_clm45 = cSoilfile_clm45.variables['TOTSOMC_1m'][clm45_time_strt_annual:clm45_time_end_annual,:,:].mean(axis=0)
        cSoilfile_clm45.close()
        #
        print( ' opening file '+clm45_datadir+clm45_input_filenames[model_i][2])
        tasfile_clm45 = Nio.open_file(clm45_datadir+clm45_input_filenames[model_i][2])
        tsa_clm45 = tasfile_clm45.variables['TSA'][clm45_time_strt_monthly:clm45_time_end_monthly,:,:].mean(axis=0)
        tasfile_clm45.close()
        #
        print( ' opening file '+clm45_datadir+clm45_input_filenames[model_i][3])
        rainfile_clm45 = Nio.open_file(clm45_datadir+clm45_input_filenames[model_i][3])
        rain_clm45 = rainfile_clm45.variables['RAIN'][clm45_time_strt_monthly:clm45_time_end_monthly,:,:].mean(axis=0) * 86400. * 365.
        rainfile_clm45.close()
        #
        print( ' opening file '+clm45_datadir+clm45_input_filenames[model_i][4])
        snowfile_clm45 = Nio.open_file(clm45_datadir+clm45_input_filenames[model_i][4])
        snow_clm45 = snowfile_clm45.variables['SNOW'][clm45_time_strt_monthly:clm45_time_end_monthly,:,:].mean(axis=0) * 86400. * 365.
        snowfile_clm45.close()
        #
        precip_clm45 = (rain_clm45 + snow_clm45)
        #
        yvar =  cSoil_clm45/(npp_clm45)
        print(yvar.mean())
        yvar_masked = np.ma.masked_array(yvar, mask=yvar < .00001)
        xvar = tsa_clm45-273.15
        map_funcs.xyplot(xvar.flatten(), yvar_masked.flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_'+case_names[model_i], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title=titles[model_i])
        #
        hist_x = xvar.flatten()
        hist_y = np.log10(yvar_masked.flatten())
        areas_clm45 = geog_funcs.gridcell_areas(lats_clm45, lons_clm45, mask=np.logical_or(tsa_clm45.mask[:], cSoil_clm45.mask[:])).flatten()
        #
        ### make common mask for data
        commonmask = np.logical_or(hist_x.mask[:], hist_y.mask[:])
        hist_x = hist_x[np.logical_not(commonmask)]
        hist_y = hist_y[np.logical_not(commonmask)]
        areas_clm45 = areas_clm45[np.logical_not(commonmask)]
        #
        #
        map_funcs.fill(yvar_masked, lats_clm45, lons_clm45, levels=log_mrt_levels, file='map_mrt_'+case_names[model_i], title=titles[model_i]+' Inferred Turnover Time (yr)')
        #
        map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_coloredbyprecip_'+case_names[model_i], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title=titles[model_i]+'; Colored by Precip', shaded_dot_data=precip_clm45[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, shuffle_shaded_dots=True)
        #
        if make_data_plots:
            map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_coloredbyprecip_'+case_names[model_i]+'_withobscurve', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', inset_title=titles[model_i], inset_title_x=25., inset_textjust="CenterRight", inset_title_y=2.5e3, shaded_dot_data=precip_clm45[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, overlay_x = x_confint, overlay_y = y_confint,overlay_color='black',overlay_linethickness=[2.5,1.,1.], shuffle_shaded_dots=True)
            #
            #
            # regrid the modispet_data
            modispet_data_regrid_clm45 = conservative_regrid(modispet_data[:], modispet_lats[:], modispet_lons, lats_clm45, lons_clm45, centers_out=True, dilute_w_masked_data=False)
            p_minus_pet_clm45 = precip_clm45 - modispet_data_regrid_clm45
            #
            p_minus_pet_mask = np.logical_or(np.logical_or(yvar_masked.mask, (p_minus_pet_clm45 < min_pmpet)), p_minus_pet_clm45.mask)
            xvar_petmask = np.ma.masked_array(xvar, mask=p_minus_pet_mask)
            yvar_petmask = np.ma.masked_array(yvar_masked, mask=p_minus_pet_mask)
            #
            xdata_linreg_filtered = xvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:]
            ydata_linreg_filtered = np.log10(yvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:])
            #
            #
            # do the quadtratic fit using numpy directly to get coefficients
            quadratic_fit_clm45 = np.polyfit(xdata_linreg_filtered,ydata_linreg_filtered,2)
            print(quadratic_fit_clm45)
            quadratic_regression_coefficients_list.append([case_names[model_i],quadratic_fit_clm45])
            quadraticfunc_clm45 = np.poly1d(quadratic_fit_clm45)
            residual_variance_clm45 =np.var(ydata_linreg_filtered - quadraticfunc_clm45(xdata_linreg_filtered))
            textfile.write(case_names[model_i]+' '+str(quadratic_fit_clm45)+', residual variance: '+str(residual_variance_clm45)+'\n')
            #
            ### next do the same regression as the data and above, but also in R to get prediction intervals
            xdata_rvect_clm45 = robjects.FloatVector(xdata_linreg_filtered)
            ydata_rvect_clm45 = robjects.FloatVector(ydata_linreg_filtered)
            #
            robjects.globalenv["xdata_rvect_clm45"] = xdata_rvect_clm45
            robjects.globalenv["ydata_rvect_clm45"] = ydata_rvect_clm45
            #
            quadreg_r_clm45 = stats.lm("ydata_rvect_clm45 ~ poly(xdata_rvect_clm45,2)")
            robjects.globalenv["quadreg_r_clm45"] = quadreg_r_clm45
            rconfint_clm45 = robjects.r['confint']
            rpredict_clm45 = robjects.r['predict']
            rsummary_clm45 = robjects.r['summary']
            thepredictint_clm45 = rpredict_clm45(quadreg_r_clm45, interval='prediction', level=0.50)
            print(rsummary_clm45(quadreg_r_clm45))
            predictint_clm45 = np.array(thepredictint_clm45)
            n_toshow_predictlines_clm45 = 50
            n_toskip_predictlines_clm45 = xdata_linreg_filtered.shape[0]/n_toshow_predictlines_clm45
            indices_toshow_clm45 = xdata_linreg_filtered.argsort()[::n_toskip_predictlines_clm45]
            x_confint_clm45 = xdata_linreg_filtered[indices_toshow_clm45]
            y_confint_clm45 = 10.**(predictint_clm45[indices_toshow_clm45,:].transpose())
            #
            select = np.logical_not(p_minus_pet_mask)
            #
            map_funcs.xyplot(xvar[select].flatten(), yvar_masked[select].flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='clm45_MRT_soilc_temp_quadraticregression_filtered_pmpet_r50pctpredint_'+case_names[model_i], dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_clm45, overlay_y = y_confint_clm45,overlay_color='red',overlay_linethickness=[2.5,1.,1.] )
            #
            np.savetxt('x_confint_'+case_names[model_i]+'.txt', x_confint)
            np.savetxt('y_confint_'+case_names[model_i]+'.txt', y_confint)
            #
            #
            #
            ### calculate the rmse of the mean model-predicted value per temperature increment
            mean_y_cmip = np.zeros(n_T_bins)
            mean_x_cmip = np.zeros(n_T_bins)
            count_Tbin_cmip = np.zeros(n_T_bins)
            for i_T in range(n_T_bins):
                T_bin_mask = np.logical_not(np.logical_and(xvar_petmask > (T_bin_centers[i_T] - Tbin_width/2.), xvar_petmask <= (T_bin_centers[i_T] + Tbin_width/2.)))
                mean_y_cmip[i_T] = np.log10(np.ma.masked_array(yvar_petmask, mask=T_bin_mask)).mean()
                mean_x_cmip[i_T] = np.ma.masked_array(xvar_petmask, mask=T_bin_mask).mean()
                count_Tbin_cmip[i_T] = np.ma.count(np.ma.masked_array(xvar_petmask, mask=T_bin_mask))
            #
            rmse = np.mean((mean_y_cmip - quadraticfunc(T_bin_centers))**2.)**(0.5)
            print(case_names[model_i] + ' RMSE: ' + str(rmse))
            rmse_list.append([case_names[model_i], rmse])
            #
            #
            #
            #
    #
    overlay_levels = np.array([.002,0.02, 0.2,  2.])
    #


if make_mimics_plots:
    mimics_filename = 'datasets/MIMICS_MRT.nc'
    print( ' opening file '+mimics_filename)
    mimics_file = Nio.open_file(mimics_filename)
    lats_mimics = mimics_file.variables['lat'][:]
    lons_mimics = mimics_file.variables['lon'][:]
    npp_mimics = mimics_file.variables['NPP'][:]
    cSoil_mimics = mimics_file.variables['MIM_SOM'][:]
    tsa_mimics = mimics_file.variables['TSA'][:]
    pr_mimics = mimics_file.variables['MAP'][:]
    #
    yvar =  cSoil_mimics/(npp_mimics)
    yvar_masked = np.ma.masked_array(yvar, mask=yvar < .00001)
    xvar = tsa_mimics
    map_funcs.xyplot(xvar.flatten(), yvar_masked.flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_mimics', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title='MIMICS')
    #
    map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_coloredbyprecip_mimics', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', title='MIMICS; Colored by Precip', shaded_dot_data=pr_mimics[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, shuffle_shaded_dots=True)
    #
    if make_data_plots:
        map_funcs.xyplot(xvar[np.logical_not(yvar_masked.mask[:])].flatten().data[:], yvar_masked[np.logical_not(yvar_masked.mask[:])].flatten().data[:], dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='modeled_MRT_soilc_temp_coloredbyprecip_mimics_withobscurve', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', inset_title='MIMICS', inset_title_x=25., inset_textjust="CenterRight", inset_title_y=2.5e3, shaded_dot_data=pr_mimics[np.logical_not(yvar_masked.mask[:])].flatten().data[:], shaded_dot_levels=rain_colorlevels, colormap=rain_colormap, overlay_x = x_confint, overlay_y = y_confint,overlay_color='black',overlay_linethickness=[2.5,1.,1.], shuffle_shaded_dots=True)
        #
        # regrid the modispet_data
        modispet_data_regrid_mimics = conservative_regrid(modispet_data[:], modispet_lats[:], modispet_lons, lats_mimics, lons_mimics, centers_out=True, dilute_w_masked_data=False)
        p_minus_pet_mimics = pr_mimics - modispet_data_regrid_mimics
        #
        p_minus_pet_mask = np.logical_or(np.logical_or(yvar_masked.mask, (p_minus_pet_mimics < min_pmpet)), p_minus_pet_mimics.mask)
        xvar_petmask = np.ma.masked_array(xvar, mask=p_minus_pet_mask)
        yvar_petmask = np.ma.masked_array(yvar_masked, mask=p_minus_pet_mask)
        #
        xdata_linreg_filtered = xvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:]
        ydata_linreg_filtered = np.log10(yvar_petmask[np.logical_not(xvar_petmask.mask[:])].flatten().data[:])
        #
        #
        # do the quadtratic fit using numpy directly to get coefficients
        quadratic_fit_mimics = np.polyfit(xdata_linreg_filtered,ydata_linreg_filtered,2)
        print(quadratic_fit_mimics)
        quadratic_regression_coefficients_list.append(['mimics',quadratic_fit_mimics])
        quadraticfunc_mimics = np.poly1d(quadratic_fit_mimics)
        residual_variance_mimics =np.var(ydata_linreg_filtered - quadraticfunc_mimics(xdata_linreg_filtered))
        textfile.write('mimics ' +str(quadratic_fit_mimics)+', residual variance: '+str(residual_variance_mimics)+'\n')
        #
        ### next do the same regression as the data and above, but also in R to get prediction intervals
        xdata_rvect_mimics = robjects.FloatVector(xdata_linreg_filtered)
        ydata_rvect_mimics = robjects.FloatVector(ydata_linreg_filtered)
        #
        robjects.globalenv["xdata_rvect_mimics"] = xdata_rvect_mimics
        robjects.globalenv["ydata_rvect_mimics"] = ydata_rvect_mimics
        #
        quadreg_r_mimics = stats.lm("ydata_rvect_mimics ~ poly(xdata_rvect_mimics,2)")
        robjects.globalenv["quadreg_r_mimics"] = quadreg_r_mimics
        rconfint_mimics = robjects.r['confint']
        rpredict_mimics = robjects.r['predict']
        rsummary_mimics = robjects.r['summary']
        thepredictint_mimics = rpredict_mimics(quadreg_r_mimics, interval='prediction', level=0.50)
        print(rsummary_mimics(quadreg_r_mimics))
        predictint_mimics = np.array(thepredictint_mimics)
        n_toshow_predictlines_mimics = 50
        n_toskip_predictlines_mimics = xdata_linreg_filtered.shape[0]/n_toshow_predictlines_mimics
        indices_toshow_mimics = xdata_linreg_filtered.argsort()[::n_toskip_predictlines_mimics]
        x_confint_mimics = xdata_linreg_filtered[indices_toshow_mimics]
        y_confint_mimics = 10.**(predictint_mimics[indices_toshow_mimics,:].transpose())
        #
        select = np.logical_not(p_minus_pet_mask)
        #
        map_funcs.xyplot(xvar[select].flatten(), yvar_masked[select].flatten(), dots=True, ylog=True, yrange=[1., 5e3], xrange=temprange, file='mimics_MRT_soilc_temp_quadraticregression_filtered_pmpet_r50pctpredint_'+'mimics', dotsize=coloreddotfig_dotsize, xtitle='Mean Air Temperature (~S~o~N~C)', ytitle='Inferred Turnover Time (yr)', overlay_x = x_confint_mimics, overlay_y = y_confint_mimics,overlay_color='red',overlay_linethickness=[2.5,1.,1.] )
        #
        np.savetxt('x_confint_'+'mimics'+'.txt', x_confint)
        np.savetxt('y_confint_'+'mimics'+'.txt', y_confint)
        #
        #
        #
        ### calculate the rmse of the mean model-predicted value per temperature increment
        mean_y_cmip = np.zeros(n_T_bins)
        mean_x_cmip = np.zeros(n_T_bins)
        count_Tbin_cmip = np.zeros(n_T_bins)
        for i_T in range(n_T_bins):
            T_bin_mask = np.logical_not(np.logical_and(xvar_petmask > (T_bin_centers[i_T] - Tbin_width/2.), xvar_petmask <= (T_bin_centers[i_T] + Tbin_width/2.)))
            mean_y_cmip[i_T] = np.log10(np.ma.masked_array(yvar_petmask, mask=T_bin_mask)).mean()
            mean_x_cmip[i_T] = np.ma.masked_array(xvar_petmask, mask=T_bin_mask).mean()
            count_Tbin_cmip[i_T] = np.ma.count(np.ma.masked_array(xvar_petmask, mask=T_bin_mask))
        #
        rmse = np.mean((mean_y_cmip - quadraticfunc(T_bin_centers))**2.)**(0.5)
        print('mimics' + ' RMSE: ' + str(rmse))
        rmse_list.append(['mimics', rmse])
        #
        #
        #
        #
    
textfile.close()
################################################################################
