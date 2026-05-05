import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point, Polygon
import glob
import pyart
import pandas as pd
import scipy
from netCDF4 import Dataset
import wrf
from metpy.calc import density
from metpy.units import units
from scipy.optimize import curve_fit
import coordinateSystems as cs
from scipy import interpolate
import sys

dirfiledict = {'morrison_noah':'mn','thom_ruc':'tr','thom_noah':'tn'}

def centers_to_edges(xcenter,ycenter):
    edgex = np.zeros(len(xcenter)+1)
    edgey = np.zeros(len(ycenter)+1)
    edgex[0] = xcenter[0]-(xcenter[1]-xcenter[0])/2
    edgey[0] = ycenter[0]-(ycenter[1]-ycenter[0])/2
    edgex[1:] = xcenter+(xcenter[1]-xcenter[0])/2
    edgey[1:] = ycenter+(ycenter[1]-ycenter[0])/2
    return edgex,edgey

def model_loc_of_station(station,run_dir='thom_ruc',return_ll=False):
    # get station latitude and longitude from station file
    stns = np.loadtxt("/nfs/turbo/seas-hutsona/NEXRAD/stn_locations.txt",
                      delimiter=",", dtype=str)
    iloc = np.argwhere(stns[:,0]==station)[0][0]
    llat,llon = stns[iloc,1:].astype(float)

    wrffile = f'/nfs/turbo/seas-hutsona/les_wrf_runs/2022jan/{run_dir}/wrfout_d01_2022-01-18_12:00:00'
    fh      = xr.open_dataset(wrffile)
    lat,lon = fh.XLAT.data[0],fh.XLONG.data[0]
    fh.close()

    stn_idx = np.argwhere((np.abs(lat-llat)<0.013)&(np.abs(lon-llon)<0.02))[0]
    if return_ll is False:
        return stn_idx
    if return_ll is True:
        return stn_idx, llat, llon
    

def cartopy_gl_map(fig,nrows,ncols,ax_idx,leftlbls=[],botlbls=[],rightlbls=[],feed_ax=None,
                   center_lat_lon=None,extent_og=None):
    ''' Create ax handles for subplots containing a map of the Great Lakes Region.
        nrows: number of rows
        ncols: number of columns
        ax_idx: specifies index of subplot (i.e., this func. needs to be called for each subplot)
        leftlbls: which subplot indices will have ticklabels on the left?
        botlbls: which subplot indices will have ticklabels on the bottom?
        center_lat_lon: array containing [lat,lon] of center of map (default is mathematical center of extent)
        extent_og: a new set of map bounds, other than just the GLR 
                   [wlon,elon,slat,nlat]
        OUTPUT: ax (the axis handle for the subplot) and datacrs (needed for the actual data plotting)
    '''
    if extent_og == None:
        extent = [-93,-76, 40, 50]
    else:
        extent = extent_og
    central_lat = np.mean(extent[2:])#42.5
    central_lon = np.mean(extent[:2]) #-85
        
    if center_lat_lon != None:
        central_lat = center_lat_lon[0]
        central_lon = center_lat_lon[1]
        if extent_og == None:
            extent = [central_lon-2,central_lon+2,central_lat-1.5,central_lat+1.5]

    # Set up plot crs (mapcrs) and the data crs, will need to transform all variables
    mapcrs = ccrs.LambertConformal(central_longitude=central_lon, central_latitude=central_lat,
                                   standard_parallels=(30, 60))
    datacrs = ccrs.PlateCarree()

    # FIGURE HANDLE
    if feed_ax == None:
        ax  = fig.add_subplot(nrows, ncols, ax_idx+1, projection=mapcrs)
    else:
        ax=feed_ax
        
    ax.set_extent(extent, crs=datacrs)
    ax.coastlines(resolution='50m')
    g = ax.gridlines(color=None,linewidth=0.5,linestyle=':',
                     draw_labels=True,x_inline=False,y_inline=False,
                     zorder=1000) #inlines=False keeps numbers out of figure
    g.right_labels =False
    g.top_labels   =False
    g.bottom_labels=False
    g.left_labels  =False

    if ax_idx in leftlbls:
        g.left_labels=True
    if ax_idx in botlbls:
        g.bottom_labels=True
    if ax_idx in rightlbls:
        g.right_labels=True

    lakes_50m=cfeature.NaturalEarthFeature('physical','lakes','50m',edgecolor='0.0',facecolor='none')
    ax.add_feature(lakes_50m,zorder=2000)
    return ax,datacrs

# find radar files closest to datetime
def find_radar_files(time,stations,nexrad_dir='/nfs/turbo/seas-ayumif/hutsona/NEXRAD/'):

    if len(stations) == 1:
        single_radar = True
    elif len(stations) > 1:
        single_radar = False
    else:
        print('no radar station sites are given')
    
    radar_files = []
    print(time)
    yr = time.year
    mth= time.month
    day= time.day
    hr = time.hour
    mn = time.minute
    for station in stations:
        stn_list = glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{hr:02d}{mn:02d}*')
        if len(stn_list)==1:
            radar_files.append(stn_list[0])

            if 'MDM.ar2v' in stn_list[0]:
                print('not the file we want, though')
        elif len(stn_list)>1:
            for f,fn in enumerate(stn_list):
                if 'MDM.ar2v' in fn:
                    continue
                else:
                    radar_files.append(stn_list[f])
        elif len(stn_list)==0:
            hourfiles = []
            if time.minute > 50:
                hourfiles = np.hstack((hourfiles,(sorted(glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{hr:02d}*')))))
                hourfiles = np.hstack((hourfiles,(sorted(glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{(hr+1):02d}*')))))
            elif time.minute < 10:
                hourfiles = np.hstack((hourfiles,(sorted(glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{(hr-1):02d}*')))))
                hourfiles = np.hstack((hourfiles,(sorted(glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{hr:02d}*')))))
            else:
                hourfiles = sorted(glob.glob(nexrad_dir+f'{station}/{station}{yr}{mth:02d}{day:02d}_{hr:02d}*'))
            if len(hourfiles)==0:
                print(f'The {station} radar was down within an hour of this time')
                continue 

            time_dist = []
            for hf in hourfiles:
                time_str = hf.split('/')[-1].split('_')[1]
                file_time= pd.to_datetime(f'{yr}{mth:02d}{day:02d} {time_str[:2]}{time_str[2:4]}{time_str[4:]}')
                time_dist.append((file_time-time).total_seconds())
            time_dist = np.asarray(time_dist)
            flag = False
            for t in time_dist:
                if t <= (60*30):
                    file_idx  = np.argmin(np.abs(time_dist))
                    radar_files.append(hourfiles[file_idx])
                    flag = True
                    break

            if flag is False:
                print('There are no radar files within 30 minutes of this time')
    return radar_files

def get_composite_radar(time,stations,radar_file=None,max_refl=False):
    
    if radar_file is not None:
        single_radar=True
        radar_files = [radar_file]
    else:
        radar_files = find_radar_files(time,stations)
        if len(radar_files)<1:
            print('cannot make composite radar plot for this time')
            return False

        if len(stations) == 1:
            single_radar = True
        elif len(stations) > 1:
            single_radar = False
        else:
            print('no radar station sites are given')
    
    gatefilters, radars, coords =[],[],[]
    for file in radar_files:
        print(file)
        if radar_file is not None:
            radar = file
        else:
            radar = pyart.io.read_nexrad_archive(file)
        coords.append([radar.latitude['data'][0],radar.longitude['data'][0]])
        radars.append(radar)
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilters.append(gatefilter)

    gatefilters = tuple(gatefilters)
    radars = tuple(radars)
    
    if single_radar is True:
        grid_limits = ((0,1000),(-250000,250000),(-300000,300000))
        grid_origin = (radar.latitude['data'][0],radar.longitude['data'][0])
    if single_radar is False:
        grid_limits = ((0,1000),(-500000,500000),(-1000000,1000000))
        grid_origin = (45,-85)
    
    if max_refl is True:
        grid = pyart.map.grid_from_radars(radars, gatefilters=gatefilters,
                                          grid_shape=(1, 351, 351),
                                          grid_limits= grid_limits,
                                          grid_origin=grid_origin,
                                          fields = ['max_refl'])
    if max_refl is False:
        grid = pyart.map.grid_from_radars(radars, gatefilters=gatefilters,
                                          grid_shape=(1, 351, 351),
                                          grid_limits= grid_limits,
                                          grid_origin=grid_origin,
                                          fields = ['reflectivity'])
    print('Radar grid object done')
    return grid

def get_wrf_for_psd(pdtime,station,run_dir,time=0,reflectivity=False, shift=None):
    '''
    This function was designed to grab the specific variables needed to calculate the particle size distributions from the microphysics schemes. Change directory locations when using this code. 

    INPUT:
    pdtime - a pandas datetime from which the variables are needed
    station - corresponds to a specific predetermined lat/lon location on the model grid
    run_dir - directory containing the wrf files. It also tells the function which microphysics scheme will be used for the PSD. 
    time - index of time output for each wrf file (default: 0)
    reflectivity - If True, the function will return the reflectivity value at the same location
    shift - if not None, will shift the results to a grid point around the station of interest

    '''
    # read in wrf variables at surface near station
    
    # get station latitude and longitude from station file
    stns = np.loadtxt("/nfs/turbo/seas-hutsona/NEXRAD/stn_locations.txt",
                      delimiter=",", dtype=str)
    iloc = np.argwhere(stns[:,0]==station)[0][0]
    llat,llon = stns[iloc,1:].astype(float)

    
    datedir_dict = {1:'2022jan', 11:'2022nov', 12:'2022dec'}
    dir = f'/nfs/turbo/seas-hutsona/les_wrf_runs/{datedir_dict[pdtime.month]}/{run_dir}/'
    wrffile = dir+f'wrfout_d01_2022-{pdtime.month:02d}-{pdtime.day:02d}_{pdtime.hour:02d}:00:00'
    print(f'wrfout_d01_2022-{pdtime.month:02d}-{pdtime.day:02d}_{pdtime.hour:02d}:00:00')
    try:
        fh      = xr.open_dataset(wrffile)
    except:
        print('Could not open the file',wrffile)
        sys.exit()
        
    lat,lon = fh.XLAT.data[0],fh.XLONG.data[0]

    stn_idx = np.argwhere((np.abs(lat-llat)<0.013)&(np.abs(lon-llon)<0.02))[0]
    if shift is not None:
        print('SHIFT',shift)
        if shift=='south':
            stn_idx = [stn_idx[0]-1,stn_idx[1]]
        if shift=='north':
            stn_idx = [stn_idx[0]+1,stn_idx[1]]
        if shift=='east':
            stn_idx = [stn_idx[0]-1,stn_idx[1]+1]
        if shift=='west':
            stn_idx = [stn_idx[0]-1,stn_idx[1]-1]
    
    timedict = {0:0,15:1,30:2,45:3}
    tidx = timedict[pdtime.minute]

    if 'thom' in run_dir:
        fh.close()

        # All of these variables are the lowest level of vertical coordinates - NOT surface variables
        ncfile = Dataset(wrffile,'r')
        if time =='all':
            temp = wrf.getvar(ncfile,'temp',units='K',timeidx=wrf.ALL_TIMES)[:,0,stn_idx[0],stn_idx[1]].data # kelvin
            pres = wrf.getvar(ncfile,'pres',units='Pa',timeidx=wrf.ALL_TIMES)[:,0,stn_idx[0],stn_idx[1]].data # Pa
            qv   = ncfile.variables['QVAPOR'][:,0,stn_idx[0],stn_idx[1]] # kg/kg
            qs   = ncfile.variables['QSNOW'][:,0,stn_idx[0],stn_idx[1]]
            if reflectivity is True:
                dbz = ncfile.variables['REFL_10CM'][:,0,stn_idx[0],stn_idx[1]]
        
        else:
            temp = wrf.getvar(ncfile,'temp',units='K',timeidx=tidx)[0,stn_idx[0],stn_idx[1]].data # kelvin
            pres = wrf.getvar(ncfile,'pres',units='Pa',timeidx=tidx)[0,stn_idx[0],stn_idx[1]].data # Pa
            qv   = ncfile.variables['QVAPOR'][tidx,0,stn_idx[0],stn_idx[1]] # kg/kg
            qs   = ncfile.variables['QSNOW'][tidx,0,stn_idx[0],stn_idx[1]]
            if reflectivity is True:
                dbz = ncfile.variables['REFL_10CM'][tidx,0,stn_idx[0],stn_idx[1]]
            
        ncfile.close()
        if reflectivity is True:
            return temp,pres,qv,qs,dbz
        else:
            return temp,pres,qv,qs

    if 'morr' in run_dir:
        if time=='all':
            ns = fh.QNSNOW.data[:,0,stn_idx[0],stn_idx[1]]
            qs = fh.QSNOW.data[:,0,stn_idx[0],stn_idx[1]]
            if reflectivity is True:
                dbz = fh.variables['REFL_10CM'][:,0,stn_idx[0],stn_idx[1]]
        else:
            ns = fh.QNSNOW.data[tidx,0,stn_idx[0],stn_idx[1]]
            qs = fh.QSNOW.data[tidx,0,stn_idx[0],stn_idx[1]]
            if reflectivity is True:
                dbz = fh.variables['REFL_10CM'][tidx,0,stn_idx[0],stn_idx[1]]

        fh.close()
        
        if reflectivity:
            return ns,qs,dbz
        else:
            return ns,qs
    
def find_moments(temp,pres,qv,qs):
    '''
    This function calculates the second and third moment of the particle size distribution as described in Thompson et al. 2008 and Field et al. 2005. These moments are needed to find the PSD for the Thompson scheme. 

    Input can come directly from the "get_wrf_for_psd" function. 
    '''
    # calculate second moment
    rho = (0.622*pres)/(287.04*temp*(0.622+qv))
    
    dens = density(pres * units.Pa, temp * units.kelvin, qv * units('kg/kg'))
#     print('rho,density',rho,dens)
    
#     qs = qs*10
    m2  = (1/0.069) * (qs*rho)

    tc = temp-273.15

    # Look-up Dictionaries
    al = {1:5.065339, 2:-0.062659, 3:-3.032362, 4:0.029469, 5:-0.000285, 6:0.312550, 
         7:0.000204, 8:0.003199, 9:0.000000, 10:-0.015952}
    bl = {1:0.476221, 2:-0.015896, 3:0.165977, 4:0.007468, 5:-0.000141, 6:0.060366, 
         7:0.000079, 8:0.000594, 9:0.000000, 10:-0.003577}

    # calculate third moment
    n = 3
    loga = al[1]+(al[2]*tc)+(al[3]*n)+(al[4]*tc*n)+(al[5]*(tc**2))+(al[6]*(n**2))+ \
            (al[7]*(tc**2)*n)+(al[8]*tc*(n**2))+(al[9]*(tc**3))+(al[10]*(n**3))
    a = 10**loga

    b = bl[1]+(bl[2]*tc)+(bl[3]*n)+(bl[4]*tc*n)+(bl[5]*(tc**2))+(bl[6]*(n**2))+ \
            (bl[7]*(tc**2)*n)+(bl[8]*tc*(n**2))+(bl[9]*(tc**3))+(bl[10]*(n**3))

    m3 = a*(m2**b)
    
    return m2,m3


def find_psd_thom(pdtime_range,station,run_dir,bins=[], reflectivity=False, shift=None):
    '''
    Calculate the snow-particle PSD as described in Thompson et al. 2008. 

    INPUT:
     - pdtime_range: a pandas datetime "date_range" for the time period of interest
     - station: a predetermined lat/lon location for the model grid 
     - run_dir: directory containing the WRF model output

    OUTPUT:
     - Ds: the bins of diameter size
     - nD: snow particle number concentration as a function of diameter
    '''
    
    # calculate diameter sizes (based on micro scheme) or use bins given
    if len(bins)==0:
        ## create bins of particle diameters (taken from Thompson microphysics scheme) ###
        nbins = 100
        D0s = 200e-6 # minimum snow diameter

        # snow particle bins from min diameter to 2cm
        Ds  = np.zeros(99)
        dts = np.copy(Ds)
        xdx = np.zeros(100)
        xdx[0] = D0s
        xdx[-1]= 0.02
        for ni in np.arange(1,nbins-1):
            xdx[ni] = np.exp(((ni-1)/nbins) * np.log(xdx[-1]/xdx[0]) + np.log(xdx[0]))
        for ni in np.arange(0,nbins-1):
            Ds[ni] = np.sqrt(xdx[ni]*xdx[ni+1])
            dts[ni]= xdx[ni+1] - xdx[ni]
    else:
        Ds = bins # bins must be in units of meters!
        
    ## Calculate PSD
    ND = np.zeros(((len(pdtime_range)*4),len(Ds)))
    if reflectivity:
        dbzs = []
    count=0
    for t,pdtime in enumerate(pdtime_range):
    
        if reflectivity is False:
            temp,pres,qv,qs = get_wrf_for_psd(pdtime,station,run_dir,time='all')
        elif reflectivity is True:
            temp,pres,qv,qs,dbz = get_wrf_for_psd(pdtime,station,run_dir,time='all',reflectivity=reflectivity,shift=shift)
            dbzs = np.hstack((dbzs,dbz))
            
        m2,m3 = find_moments(temp,pres,qv,qs)

        ## calculate number concentration as a function of diameter (N(D))
        for i in range(len(m2)):
            ND[count] = ((m2[i]**4)/(m3[i]**3)) * ( (490.6*np.exp((-m2[i]/m3[i])*20.78*Ds)) + \
                        ((17.46*((m2[i]/m3[i])*Ds)**0.6357)*np.exp((-m2[i]/m3[i])*3.29*Ds)))
            count = count+1
    
    if reflectivity is True:
        return Ds,ND,dbzs
    else:
        return Ds,ND # Diameter sizes, Number concentration

def find_psd_morr(pdtime_range,station,run_dir,bins=[],return_slope=False, reflectivity=False, shift=None):
    '''
    Calculate the snow-particle PSD as described in Morrison et al. 2009. 

    INPUT:
     - pdtime_range: a pandas datetime "date_range" for the time period of interest
     - station: a predetermined lat/lon location for the model grid
     - run_dir: directory containing the WRF model output

    OUTPUT:
     - Ds: the bins of diameter size
     - nD: snow particle number concentration as a function of diameter
    '''

    # Calculate diameter size bins, or use bins given
    if len(bins)==0:
        ## create bins of particle diameters (taken from Thompson microphysics scheme) ###
        nbins = 100
        D0s = 200e-6 # minimum snow diameter

        # snow particle bins from min diameter to 2cm
        Ds  = np.zeros(99)
        dts = np.copy(Ds)
        xdx = np.zeros(100)
        xdx[0] = D0s
        xdx[-1]= 0.02
        for ni in np.arange(1,nbins-1):
            xdx[ni] = np.exp(((ni-1)/nbins) * np.log(xdx[-1]/xdx[0]) + np.log(xdx[0]))
        for ni in np.arange(0,nbins-1):
            Ds[ni] = np.sqrt(xdx[ni]*xdx[ni+1])
            dts[ni]= xdx[ni+1] - xdx[ni]
    else:
        Ds = bins # bins must have units of meters!!
        
    ## Calculate PSD
    ND    = np.zeros(( (len(pdtime_range)*4),len(Ds) ))
    n0_all= np.zeros(( (len(pdtime_range)*4) ))
    slope_all = np.zeros(( (len(pdtime_range)*4) ))
    if reflectivity:
        dbzs = []
    count = 0
    for t,pdtime in enumerate(pdtime_range):
    
        if reflectivity is True:
            ns,qs,dbz = get_wrf_for_psd(pdtime,station,run_dir,time='all',reflectivity=True, shift=shift)
            dbzs = np.hstack((dbzs,dbz))
        else:
            ns,qs = get_wrf_for_psd(pdtime,station,run_dir,time='all')
#         print(ns.shape)

        # constants
        c = 100*(np.pi/6.)
        d = 3.0

        slope = ((c*ns*scipy.special.gamma(d+1))/qs)**(1./d)
        n_0   = ns * slope

        ## calculate number concentration as a function of diameter (N(D))
        for i in range(len(slope)):
            ND[count] = n_0[i]*np.exp(-slope[i]*Ds)
            n0_all[count] = n_0[i]
            slope_all[count] = slope[i]
            count = count+1
    
    # I belive ND is in units of m^-3 m^-1!!!
    if return_slope is True:
        if reflectivity is True:
            return ND,n0_all,slope_all,dbzs
        else:
            return ND, n0_all, slope_all 
    if return_slope is False:
        if reflectivity is True:
            return Ds,ND,dbzs
        else:
            return Ds,ND
    
def ob_n0_lambda(dsdperminute,bins,running_mean=False):
    # dsdperminute is the PIP PSD (or DSD) data that you want to get N0 and Lambda values for         
    if running_mean is False:
        dsdperminute_avg = np.empty([int(np.shape(dsdperminute)[0]/15),np.shape(dsdperminute)[1]])
        for i in range(int(np.shape(dsdperminute)[0]/15)):
            dsdperminute_avg[i] = np.nanmean(dsdperminute[int(i*15):int(i*15)+15])
            
    if running_mean is True:
        dsdperminute_avg = np.empty(np.shape(dsdperminute))
        for i in range(np.shape(dsdperminute)[0]):
            if i < 8:
                curr_avg = np.nanmean(dsdperminute[0:15, :], axis=0)
            elif i >= 8 and i <= np.shape(dsdperminute)[0] - 7:
                curr_avg = np.nanmean(dsdperminute[int(i-(15/2)):int(i+(15/2)), :], axis=0)
            elif i > np.shape(dsdperminute)[0] - 7:
                curr_avg = np.nanmean(dsdperminute[np.shape(dsdperminute)[0]-15:, :], axis=0)
            dsdperminute_avg[i] = curr_avg

    N0s = np.empty(np.shape(dsdperminute_avg)[0])
    Lambdas = np.empty(np.shape(dsdperminute_avg)[0])

    # for each of the 15-minute running window averaged PSDs, calculate the intercept paramter (N0) and exponential
    #     slope parameter (lambda) using the scipy.optimize curve_fit function
    for i in range(np.shape(dsdperminute_avg)[0]):
        try:
            if np.isnan(np.mean(dsdperminute_avg[i])):
                print('NaN')
                N0s[i] = np.nan
                Lambdas[i] = np.nan

            else:
                params = curve_fit(lambda x,N0,Lambda: N0*np.exp(-Lambda*x), bins, dsdperminute_avg[i], p0 = [1e4, 2], maxfev=600)

                # below is a check to make sure the calculated params aren't erroneous/unphysical
                # params[0][0] is N0, params[0][1] is Lambda
                if params[0][0] > 0 and params[0][0] < 10**7 and params[0][1] > 0 and params[0][1] < 10:
                    N0s[i] = params[0][0]
                    Lambdas[i] = params[0][1]
        except Exception as e:
            print(e)
            continue

    return N0s,Lambdas


def get_observed_lwe(station,start_time,end_time):
    
    stn_code = {'KMQT':'006','KAPX':'007'}
    
    pdtime_range = pd.date_range(start=start_time, end=end_time, freq='D')
    if end_time.day == start_time.day+1:
        print('This is an especially complicated time period')
        pdtime_range = pd.date_range(start=start_time, end=(end_time+pd.Timedelta('1D')), freq='D')
        
    nrr = []
    hours = []
    flag = False
    for i,pdt in enumerate(pdtime_range):
        if station=='KBUF':
            filename = f'/nfs/turbo/seas-hutsona/les_data/PARSIVEL/{pdt.strftime("%Y%m%d")}_BUF.nc'
        else:
            filename = f'/nfs/turbo/seas-hutsona/les_data/PIP/{start_time.year}_{station[1:]}/adjusted_edensity_lwe_rate/{stn_code[station]}{pdt.strftime("%Y%m%d")}_min.nc'
        try:
            fh = xr.open_dataset(filename)
            piptime   = pd.to_datetime(fh.time.data)
            if station=='KBUF':
                nrrs_full = fh.Heymsfield_preciprate_total.data
#                 nrrs_full = fh.Manual_preciprate_total.data
#                 nrrs_full[nrrs_full<0.0] = np.nan
            else:
                nrrs_full = fh.nrr_adj.data
            obtime  = np.copy(piptime)
            if (piptime[0]<start_time):
                ot_start = np.argwhere(piptime==start_time)[0][0]
            else:
                ot_start = 0

            if (piptime[-1]>end_time):
                ot_end = np.argwhere(obtime==end_time)[0][0]+1
            else:
                ot_end = np.shape(obtime)[0]

    #         print(ot_start,ot_end,'<<Test')
            obtime= obtime[ot_start:ot_end]
            nrrs  = nrrs_full[ot_start:ot_end]

            ### calculate hourly mean
            special_flag = False
            if (obtime[-1]-obtime[0])<pd.Timedelta('1H'):
#                 if (obtime[-1]==obtime[0]):
#                     special_flag=False
#                 else:
#                 print('This 1 hour period straddles two PIP files')
                special_flag = True

            hour_flag = False
            for h,hour in enumerate(pd.date_range(start=obtime[0],end=obtime[-1],freq='H')):
                if flag is False:
                    nrr.append(np.nan)
                    flag = True
                    if hour.hour==23:
                        index   = np.argwhere(obtime==hour)[0][0]
                        next_mean = np.nanmean(nrrs[index:])

                elif (hour.hour==0):
                    if special_flag ==True:
                        index_end = np.argwhere(obtime==obtime[-1])[0][0]
                        nnext_mean = np.nansum(np.asarray([next_mean,np.nanmean(nrrs[:index_end])]))
                        nrr.append(nnext_mean)
                        hours.append(pd.to_datetime(obtime[-1]))
                        hour_flag = True
                    else:
                        nrr.append(next_mean)

                else:
                    index   = np.argwhere(obtime==hour)[0][0]
                    nrr.append(np.nanmean(nrrs[index-60:index]))
                    if hour.hour==23:
                        print('last hour of the day')
                        next_mean = np.nanmean(nrrs[index:])

                if hour_flag == False:
                    hours.append(hour)
            
        except:
            print(filename,'does not exist...returning NaNs for this day')
            for h,hour in enumerate(pd.date_range(start=pdt,end=pdt+pd.Timedelta('23H'),freq='H')):
                hours.append(hour)
                nrr.append(np.nan)

    return np.array(nrr),hours

def get_observed_density(var,station,start_time,end_time,hourly=False):
    stn_code = {'KMQT':'006','KAPX':'007'}
    
    pdtime_range = pd.date_range(start=start_time, end=end_time, freq='D')
    denss = []
    obtime_all = []
    for i,pdt in enumerate(pdtime_range):
        if station=='KBUF':
            filename = f'/nfs/turbo/seas-hutsona/les_data/PARSIVEL/{pdt.strftime("%Y%m%d")}_BUF.nc'
        else:
            if var=='rho':
                filename = f'/nfs/turbo/seas-hutsona/les_data/PIP/{start_time.year}_{station[1:]}/edensity_distributions/{stn_code[station]}{pdt.strftime("%Y%m%d")}_rho.nc'
            elif var=='ed_adj':
                filename = f'/nfs/turbo/seas-hutsona/les_data/PIP/{start_time.year}_{station[1:]}/adjusted_edensity_lwe_rate/{stn_code[station]}{pdt.strftime("%Y%m%d")}_min.nc'
        
        try:
            fh = xr.open_dataset(filename)
        except:
            print('PIP file does not exist...appending NaNs')
            if hourly is True:
                for hour in pd.date_range(start=pdt, end=pdt+pd.Timedelta(hours=23),freq='1H'):
                    denss.append(np.nan)
                    obtime_all.append(hour)
            if hourly is False:
                for hour in pd.date_range(start=pdt, end=pdt+pd.Timedelta(hours=24),freq='15min')[:-1]:
                    denss.append(np.nan)
                    obtime_all.append(hour)
            continue
            
        piptime = pd.to_datetime(fh.time.data)
        if station=='KBUF':
            var = 'Heymsfield_effective_density'
        dens    = fh[var].data
        obtime  = np.copy(piptime)
        if (piptime[0]<start_time):
            obtime= obtime[np.argwhere(piptime==start_time)[0][0]:]
            dens = dens[np.argwhere(piptime==start_time)[0][0]:]
        if (piptime[-1]>end_time):
            obtime= obtime[:np.argwhere(obtime==end_time)[0][0]+1]
            dens = dens[:np.argwhere(obtime==end_time)[0][0]+1]
            
        if hourly is False:
            index_prev = 0
            for hour in pd.date_range(start=obtime[0],end=obtime[-1],freq='15min')[1:]:
                obtime_all.append(hour)
                index = np.argwhere(obtime==hour)[0][0]
                denss.append(np.nanmean(dens[index_prev:index]))
                index_prev = index
                if (hour.hour==23)&(hour.min==45):
                    denss.append(np.nanmean(dens[index_prev:]))
                
        if hourly is True:
            index_prev = 0
    #         print(pd.date_range(start=obtime[0],end=obtime[-1],freq='15min'))
            for hour in pd.date_range(start=obtime[0],end=obtime[-1],freq='1H')[1:]:
                obtime_all.append(hour)
                index = np.argwhere(obtime==hour)[0][0]
                denss.append(np.nanmean(dens[index_prev:index]))
                index_prev = index
                if (hour.hour==23):
                    obtime_all.append(hour+pd.Timedelta(hours=1))
                    denss.append(np.nanmean(dens[index_prev:]))
#         print(pdt,np.array(denss).shape,np.shape(obtime_all))
                    
    return np.array(denss),obtime_all

