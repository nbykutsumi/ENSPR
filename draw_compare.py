#! /opt/local/bin/python
import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset
from math import sin, cos, sqrt, atan2, radians
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import h5py
import sys
import os



"""
prj     = 'GPM.KuPR'
prdLv   = 'L2'
#prdVer  = '03'
prdVer  = '05'
var     = 'NS/SLV/precipRateESurface'

#prj     = 'GPM.KaPR'
#prdLv   = 'L2'
#prdVer  = '03'
##var     = 'HS/SLV/precipRateESurface'
#var     = 'MS/SLV/precipRateESurface'

#prj     = 'GPM.GMI'
#prdLv   = 'L2'
#prdVer  = '05'
#var     = 'S1/surfacePrecipitation'
"""


def fig_precip(Lat, Lon, Dat, stitle, figPath, vmin, vmax):
    figmap = plt.figure(figsize=(3.5,3.5),dpi=80)
    axmap  = figmap.add_axes([0.1,0.12,0.8,0.75])
    
    dlat  = 10
    dlon  = 10
    lllat = clat-dlat
    urlat = clat+dlat
    lllon = clon-dlon
    urlon = clon+dlon
    
    M  = Basemap(resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
    
    Dat = ma.masked_less(Dat,0)
    im = M.contourf(Lon, Lat, Dat, extend="both", cmap="jet", vmin=0, vmax=10)
    
    print precip1.min(), precip1.max()
    print precip1
    # meridians and parallels
    dlatlon = 3
    parallels = np.arange(-90, 90, dlatlon)
    meridians = np.arange(0, 360, dlatlon)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt="%d")
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, rotation=60, fmt="%d")
    
    # coastlines
    M.drawcoastlines()
    
    # title
    plt.title(stitle)
    
    # colorbar
    plt.colorbar(im, orientation="horizontal") 
    
    # save
    plt.savefig(figPath)
    print figPath

    


#  Map of all four retrievals from EPC

# Retrieved file
ncfile = "/home/utsumi/bin/ENSPR/GPM_EPC_010787_20160121_2335.NS_MS.nc"
f = Dataset(ncfile)
for attr in f.ncattrs():
    print attr, '=', getattr(f, attr)

Lat = f['latitude'][:]
Lon = f['longitude'][:]
emis= f['emis'][:]
tb= f['Tb'][:]

precip1= f['NS/precip'][:]
precip2= f['NS/precip2'][:]
precip3= f['MS/precip'][:]
precip4= f['MS/precip2'][:]
#precipg= 1.0*np.array(f['GPROF/precip'])
#sfc= np.array(f['GPROF/sfc_class'])

clat = float(getattr(f,'center_lat'))
clon = float(getattr(f,'center_lon'))
cdate = getattr(f,'center_date')
#nsms = getattr(f,'NS_or_MS')
nsms = 'NS'
f.close()

nx = Lat.shape[0]
ny = Lon.shape[1]
print nx,ny

cdate2 = cdate.split()[0]
ctime2 = cdate.split()[1].split(':')
cdate = cdate2 + ' ' + ctime2[0] + ctime2[1] + ' UTC' 
print cdate


# GPROF
ncfile = "/home/utsumi/bin/ENSPR/GPM_EPC_010787_20160121_2335.NS_MS.nc"
f = Dataset(ncfile)
for attr in f.ncattrs():
    print attr, '=', getattr(f, attr)

Lat = f['latitude'][:]
Lon = f['longitude'][:]
emis= f['emis'][:]
tb= f['Tb'][:]

precip1= f['NS/precip'][:]
precip2= f['NS/precip2'][:]
precip3= f['MS/precip'][:]
precip4= f['MS/precip2'][:]
#precipg= 1.0*np.array(f['GPROF/precip'])
#sfc= np.array(f['GPROF/sfc_class'])

clat = float(getattr(f,'center_lat'))
clon = float(getattr(f,'center_lon'))
cdate = getattr(f,'center_date')
#nsms = getattr(f,'NS_or_MS')
nsms = 'NS'
f.close()

nx = Lat.shape[0]
ny = Lon.shape[1]
print nx,ny



# create figures




# create polar stereographic Basemap instance.
# 2400-km per side is OK for 5-min GMI
#m = Basemap(width=3600000,height=3600000,
#m = Basemap(width=1400000,height=1400000,
#            resolution='l',projection='stere',\
#            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(width=6000000,height=4000000,
#            resolution='l',projection='stere',\
#            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(projection='cyl',llcrnrlat=27,urcrnrlat=35,\
#            llcrnrlon=-101,urcrnrlon=-93,resolution='i')




#ncfile = sys.argv[2]
#f = h5py.File(ncfile,'r')
#lats_cmb_NS = 1.0*np.array(f['NS/Latitude'])
#lons_cmb_NS = 1.0*np.array(f['NS/Longitude'])
#precip_cmb_NS = 1.0*np.array(f['NS/surfPrecipTotRate'])
#lats_cmb_MS = 1.0*np.array(f['MS/Latitude'])
#lons_cmb_MS = 1.0*np.array(f['MS/Longitude'])
#precip_cmb_MS = 1.0*np.array(f['MS/surfPrecipTotRate'])
#f.close()
#
#ncfile = sys.argv[3]
#f = h5py.File(ncfile,'r')
#lats_dpr_NS = 1.0*np.array(f['NS/Latitude'])
#lons_dpr_NS = 1.0*np.array(f['NS/Longitude'])
#precip_dpr_NS = 1.0*np.array(f['NS/SLV/precipRateESurface'])
#lats_dpr_MS = 1.0*np.array(f['MS/Latitude'])
#lons_dpr_MS = 1.0*np.array(f['MS/Longitude'])
#precip_dpr_MS = 1.0*np.array(f['MS/SLV/precipRateESurface'])
#f.close()
#
#
#
##---------------------------------------------
#
#
#nx = tb.shape[0]
#ny = tb.shape[1]
#nc = tb.shape[2]
#print "GMI dims=",nx,ny,nc
# 
#tb89h = tb[:,:,8]
#tb89h[tb89h < 0]= np.nan
#
#precip1[precip1 < 0.1]= np.nan
#precip2[precip2 < 0.1]= np.nan
#precip3[precip3 < 0.1]= np.nan
#precip4[precip4 < 0.1]= np.nan
#precipg[precipg < 0.1]= np.nan
#precip_cmb_NS[precip_cmb_NS < 0.1]= np.nan
#precip_cmb_MS[precip_cmb_MS < 0.1]= np.nan
#precip_dpr_NS[precip_dpr_NS < 0.1]= np.nan
#precip_dpr_MS[precip_dpr_MS < 0.1]= np.nan
#
#e10h = emis[:,:,1]
#ts = 1.0*emis[:,:,9]  
#tqv = 1.0*emis[:,:,10] 
#e10h[e10h < 0]= np.nan
#ts[ts < 0]= np.nan
#tqv[tqv < 0]= np.nan
#
#edge_lat1 = lats[:,0]
#edge_lon1 = lons[:,0]
#edge_lat2 = lats[:,ny-1]
#edge_lon2 = lons[:,ny-1]
#edge_lat3 = lats[0,:]
#edge_lon3 = lons[0,:]
#edge_lat4 = lats[nx-1,:]
#edge_lon4 = lons[nx-1,:]
#
#nx2 = lats_cmb_NS.shape[0]
#ny2 = lats_cmb_NS.shape[1]
#print "DPR NS dims=",nx2,ny2
#edge2_lat1 = lats_cmb_NS[:,0]
#edge2_lon1 = lons_cmb_NS[:,0]
#edge2_lat2 = lats_cmb_NS[:,ny2-1]
#edge2_lon2 = lons_cmb_NS[:,ny2-1]
#
#nx3 = lats_cmb_MS.shape[0]
#ny3 = lats_cmb_MS.shape[1]
#print "DPR MS dims=",nx3,ny3
#edge3_lat1 = lats_cmb_MS[:,0]
#edge3_lon1 = lons_cmb_MS[:,0]
#edge3_lat2 = lats_cmb_MS[:,ny3-1]
#edge3_lon2 = lons_cmb_MS[:,ny3-1]
#
#
## create figure and axes instances
#fig = plt.figure(figsize=(21,15),dpi=80)
#
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
## create polar stereographic Basemap instance.
## 2400-km per side is OK for 5-min GMI
##m = Basemap(width=3600000,height=3600000,
#m = Basemap(width=1400000,height=1400000,
#            resolution='l',projection='stere',\
#            lat_ts=clat,lat_0=clat,lon_0=clon)
##m = Basemap(width=6000000,height=4000000,
##            resolution='l',projection='stere',\
##            lat_ts=clat,lat_0=clat,lon_0=clon)
##m = Basemap(projection='cyl',llcrnrlat=27,urcrnrlat=35,\
##            llcrnrlon=-101,urcrnrlon=-93,resolution='i')
#
#
#x, y = m(lons, lats)     # EPC map proj coordinates.
#xn, yn = m(lons_cmb_NS, lats_cmb_NS) # NS map proj coordinates.
#xm, ym = m(lons_cmb_MS, lats_cmb_MS) # MS map proj coordinates.
#
#x1, y1 = m(edge_lon1, edge_lat1) 
#x2, y2 = m(edge_lon2, edge_lat2) 
#x3, y3 = m(edge_lon3, edge_lat3) 
#x4, y4 = m(edge_lon4, edge_lat4) 
#
#dx1, dy1 = m(edge2_lon1, edge2_lat1) 
#dx2, dy2 = m(edge2_lon2, edge2_lat2) 
#
#dx3, dy3 = m(edge3_lon1, edge3_lat1) 
#dx4, dy4 = m(edge3_lon2, edge3_lat2) 
#
#dlat = 2
#dlon = 3
#
#
#ax=plt.subplot(3, 4, 1)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'EPC Precip (CMB-NS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
##clevs = np.arange(0,15,1)
##clevs2 = np.arange(0,15,1)
##levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#cs = m.contourf(x,y,precip1,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx1,dy1,color='black',linewidth=1.0)
#plt.plot(dx2,dy2,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 1"
#
#
#ax=plt.subplot(3, 4, 2)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'EPC Precip (DPR-NS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(x,y,precip2,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx1,dy1,color='black',linewidth=1.0)
#plt.plot(dx2,dy2,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 2"
#
#
#ax=plt.subplot(3, 4, 3)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'EPC Precip (CMB-MS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(x,y,precip3,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 3"
#
#
#
#ax=plt.subplot(3, 4, 4)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'EPC Precip (DPR-MS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(x,y,precip4,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 4"
#
#
#ax=plt.subplot(3, 4, 5)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'Precip (CMB-NS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(xn,yn,precip_cmb_NS,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx1,dy1,color='black',linewidth=1.0)
#plt.plot(dx2,dy2,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 5"
#
#
#
#
#ax=plt.subplot(3, 4, 6)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'Precip (DPR-NS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(xn,yn,precip_dpr_NS,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx1,dy1,color='black',linewidth=1.0)
#plt.plot(dx2,dy2,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 6"
#
#
#ax=plt.subplot(3, 4, 7)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'Precip (CMB-MS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(xm,ym,precip_cmb_MS,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 7"
#
#
#
#ax=plt.subplot(3, 4, 8)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'Precip (DPR-MS)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(xm,ym,precip_dpr_MS,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 8"
#
#
#ax=plt.subplot(3, 4, 9)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'Precip (GPROF-V5)'
#clevs = np.arange(0,25.0,0.2)
#clevs2 = np.arange(0,25.0,5.0)
#cs = m.contourf(x,y,precipg,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 9"
#
#
#ax=plt.subplot(3, 4, 10)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90,90,dlat)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0,360,dlon)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#
#svar = 'GMI 89H (K)'
#clevs = np.arange(160,290,1)
#clevs2 = np.arange(160,290,20)
#cs = m.contourf(x,y,tb89h,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=1.5)
#plt.plot(x2,y2,color='black',linewidth=1.5)
#plt.plot(x3,y3,color='black',linewidth=1.5)
#plt.plot(x4,y4,color='black',linewidth=1.5)
#plt.plot(dx3,dy3,color='black',linewidth=1.0)
#plt.plot(dx4,dy4,color='black',linewidth=1.0)
#
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 10"
#
#
#
#
#
#
#
#plt.suptitle(cdate, fontsize=16, fontweight='bold')
#
#
#
#plt.show()
##plt.savefig(outfile, dpi=80, transparent='True', bbox_inches='tight', pad_inches=0.1)
#sys.exit()
#
#
#
#
#
#
#
#
#
#ax=plt.subplot(3, 3, 8)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#parallels = np.arange(-90.,90,2.)
#m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
#meridians = np.arange(0.,360.,2.)
#m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
#m.drawlsmask(land_color = 'lightgray')
#
#svar = 'Ku PIA (using CMB Ku)'
#clevs = np.arange(0,15.0,0.2)
#clevs2 = np.arange(0,15.0,5.0)
##cs = m.contourf(x,y,precipg,clevs,extend='both')
#cs = m.contourf(x,y,pia,clevs,extend='both',cmap='jet')
#plt.plot(x1,y1,color='black',linewidth=0.5)
#plt.plot(x2,y2,color='black',linewidth=0.5)
#cbar = m.colorbar(cs,location='bottom',pad="8%")
#cbar.set_ticks(clevs2)
##cbar.set_label(svar + ' (K)', fontsize=14, family='Arial', fontweight='bold')
#plt.title(svar, fontsize=14, fontweight='bold')
#
#print "frame 8"
#
#
#
#
#
#
##plt.show()
##plt.savefig(ncfile + '.tb.png', transparent='True', bbox_inches='tight', pad_inches=0.25)
#plt.savefig('gmi_precip.png', transparent='True', bbox_inches='tight', pad_inches=0.1)
##print imgname
#
#sys.exit()
