import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from PIL import Image

import h5py
import glob
from numpy import *
from datetime import datetime, timedelta
from bisect import bisect_right
import JPLDB

"""
Year = 2014
Mon  = 9
Day  = 6
Hour = 20
Mnt  = 34
lat,lon = 15.1192, 8.47616
granule = 2974
"""
"""
Year = 2014
Mon  = 9
Day  = 10
Hour = 1
Mnt  = 33
lat,lon = 40.8904, -96.5064
granule = 3024
"""
"""
Year = 2014
Mon  = 9
Day  = 19
Hour = 8
Mnt  = 30
lat,lon = -24.621, -50.3465
granule = 3169
"""
Year = 2014
Mon  = 9
Day  = 30
Hour = 3
Mnt  = 50
lat,lon = 9.40949, -11.0404
granule = 3337



dlatlon = 0.8
lllat = lat-dlatlon
lllon = lon-dlatlon
urlat = lat+dlatlon
urlon = lon+dlatlon
BBox  = [[lllat,lllon],[urlat,urlon]]



jpl  = JPLDB.JPLDB()

def single_fig_KuPR(Dat,Lat,Lon, BBox, figPath, cbarPath, stitle, scatterFlag=True, DPR=None):

    #-- figure -
    figmap = plt.figure(figsize=(3.5,3.5))
    axmap  = figmap.add_axes([0.18, 0.12, 0.75, 0.75])
    [[lllat,lllon],[urlat,urlon]] = BBox
    M      = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

    #- plot -
    vmin, vmax=0, 30
    if scatterFlag==True:
        im = M.scatter(Lon,Lat,c=Dat, cmap="jet", vmin=vmin, vmax=vmax,s=20)
    else:
        clevels = np.arange(0,30+1,1)
        im = M.contourf(Lon,Lat,Dat, extend="max",cmap="jet", vmin=vmin, vmax=vmax)

    #- coastlines
    M.drawcoastlines()


    #-- draw center cross
    M.plot((lon,lon), (lat-30,lat+30), linewidth=0.5, color="k")
    M.plot((lon-30,lon+30), (lat,lat), linewidth=0.5, color="k")

    #-- meridians and parallels -
    if urlat-lllat > 5:
        dlatlon=3
    else:
        dlatlon   = 0.2

    parallels = arange(-90, 90, dlatlon)
    meridians = arange(0, 360, dlatlon)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8,linewidth=0.5, fmt="%.1f")
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8,linewidth=0.5,rotation=60, fmt="%.1f")

    #-- title --
    plt.title(stitle, fontsize=10)
    #-- save
    figmap.savefig(figPath)

    print figPath
    #-- colorbar -
    if cbarPath != None:
        figcbar = plt.figure(figsize=(3.5,0.5))
        axcbar   = figcbar.add_axes([0.05, 0.5, 0.9, 0.4])
        plt.colorbar(im, extend="max", cax = axcbar, orientation="horizontal")
        figcbar.savefig(cbarPath)
        print cbarPath


#-- original KuPR L2 data ---

kuDir  = "/work/a01/utsumi/data/GPM/GPM.KuPR/L2/05/%04d/%02d"%(Year,Mon)
print kuDir
kuPath  = glob.glob(kuDir + "/GPMCOR_KUR_*_*_%06d_*.h5"%(granule))[0]
print kuPath

varName = 'NS/SLV/precipRateESurface'
h5 = h5py.File(kuPath, "r")
Dat = h5[varName][:]
Lat = h5["NS/Latitude"][:]
Lon = h5["NS/Longitude"][:]

h5Grp   = "NS"
aYear    = h5['%s/ScanTime/Year'%h5Grp        ][:].astype('int')
aMonth   = h5['%s/ScanTime/Month'%h5Grp       ][:].astype('int')
aDay     = h5['%s/ScanTime/DayOfMonth'%h5Grp  ][:].astype('int')
aHour    = h5['%s/ScanTime/Hour'%h5Grp        ][:].astype('int')
aMinute  = h5['%s/ScanTime/Minute'%h5Grp      ][:].astype('int')
aSecond  = h5['%s/ScanTime/Second'%h5Grp      ][:].astype('int')
aMicSec  = h5['%s/ScanTime/MilliSecond'%h5Grp ][:].astype('int')*1000

aTime =  array( [aYear, aMonth, aDay, aHour, aMinute, aSecond, aMicSec] ).T
Time  = [datetime(year,month,day,hour,minute,second,micsec) for (year,month,day,hour,minute,second,micsec) in aTime]

figDir   = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
figPath  = figDir + "/temp.GPM.Ku.org.%05d.png"%(granule)
cbarPath = figDir + "/temp.cbar.png"
stitle   = "gra=%06d   %04d/%02d/%02d %02d:%02d\n"%(granule, Year,Mon,Day,Hour,Mnt)
stitle   = stitle + "KuPR(Org) lat,lon=(%.2f, %.2f)"%(lat, lon )
scatterFlag = True
single_fig_KuPR(Dat,Lat,Lon, BBox, figPath, cbarPath, stitle, scatterFlag=True, DPR=None)


#-- Ku from JPLDB --
jplDir =  "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5"
jplPath = glob.glob(jplDir + "/%04d%02d*%06d.strm"%(Year,Mon,granule))[0]

jpl.set_file(jplPath)
Dat = jpl.get_var("precip_NS")
gLat = jpl.get_var("glat1")
gLon = jpl.get_var("glon1")

sLat = jpl.get_var("slat")
sLon = jpl.get_var("slon")
timediff = jpl.get_var("timediff")

# based on gLat, gLon

figDir   = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
figPath  = figDir + "/temp.JPL.Ku.glatlon.%05d.png"%(granule)
cbarPath = figDir + "/temp.cbar.png"
stitle   = "gra=%06d %04d/%02d/%02d %02d:%02d\n KuPR(JPLDB) based on glat1 & glon1"%(granule, Year,Mon,Day,Hour,Mnt)
scatterFlag = True
single_fig_KuPR(Dat,gLat,gLon, BBox, figPath, cbarPath, stitle, scatterFlag=True, DPR=None)

