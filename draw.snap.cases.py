import matplotlib
matplotlib.use("Agg")
from numpy import *
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import myfunc.util as util
import myfunc.IO.MERGIR     as MERGIR
import os, sys
from bisect import bisect_right
from mpl_toolkits.basemap import Basemap
from PIL import Image
import glob
from myfunc.IO import GPyM
import Image

#-- KuPR -
prj     = 'GPM.KuPR'
prdLv   = 'L2'
prdVer  = '05'

kupr  = GPyM.GPM(prj, prdLv, prdVer)

#-- GMI-TB -
prj     = 'GPM.GMI'
prdLv   = 'L1B'
prdVer  = '05'

tb      = GPyM.GPM(prj, prdLv, prdVer)
#-----------
iYM   = [2014,9]
eYM   = [2014,9]
lYM   = util.ret_lYM(iYM, eYM)


ir    = MERGIR.MERGIR()
LatIR = ir.Lat
LonIR = ir.Lon
missIR= 75

LatBndIR  = (LatIR[1:]+LatIR[:-1])*0.5
LatBndIR  = r_[-60, LatBndIR, 60]

LonBndIR  = (LonIR[1:]+LonIR[:-1])*0.5
LonBndIR  = r_[-180, LonBndIR, 180]


#----------------------------
def load_img(figPath):
    iimg = Image.open(figPath)
    return asarray( iimg )
#----------------------------
def latlon2yx(lat,lon, LatBnd, LonBnd):
    y = bisect_right(LatBnd, lat)-1
    x = bisect_right(LonBnd, lon)-1
    return y,x

#----------------------------
def single_fig_IR(DTime, lat, lon, BBox, figPath, cbarPath):
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour
    Mnt     = DTime.minute

    #a2ir    = ma.masked_equal(ir.load_30min(DTime), missIR)
    a2ir    = ir.load_30min(DTime)
    [[lllat,lllon],[urlat,urlon]] = BBox
    y0      = bisect_right(LatBndIR, lllat)-1
    y1      = bisect_right(LatBndIR, urlat)-1
    x0      = bisect_right(LonBndIR, lllon)-1
    x1      = bisect_right(LonBndIR, urlon)-1
    
    a2in    = a2ir[y0:y1+1, x0:x1+1]

    LatMap  = LatBndIR[y0:y1+1]
    LonMap  = LonBndIR[x0:x1+1]
    LON,LAT = meshgrid(LonMap, LatMap)

    #-- figure -
    figmap = plt.figure(figsize=(3.5,3.5))
    axmap  = figmap.add_axes([0.18,0.18,0.75,0.75])
    M   = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

    #-- draw -
    cmap="gnuplot2_r"
    im = M.pcolormesh(LON,LAT,a2in, vmin=190,vmax=300,cmap=cmap)

    #-- draw center cross
    M.plot((lon,lon), (lat-30,lat+30), linewidth=1, color="w")
    M.plot((lon-30,lon+30), (lat,lat), linewidth=1, color="w")

    #-- coastline -
    M.drawcoastlines()

    #-- meridians and parallels -
    if urlat-lllat > 5:
        dlatlon=3
    else:
        dlatlon   = 1.0

    parallels = arange(-90, 90, dlatlon)
    meridians = arange(0, 360, dlatlon)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8,linewidth=0.5, fmt="%d")
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8,linewidth=0.5,rotation=60, fmt="%d")

    #-- title -
    stitle = "idx=%d IR %04d/%02d/%02d %02d:%02d"%(idx, Year,Mon,Day,Hour,Mnt)
    plt.title(stitle, fontsize=10)
    #-- save -
    plt.savefig(figPath)
    print figPath

    #-- colorbar -
    if cbarPath != None:
        figcbar = plt.figure(figsize=(3.5,0.5))
        axcbar   = figcbar.add_axes([0.05, 0.5, 0.9, 0.4])
        plt.colorbar(im, extend="both", cax = axcbar, orientation="horizontal")
        figcbar.savefig(cbarPath)
        print cbarPath

#----------------------------
def single_fig_KuPR(DTime, lat, lon, BBox, figPath, cbarPath, scatterFlag=True, DPR=None):
    varKu   = 'NS/SLV/precipRateESurface'
    dDTime  = timedelta(seconds=60*30)
    kuprObt = kupr(varKu, DTime-dDTime, DTime+dDTime, BBox=BBox)
    Dat     = ma.masked_less_equal(kuprObt.data,0)
    Lat     = kuprObt.lat
    Lon     = kuprObt.lon
    
    #-- figure -
    figmap = plt.figure(figsize=(3.5,3.5))
    axmap  = figmap.add_axes([0.18, 0.18, 0.75, 0.75])
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
        dlatlon   = 1.0

    parallels = arange(-90, 90, dlatlon)
    meridians = arange(0, 360, dlatlon)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8,linewidth=0.5, fmt="%d")
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8,linewidth=0.5,rotation=60, fmt="%d")
    
    #-- title --
    stitle = "KuPR %04d/%02d/%02d %02d:%02d NS=%.1f Dist=%.2f"%(Year,Mon,Day,Hour,Mnt, DPR, dist)
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
    
#----------------------------
def single_fig_TB(DTime, lat, lon, BBox, figPath, cbarPath, chName, scatterFlag=True):
    #S1: 10V 10H 19V 19H 23V 37V 37H 89V 89H
    dich    = {"10V":0, "10H":1, "19V":2, "19H":3, "23V":4, "37V":5, "37H":6, "89V":7, "89H":8}
    ich     = dich[chName]
    dDTime  = timedelta(seconds=60*30)
    tbObt   = tb(vartb, DTime-dDTime, DTime+dDTime, BBox=BBox)
    Dat     = tbObt.data[:,:,ich]
    Lat     = tbObt.lat
    Lon     = tbObt.lon
    
    #-- figure -
    figmap = plt.figure(figsize=(3.5,3.5))
    axmap  = figmap.add_axes([0.18, 0.18, 0.75, 0.75])
    [[lllat,lllon],[urlat,urlon]] = BBox
    M      = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
    
    #- plot -
    if chName in ["10V","10H"]:
        vmin, vmax = 80,150
    elif chName in ["89V","89H"]:
        vmin, vmax = 120, 300
    if scatterFlag==True:
        im = M.scatter(Lon,Lat,c=Dat, cmap="jet",s=28, vmin=vmin, vmax=vmax)
        #im = M.scatter(Lon,Lat,c=Dat, cmap="jet",s=28)
    else:
        clevels = np.arange(120,300+1,5)
        im = M.contourf(Lon,Lat,Dat, extend="both",cmap="jet", vmin=vmin, vmax=vmax)
 
    
    #- coastlines
    M.drawcoastlines()
    
    
    #-- draw center cross
    M.plot((lon,lon), (lat-30,lat+30), linewidth=0.5, color="k")
    M.plot((lon-30,lon+30), (lat,lat), linewidth=0.5, color="k")
    
    #-- meridians and parallels -
    if urlat-lllat > 5:
        dlatlon=3
    else:
        dlatlon   = 1.0

    parallels = arange(-90, 90, dlatlon)
    meridians = arange(0, 360, dlatlon)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8,linewidth=0.5, fmt="%d")
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8,linewidth=0.5,rotation=60, fmt="%d")
    
    #-- title --
    stitle = "TB(%s) %04d/%02d/%02d %02d:%02d +-30min"%(chName, Year,Mon,Day,Hour,Mnt)
    plt.title(stitle, fontsize=10)
    #-- save
    figmap.savefig(figPath)
    
    print figPath
    #-- colorbar -
    if cbarPath != None:
        figcbar = plt.figure(figsize=(3.5,0.5))
        axcbar   = figcbar.add_axes([0.05, 0.5, 0.9, 0.4])
        plt.colorbar(im, extend="both", cax = axcbar, orientation="horizontal")
        figcbar.savefig(cbarPath)
        print cbarPath
    


#----------------------------
for YM in lYM:
    Year,Mon = YM
    csvDir   = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    filtPath = csvDir + "/cases.filt.%04d%02d.csv"%(Year,Mon)
    f = open(filtPath); lines1=f.readlines(); f.close()
    for line1 in lines1:
        line1   = line1.strip().split(",")
        granule = int(line1[0])
        idx_tgt = int(line1[1])
        #for line in lines[1+11:]:
        casePath = glob.glob(csvDir + "/cmp.case.db.*.g.%06d.idx.%05d.csv"%(granule, idx_tgt))[0]
        f = open(casePath,"r"); lines2=f.readlines(); f.close() 

        caseName = casePath.split("/")[-1]
        idx_db   = int(caseName.split(".")[3])
        #for line2 in lines2[1:]:
        for line2 in lines2[1:2]:
            line2 = line2.strip().split(",")
            idx  = int(line2[0])
            KuNS = float(line2[1])
            dist = float(line2[3])
            Year = int(line2[4])
            Mon  = int(line2[5])
            Day  = int(line2[6])
            Hour = int(line2[7])
            Mnt  = int(line2[8])
            lat  = float(line2[10])
            lon  = float(line2[11])
        
            MntIR    = int(Mnt/30)*30
            DTime    = datetime(Year,Mon,Day,Hour,Mnt)
            DTimeIR  = datetime(Year,Mon,Day,Hour,MntIR)
            dlatlon  = 10
            BBoxWide     = [[lat-dlatlon, lon-dlatlon],[lat+dlatlon, lon+dlatlon]]
            
            dlatlon  = 0.7
            BBoxZoom     = [[lat-dlatlon, lon-dlatlon],[lat+dlatlon, lon+dlatlon]]
            
            figDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
            #figDir  = "/home/utsumi/fig"
            
            #- IR --
            figPath = figDir + "/IR.wide.%d.png"%(idx)
            cbarPath = figDir + "/cbar.IR.png"
            single_fig_IR(DTimeIR, lat, lon, BBoxWide, figPath, cbarPath)
            
            figPath = figDir + "/IR.zoom.%d.png"%(idx)
            cbarPath = figDir + "/cbar.IR.png"
            single_fig_IR(DTimeIR, lat, lon, BBoxZoom, figPath, cbarPath)
            
            #- KuPR --
            figPath = figDir + "/KuPR.wide.%d.png"%(idx)
            cbarPath = figDir + "/cbar.KuPR.png"
            single_fig_KuPR(DTime, lat, lon, BBoxWide, figPath, cbarPath, scatterFlag=True, DPR=KuNS)
            
            figPath = figDir + "/KuPR.zoom.%d.png"%(idx)
            cbarPath = figDir + "/cbar.KuPR.png"
            single_fig_KuPR(DTime, lat, lon, BBoxZoom, figPath, cbarPath, scatterFlag=True, DPR=KuNS)
            
            #- Tb --
            vartb   = 'S1/Tb'
            #S1: 10V 10H 19V 19H 23V 37V 37H 89V 89H
            chName  = "10H"
            figPath = figDir + "/TB.%s.wide.%d.png"%(chName,idx)
            cbarPath= figDir + "/cbar.TB.%s.png"%(chName)
            single_fig_TB(DTime, lat, lon, BBoxWide, figPath, cbarPath, chName, scatterFlag=True)
            
            figPath = figDir + "/TB.%s.zoom.%d.png"%(chName,idx)
            cbarPath= figDir + "/cbar.TB.%s.png"%(chName)
            single_fig_TB(DTime, lat, lon, BBoxZoom, figPath, cbarPath, chName, scatterFlag=True)
            
            chName  = "89H"
            figPath = figDir + "/TB.%s.wide.%d.png"%(chName,idx)
            cbarPath = figDir + "/cbar.TB.%s.png"%(chName)
            single_fig_TB(DTime, lat, lon, BBoxWide, figPath, cbarPath, chName, scatterFlag=True)
            
            figPath = figDir + "/TB.%s.zoom.%d.png"%(chName,idx)
            cbarPath = figDir + "/cbar.TB.%s.png"%(chName)
            single_fig_TB(DTime, lat, lon, BBoxZoom, figPath, cbarPath, chName, scatterFlag=True)
            
            #- Join --
            
            figPath0 = figDir + "/IR.wide.%d.png"%(idx)
            figPath1 = figDir + "/KuPR.wide.%d.png"%(idx)
            figPath2 = figDir + "/TB.%s.wide.%d.png"%("10H",idx)
            figPath3 = figDir + "/TB.%s.wide.%d.png"%("89H",idx)
            
            figPath4 = figDir + "/IR.zoom.%d.png"%(idx)
            figPath5 = figDir + "/KuPR.zoom.%d.png"%(idx)
            figPath6 = figDir + "/TB.%s.zoom.%d.png"%("10H",idx)
            figPath7 = figDir + "/TB.%s.zoom.%d.png"%("89H",idx)
            
            cbarPath0 = figDir + "/cbar.IR.png"
            cbarPath1 = figDir + "/cbar.KuPR.png"
            cbarPath2 = figDir + "/cbar.TB.10H.png"
            cbarPath3 = figDir + "/cbar.TB.89H.png"
            
            a2img0 = load_img(figPath0)
            a2img1 = load_img(figPath1)
            a2img2 = load_img(figPath2)
            a2img3 = load_img(figPath3)
            a2img4 = load_img(figPath4)
            a2img5 = load_img(figPath5)
            a2img6 = load_img(figPath6)
            a2img7 = load_img(figPath7)
            
            a2cbar0 = load_img(cbarPath0)
            a2cbar1 = load_img(cbarPath1)
            a2cbar2 = load_img(cbarPath2)
            a2cbar3 = load_img(cbarPath3)
            
            a2line0 = concatenate([a2img0,a2img1, a2img2, a2img3], axis=1)
            a2line1 = concatenate([a2img4,a2img5, a2img6, a2img7], axis=1)
            a2linecbar = concatenate([a2cbar0,a2cbar1, a2cbar2, a2cbar3], axis=1)
            
            
            a2out   = concatenate([a2line0,a2linecbar,a2line1,a2linecbar], axis=0)
            oimg    = Image.fromarray(a2out)
            oPath   = figDir + "/join.snap.%04d.%02d.db.%05d.g.%06d.idx.%05d.png"%(Year,Mon,idx_db, granule, idx_tgt)
            oimg.save(oPath)
            print oPath
            os.remove(figPath0)
            os.remove(figPath1)
            os.remove(figPath2)
            os.remove(figPath3)
            os.remove(figPath4)
            os.remove(figPath5)
            os.remove(figPath6)
            os.remove(figPath7)
            os.remove(cbarPath0)
            os.remove(cbarPath1)
            os.remove(cbarPath2)
            os.remove(cbarPath3)
