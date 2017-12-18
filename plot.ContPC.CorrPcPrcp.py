import matplotlib
matplotlib.use("Agg")
import numpy as np
from collections import deque
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util

srcDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
lthpr  = [-9999,0]
#lthpr  = [-9999]

for thpr in lthpr:
    lnrec   = deque([])
    lcor_pc = deque([])
    lcor_prcp=deque([])
    lcont   = deque([])
    lprcp_mean= deque([])
    lprcp_std = deque([])
    
    
    #lldbNum =[range(k,k+1000) for k in range(0,9999,1000)]
    lldbNum =[range(k,k+1000) for k in [9000]]
    for ldbNum in lldbNum:
        srcPath = srcDir + "/pc_prcp.thpr.%d.%05d-%05d.csv"%(thpr, ldbNum[0], ldbNum[-1])
        f=open(srcPath,"r")
        lines = f.readlines()
        f.close()
    
        for line in lines[1:]:
            line = map(float, line.split(","))
            num,nrec,cor_pc, cor_prcp, cont, prcp_mean, prcp_std = line
    
    
            if (np.isnan(cor_pc)or(np.isnan(cor_prcp))or(np.isnan(cont))):
                continue
            lcor_pc.append(cor_pc)
            lcor_prcp.append(cor_prcp)
            lcont.append(cont)
            lprcp_mean.append(prcp_mean)
            lprcp_mean.append(prcp_std)
    
    #-- draw plot ------------
    x     = map(abs, lcont)
    y     = map(abs, lcor_prcp)
    fig   = plt.figure(figsize=(4,4))
    ax    = fig.add_axes([0.2,0.2,0.7,0.7])
    
    
    ax.plot(x,y,"o",markersize="1",color="k")
    
    # axis label
    ax.set_xlabel("Cont(PC'1)")
    ax.set_ylabel("Corr(PC'1, Prcp)")

    # title
    stitle = "prcp > %d [mm/hour]"%(thpr)
    plt.title(stitle)

    
    outDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
    util.mk_dir(outDir)
    outPath= outDir + "/plot.2D.Cont.vs.CorPrcp.thpr.%d.png"%(thpr)
    fig.savefig(outPath)
    print outPath
    plt.close()
    
    #-- draw 2D-histogram ------- 
    x     = map(abs, lcont)
    y     = map(abs, lcor_prcp)
    xbins  = arange(50,100+0.01,1)
    ybins  = arange(0,1+0.001,0.02)
    
    l = np.histogram2d(x, y, bins=[xbins,ybins])
    H = l[0].T
    
    xbnd= l[1]
    ybnd= l[2]
    X,Y = np.meshgrid(xbnd,ybnd)
    
    
    fig   = plt.figure(figsize=(4,4))
    ax    = fig.add_axes([0.2,0.2,0.7,0.7])
    
    # axis label
    ax.set_xlabel("Cont(PC'1)", fontsize=15)
    ax.set_ylabel("Corr(PC'1, Prcp)",fontsize=15)
    
    H = ma.masked_equal(H,0)
    im  = ax.pcolormesh(X,Y,H)

    # title
    stitle = "prcp > %d [mm/hour]"%(thpr)
    plt.title(stitle)

    outDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
    util.mk_dir(outDir)
    outPath= outDir + "/hist.2D.Cont.vs.CorPrcp.thpr.%d.png"%(thpr)
    fig.savefig(outPath)
    print outPath
    plt.close()
    
    # colorbar
    figcbar = plt.figure(figsize=(0.5,1.6))
    axcbar  = figcbar.add_axes([0.1,0.1,0.2,0.8])
    plt.colorbar(im, cax=axcbar, orientation="vertical")
    cbarPath = outDir + "/cbar.hist.2D.Cont.vs.CorPrcp.thpr.%d.png"%(thpr)
    figcbar.savefig(cbarPath)
    plt.close()
