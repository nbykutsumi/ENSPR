import matplotlib
matplotlib.use("Agg")
from numpy import *
from sklearn.decomposition import PCA
from sklearn import linear_model
from operator import itemgetter
from bisect import bisect_right
import JPLDB as JPLDB
import matplotlib.pyplot as plt
import numpy as np
import myfunc.util as util
import sys, os
from collections import deque



#-------------------------------------------------------
def ret_ldbID():
    csvDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    csvPath= csvDir + "/EMIS.all.csv"
    f = open(csvPath,"r")
    lines = f.readlines()
    f.close()

    ldbID = []
    for line in lines:
        line = line.split(",")
        dbID = int(line[0])
        num  = int(line[1])
        corr = float(line[2])
        #if ((0.27<corr) & (corr<0.3) & (num<1000)&(num>300)):
        #if ((corr >0.7) & (num>1000) & (num<10000) & (dbID>325)):
        if ((corr >0.7) & (800>num)&(num>700)):
            ldbID.append(dbID)
    return ldbID

#-------------------------------------------------------
def findMaxCorr(a2m, a1x):
    # shape of a2m = (nobs,11)
    # shape of a1x = (nobs)

    #----------------------
    # A = inv(Cm) (Cm,x) inv(Cx) (Cm,x).T  # npc x npc matrix
    #----------------------
    nobs = (a2m.shape)[0]
    npc  = (a2m.shape)[1]
    if nobs <=1:
        print "nobs=1, exit"
        return None,None,None,None

    Cmx = array([ sum( (a2m[:,i] - mean(a2m[:,i]))*(a1x-mean(a1x)))/(nobs-1) for i in range(npc)]).reshape(npc,-1)   # (npc x 1)

    Cm    = np.cov(a2m, rowvar=0, bias=0)     # (npc x npc)
    Cx  = np.var(a1x, ddof=1)   # scalar


    #-- check if Cm is invertable --
    if not linalg.cond(Cm) < 1/sys.float_info.epsilon:
        print "Cm is not invertable: skip"
        return None,None,None,None
    #-------------------------------

    A   = np.dot(linalg.inv(Cm), Cmx)/Cx
    A   = np.dot(A, Cmx.T)
    leval, levect= linalg.eig(A)

    corrmax = 0
    for i in range(npc):
        a1ktmp  = levect[:,i]
        if np.iscomplex(a1ktmp).sum()>0: continue

        a1kmtmp = np.dot(a1ktmp.real, a2m.T)
        corr= np.corrcoef( a1kmtmp, a1x)[0,1]
        if abs(corrmax) < abs(corr):
            corrmax = corr
            a1k     = np.real(a1ktmp)
            a1km    = np.real(a1kmtmp)
            ikmax   = argmax(abs(a1k))
    #print "%.2f"%(corrmax),nobs
    return a1k, a1km, ikmax, corrmax
#-------------------------------------------------------
def ret_lmean(a1x,a1y):
    bins = linspace(0,a1x.max()*1.2,20)
    BINS = zip(bins[:-1],bins[1:])
    lmean = [ma.masked_where(logical_or(a1x<binmin, binmax<=a1x), a1y).mean() for (binmin,binmax) in BINS]
    lxmean= [mean(BIN) for BIN in BINS]
    return lxmean, lmean
#-------------------------------------------------------
def CDF_matching(x, a1x, a1y, sortedFlag=False, db_maxrec=None):
    if db_maxrec !=None:
        a1x = a1x[:db_maxrec]
        a1y = a1y[:db_maxrec]

    if sortedFlag == False:
        a1x.sort()
        a1y.sort()
    else:
        pass

    lm = linear_model.LinearRegression()
    ndat = a1x.shape[0]
    ndat_regres = min(5, ndat)
    if x <= a1x[0]:
        a1x = a1x.reshape(-1,1)
        a1y = a1y.reshape(-1,1)
        lm.fit(a1x[:ndat_regres], a1y[:ndat_regres])
        a = lm.coef_
        b = lm.intercept_
        y = a*x+b
        y = y[0][0]
    elif a1x[-1] <= x:
        a1x = a1x.reshape(-1,1)
        a1y = a1y.reshape(-1,1)
        lm.fit(a1x[-ndat_regres:], a1y[-ndat_regres:])
        a = lm.coef_
        b = lm.intercept_
        y = a*x+b
        y = y[0][0]
    else:
        idx = bisect_right(a1x, x)-1
        y0  = a1y[idx]
        y1  = a1y[idx+1]
        x0  = a1x[idx]
        x1  = a1x[idx+1]

        if x0==x1:
            y = a1y[idx]
        else:
            y = (y1-y0)*(x-x0)/(x1-x0) +y0

    return y


#-------------------------------------------------------
def load_ave_std_pc():
    srcPath = "./ave_pc.txt"
    f=open(srcPath,"r"); lines=f.readlines(); f.close()
    a1ave_pc = empty(len(lines))
    a1std_pc = empty(len(lines)) 

    for i,line in enumerate(lines):
        line = line.split()
        a1ave_pc[i] = float(line[1])
        a1std_pc[i] = float(line[2])
    return a1ave_pc, a1std_pc

#-------------------------------------------------------

jpldb  = JPLDB.JPLDB()
jpldb()

a1ave_pc, a1std_pc = load_ave_std_pc()

#ldbID = [1159]
#ldbID = ret_ldbID()

for dbID in ldbID:
    print dbID
    db = jpldb.loadDBsorted(dbID)

    a1precip    = db.a1precip_NS
    a1precip_gprof= db.a1precip_GPROF
    ndb_org     = len(a1precip)
    a1idx       = range(ndb_org)
    ndb         = ndb_org - 1
    #- pickup higher value --
    prmax = a1precip.max()
    a1idx_hi = ma.masked_where(a1precip < prmax*0.5, a1idx).compressed()
    a1idx_lo = ma.masked_where(a1precip > prmax*0.5, a1idx).compressed()

    #----------------
    # screen zero-precip
    #----------------
    a1idx_sc    = range(len(a1precip))
    a1idx_sc    = ma.masked_where(a1precip <=0, a1idx_sc).compressed()
    a2emis_sc   = db.a2emis[a1idx_sc]
    a1precip_sc = a1precip[a1idx_sc]
    a1precip_gprof_sc = a1precip_gprof[a1idx_sc]


    #----------------
    # Linear fitting
    #----------------
    nsamples = len(a1precip_sc)
    a1std    = a2emis_sc.std(axis=0)
    a2x      = a2emis_sc/a1std
    a1y      = a1precip_sc
    lm = linear_model.LinearRegression()
    lm.fit(a2x, a1y)
    
    a1coef   = lm.coef_
    intersept= lm.intercept_
    score    = lm.score(a2x, a1y)
    a1est_LF = np.dot(a2x, a1coef.reshape(-1,1)) + intersept
    a1est_LF = a1est_LF.reshape(-1,)
    a1est_LF = ma.masked_less(a1est_LF,0).filled(0.0) 

    #****************************
    a1est_epc = deque([])
    a1est_cdf = deque([])
 
    #for idx in a1idx_hi:
    #for idx in a1idx_lo:
    for idx in a1idx_sc:
        #- extract sample OBS ----

        precip_obs = a1precip[idx]
        a1pc_emis_obs = db.a2pc_emis[idx,:]
        a1emis_obs = db.a2emis[idx,:]
        
        #- Mask sample-obs data --
        a1idx_use = ma.masked_equal(a1idx_sc, idx).compressed()


        #"""
        #----------------
        # EPC-distance-based
        #----------------
        a2pc_emis_use  = db.a2pc_emis[a1idx_use]
        a1dist  = ((a2pc_emis_use - a1pc_emis_obs)/a1std_pc)**2
        a1dist  = sqrt(a1dist.sum(axis=1))

        #- weight (no sorting)
        dist_min = a1dist.min()
        wt = a1dist /dist_min
        wt = np.exp(-0.5*wt**2)
        swt= wt.sum()
        precip = a1precip[a1idx_use] * wt
        precip = precip.sum()/swt
        a1est_epc.append(precip)
        #print ""
        #print "%.2f"%precip, "%.2f"%precip_obs, "%.2f"%(precip - precip_obs)
        #"""
        #""" 
        #----------------
        # CDF-matching (Corr.Max)
        #---------------- 
        a2emis_use    = db.a2emis[a1idx_use]
        a1precip_use  = a1precip[a1idx_use]
        a1k, a1kx, ikmax, corrmax = findMaxCorr(a2emis_use, a1precip_use)
        if corrmax <0:
            a1k  = -a1k
            a1kx = -a1kx

        x      = (a1k*a1emis_obs).sum()
        precip = CDF_matching(x, a1kx, a1precip_use, sortedFlag=False)
        a1est_cdf.append(precip)
        #print "%.2f"%precip, "%.2f"%precip_obs, "%.2f"%(precip - precip_obs)
        #""" 
    #"""
    #************************
    # Figures
    #------------------------
    #-- figure -
    fig = plt.figure(figsize=(3,3))
    ax  = fig.add_axes([0.25,0.2, 0.7, 0.7])

    #-- scatter -
    a1obs  = a1precip_sc
    ax.plot(a1obs, a1est_LF, "o",color="limegreen",markersize=1)
    ax.plot(a1obs, a1est_epc, "o",color="orange",markersize=1)
    ax.plot(a1obs, a1est_cdf, "o",color="b",markersize=1)
    ax.plot(a1obs, a1precip_gprof_sc, "o",color="magenta",markersize=1)

    #-- mean lines-
    a1xmean, a1mean_LF  = ret_lmean(a1obs, a1est_LF )
    a1xmean, a1mean_epc = ret_lmean(a1obs, a1est_epc)
    a1xmean, a1mean_cdf = ret_lmean(a1obs, a1est_cdf)
    a1xmean, a1mean_gprof = ret_lmean(a1obs, a1precip_gprof_sc)
    ax.plot(a1xmean, a1mean_LF , "-", color="lime", linewidth=2)
    ax.plot(a1xmean, a1mean_epc, "-", color="orange", linewidth=2)
    ax.plot(a1xmean, a1mean_cdf, "-", color="royalblue", linewidth=2)
    ax.plot(a1xmean, a1mean_gprof, "-", color="magenta", linewidth=2)


    #-- 1:1 line -
    x1 = a1obs.max()*1.2
    ax.plot([0,x1],[0,x1],"--",color="k", linewidth=0.5)

    #-- axis-lim---
    ax.set_xlim([-0.2, x1])
    ax.set_ylim([-0.2, x1])

    #-- num -
    plt.text(0.3,0.91, "Num=%d"%(len(a1y)),fontsize=12, transform=ax.transAxes)
    #-- title -
    stitle = "ID:%04d"%(dbID)
    plt.title(stitle)

    #-- axis-label -
    plt.xlabel("Precip (Obs) [mm/h]", fontsize=12)
    plt.ylabel("Precip (Est) [mm/h]", fontsize=12)

    figDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
    figPath= figDir + "/scatter.comp.%04d.png"%(dbID)
    plt.savefig(figPath)
    print figPath
    plt.clf()
    #"""


 
