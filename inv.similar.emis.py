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

a1ave_pc_emis_all, a1std_pc_emis_all = load_ave_std_pc()

#ldbID = [3128]
ldbID = [1159]
#ldbID = ret_ldbID()

for dbID in ldbID:
    print dbID
    db = jpldb.loadDBsorted(dbID)

    a1precip    = db.a1precip_NS
    a1precip_GPROF = db.a1precip_GPROF
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
    #a2emis_sc   = db.a2emis[a1idx_sc]
    #a1precip_sc = a1precip[a1idx_sc]
    #----------------
    # Pickup a target case
    #----------------
    for [i, idx] in enumerate(a1idx_sc):
        #-- check LF-precip --
        kupr  = a1precip[idx]
        gprof = a1precip_GPROF[idx]
        if ( ((17<kupr) & (kupr<19))
            &((10<gprof) & (gprof<15):

            idx_target = idx

    #----------------
    # candidates
    #----------------
    lidx_check_cand = []

    a1pc_emis_obs  = db.a2pc_emis[idx_target]

    # idx_target should be masked out from a1idx_sc here.
    a2pc_emis_use  = db.a2pc_emis[a1idx_sc]
    a1dist  = ((a2pc_emis_use - a1pc_emis_obs)/a1std_pc_emis_all)**2
    a1dist  = sqrt(a1dist.sum(axis=1))

    a2tobesorted = zip(a1dist,a1idx_sc)

    

    #****************************
    #----------------
    # EPC-distance-based
    #----------------
    """
    idx = lidx_check_cand[0]
    a1pc_emis_obs  = db.a2pc_emis[idx]
    a2pc_emis_use  = db.a2pc_emis[a1idx_sc]
    a1dist  = ((a2pc_emis_use - a1pc_emis_obs)/a1std_pc_emis_all)**2
    a1dist  = sqrt(a1dist.sum(axis=1))

    a2tobesorted = zip(a1dist,a1idx_sc)
    a2tobesorted 
    lidx  = ma.masked_where(a1dist
    """



    #--- std ----
    a1std_pc_emis = db.a2pc_emis.std(axis=0)
    a1std_emis    = db.a2emis.std(axis=0)
    a1ave_pc_emis = db.a2pc_emis.mean(axis=0)
    a1ave_emis    = db.a2emis.mean(axis=0)

    #--- check variables ------
    NLEV_NS = 88

    lines = []
    line  = ["idx","precip_NS","precip_GPROF","Year","Mon","Day","Hour","Mnt","SfcClass","Lat","Lon","Elev","Geoid"]
    line  = line + ["PC_emis %d"%(i) for i in range(1,11+1)]
    line  = line + ["Eemis %d"%(i) for i in range(1,11+1)]
    line  = line + ["Ts","T2m","qtot","Ps"]
    line  = line + ["Ta %d"%p for p in db.a1p_prof[0::4]]
    line  = line + ["Qv %d"%p for p in db.a1p_prof[0::4]]
    line  = line + ["Z_ku %.1f"%(0.25*(NLEV_NS-z)) for z in range(1,NLEV_NS+1)[::1]]
    lines.append(line)
    for idx in lidx_check_cand:
        print idx,"timefiff",db.a1timediff[idx]
        line = []
        line = line+ [idx]
        line = line+ [db.a1precip_NS[idx]]
        line = line+ [db.a1precip_GPROF[idx]]
        line = line+ [db.a1Year[idx]]
        line = line+ [db.a1Mon[idx]]
        line = line+ [db.a1Day[idx]]
        line = line+ [db.a1Hour[idx]]
        line = line+ [db.a1Minute[idx]]
        line = line+ [db.a1sfc_class[idx]]
        line = line+ [db.a1glat1[idx]]
        line = line+ [db.a1glon1[idx]]
        line = line+ [db.a1elev[idx]]
        line = line+ [db.a1hs[idx]]
        line = line+ ((db.a2pc_emis[idx]-a1ave_pc_emis) / a1std_pc_emis).tolist()
        line = line+ ((db.a2emis[idx]-a1ave_emis) / a1std_emis).tolist()

        line = line+ [db.a1ts[idx]]
        line = line+ [db.a1t2m[idx]]
        line = line+ [db.a1tqv[idx]]
        line = line+ [db.a1ps[idx]]
        line = line+ db.a2t_prof[idx,0::4].tolist()
        line = line+ db.a2qv_prof[idx,0::4].tolist()
        line = line+ db.a2z_ku[idx,0::1].tolist()
        lines.append(line)

    sout = util.list2csv(lines)
    outDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    outPath = outDir + "/cmpEntries.%05d.csv"%(dbID)
    f=open(outPath,"w"); f.write(sout); f.close()
    print outPath
