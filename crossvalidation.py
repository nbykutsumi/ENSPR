from numpy import *
from sklearn.decomposition import PCA
from operator import itemgetter
import JPLDB as JPLDB
import matplotlib.pyplot as plt
import numpy as np
import myfunc.util as util
import sys, os

ldbID = [5894]
nem_compare = 6
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
            a1k     = a1ktmp
            a1km    = a1kmtmp
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

a1ave_pc, a1std_pc = load_ave_std_pc()


for dbID in ldbID:
    db = jpldb.loadDBsorted(dbID)

    a1precip_NS = db.a1precip_NS
    a1idx       = range(len(a1precip_NS))
    ndb_org     = len(a1idx)
    ndb         = ndb_org - 1
    #- pickup higher value --
    prmax = a1precip_NS.max()
    a1idx_hi = ma.masked_where(a1precip_NS < prmax*0.5, a1idx).compressed()
    a1idx_lo = ma.masked_where(a1precip_NS > prmax*0.5, a1idx).compressed()
    
    for idx in a1idx_hi:
    #for idx in a1idx_lo:
        #- Sample-OBS ---
        
        precip_NS_obs = a1precip_NS[idx]
        a1pc_emis_obs = db.a2pc_emis[idx,:]
        
        #- Mask sample-obs data --
        if idx==0:
            a1idx_use  = range(1,ndb_org)
        elif idx==ndb_org-1:
            a1idx_use  = range(ndb_org-1)
        else:
            a1idx_use  = range(idx)+range(idx+1,ndb_org) 


        #----------------
        # EPC based
        #----------------
        a2pc_emis_db  = db.a2pc_emis[a1idx_use]
        a1idx_db      = range(ndb)
        a1dist  = ((a2pc_emis_db - a1pc_emis_obs)/a1std_pc)**2
        a1dist  = sqrt(a1dist.sum(axis=1))

        ##- weight with sorted -   # No need to sort?
        #a2tobesorted = zip(a1dist, a1idx_db)
        #a2tobesorted.sort(key=itemgetter(0), reverse=True)
        #a1dist_sorted, a1idx_db_sorted = zip(*a2tobesorted)
        #a1dist_sorted    = array(a1dist_sorted)
        #a1idx_db_sorted  = array(a1idx_db_sorted)
        #wt = a1dist_sorted / a1dist_sorted[0]
        #wt = exp(-0.5*wt**2)
        #print a1dist_sorted[:7]
        #print wt[1:7]

        #- weight (no sorting)
        dist_min = a1dist.min()
        wt = a1dist /dist_min
        wt = np.exp(-0.5*wt**2)
        swt= wt.sum()
        precip = a1precip_NS[a1idx_use] * wt
        precip = precip.sum()/swt
        print precip, precip_NS_obs, precip - precip_NS_obs
 
        

