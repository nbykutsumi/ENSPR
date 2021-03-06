import numpy as np
import myfunc.util as util
from   numpy import *
#from f_read_db import *
import os, sys
from numpy import linalg
from collections import deque
import JPLDB
from sklearn.decomposition import PCA

#lvarType = ["PCPC","PCEMIS"]
#lvarType = ["EMIS"]
lvarType = ["EMIS_QPRO"]

ldbID = range(0,9999+1)
#ldbID = [7992]
#ldbID = range(1042,1042+10)
#ldbID = range(8062,8062+1000)

jpldb = JPLDB.JPLDB()
jpldb()

pca   = PCA()
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


for varType in lvarType:
    lcorr = deque([])
    lnobs = deque([])
    lid   = deque([])
    limax = deque([])
    #print ldbNum


    if varType in ["EMIS_QPRO"]:
        db      = jpldb.loadDBsorted(1)
        a1dp_tmp = db.a1p_prof[1:] -db.a1p_prof[:-1]
        a1dp     = (a1dp_tmp[:-1]+a1dp_tmp[1:])*0.5
        a1dp     = r_[a1dp_tmp[0], a1dp, a1dp_tmp[-1]]

    for dbID in ldbID:
        #----------------------
        # load JPLDB
        #----------------------
        print dbID
        srcPath = jpldb.srcDir + "/db_%05d.bin"%(dbID)
        if not os.path.exists(srcPath):
            print "NO File"
            continue
    
        db      = jpldb.loadDBsorted(dbID)
        #----------------------
        # m = a2pc,  x = pr 
        #----------------------
        a1x = db.a1precip_NS
    
        if   varType == "PC":
            a2m = db.a2pc_emis
    
        elif varType in ["PCPC","PCPC.Nrm","PCEMIS","PCEMIS.Nrm","EMIS","EMIS_QTOT","EMIS_QPRF","EMIS_QPRO"]:
            if varType.split(".")[0]=="PCPC":
                a2in   = db.a2pc_emis    # (nobs, nemis)
            elif varType.split(".")[0] in ["PCEMIS","EMIS"]:
                a2in   = db.a2emis       # (nobs, nemis)
            elif varType.split(".")[0] in ["EMIS_QTOT"]:
                a2in   = concatenate([db.a2emis, db.a1tqv.reshape(-1,1)],axis=1)

            elif varType.split(".")[0] in ["EMIS_QPRO"]:
                a2qv   = ma.masked_less(db.a2qv_prof,0.)*a1dp
                a1qv_h = a2qv[:,:14].mean(axis=1).reshape(-1,1)
                a1qv_m = a2qv[:,14:28].mean(axis=1).reshape(-1,1) # may be masked at high elevations
                a1qv_l = a2qv[:,28:].mean(axis=1).reshape(-1,1)   # may be masked at high elevations
                a2in   = concatenate([db.a2emis, a1qv_h, a1qv_m, a1qv_l], axis=1)
            else:
                print "Check varType",varType
                sys.exit()
                

        else:
            print "Check varType",varType
            sys.exit()

        #----------------------
        # remove a1x (a1pr) ==0
        #----------------------
        a1idx = arange(len(a1x)).astype(int32)
        a1idx = ma.masked_where(a1x==0, a1idx).compressed() 
    
        a2in = a2in[a1idx,:]
        a1x  = a1x [a1idx]
    
        if a1x.shape[0]==0:
            print "NO Precipitation",dbID
            continue

        #----------------------
        # PCA
        #----------------------
        if varType[:2]=="PC":
            pca.fit(a2in)
            a2m    = pca.fit_transform(a2in - a2in.mean(axis=0))   # (nobs, nemis)
        else:
            a2m    = a2in
      
        #----------------------
        # Normalize
        #----------------------
        if varType.split(".")[-1] == "Nrm": 
            a1cont = pca.explained_variance_ratio_    # (nemis)
            print "a1cont.shape",a1cont.shape
            print "a2m.shape",a2m.shape
            a2m    = a2m / a1cont
    
            print varType, a2m.sum()

        #sys.exit()
        #----------------------
        a1k, a1km, ikmax, corrmax = findMaxCorr(a2m, a1x)
        if a1k == None:continue
    
        nobs = len(a1x)
    
        lid.append(dbID)
        lnobs.append(nobs)
        lcorr.append(corrmax)
        limax.append(ikmax)
    
    #li   = list(lid)
    #lnobs = list(lnobs)
    #lcorr = list(lcorr)
    #limax = list(limax)
    
    slabel= "id,Nobs,Corr,imax" 
    lout  = map(list, zip(lid, lnobs, lcorr, limax))
    sout  = util.list2csv(lout)
    oDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    oPath = oDir + "/%s.all.csv"%(varType)
    f= open(oPath,"w"); f.write(sout); f.close()
    print oPath
    
    
    
