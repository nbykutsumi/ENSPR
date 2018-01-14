import numpy as np
import myfunc.util as util
from   numpy import *
#from f_read_db import *
import os, sys
from numpy import linalg
from collections import deque
import JPLDB
from sklearn.decomposition import PCA

lvarType = ["PCPC","PCEMIS"]
#lvarType = ["PCPC"]

#ldbID = range(0,9999+1)
ldbID = [56]
#ldbID = [68]


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
    
        elif varType in ["PCPC","PCPC.Nrm","PCEMIS","PCEMIS.Nrm"]:
            if varType.split(".")[0]=="PCPC":
                a2in   = db.a2pc_emis
            elif varType.split(".")[0]=="PCEMIS":
                a2in   = db.a2emis
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
        pca.fit(a2in)
        a2m    = pca.fit_transform(a2in - a2in.mean(axis=0))   # (nobs, nemis)
      
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
    
    
    
