import numpy as np
import myfunc.util as util
from   numpy import *
#from f_read_db import *
import os, sys
from numpy import linalg
from collections import deque
import JPLDB

#lvarType = ["TB","PC"]
lvarType = ["PCEMIS"]

ldbID = range(0,9999+1)
#ldbID = [7992]
#ldbID = range(1042,1042+10)
#ldbID = range(8062,8062+1000)

jpldb = JPLDB.JPLDB()
jpldb()

dNCH  = {"Tb":13, "EMIS":11, "PCEMIS":11}
#-------------------------------------------------------
for varType in lvarType:
    dcorr = {(varType,i):deque([]) for i in range(dNCH[varType])}
dnobs = deque([])
did   = deque([])

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

    elif varType in ["PCPC","PCPC.Nrm","PCEMIS","PCEMIS.Nrm","EMIS"]:
        if varType.split(".")[0]=="PCPC":
            a2in   = db.a2pc_emis    # (nobs, nemis)
        elif varType.split(".")[0] in ["PCEMIS","EMIS"]:
            a2in   = db.a2emis       # (nobs, nemis)

    else:
        print "Check varType",varType
        sys.exit()

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

#lid   = list(lid)
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
    
    
    
