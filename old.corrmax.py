import numpy as np
import myfunc.util as util
from   numpy import *
from f_read_db import *
import os, sys
from numpy import linalg
from collections import deque

#ldbNum = range(0,9999,1000)]
#ldbNum = [7992]
ldbNum = range(799,799+100)
#ldbNum = range(8062,8062+1000)

lcorr = deque([])
lnobs = deque([])
lid   = deque([])
limax = deque([])
#print ldbNum
for dbNum in ldbNum:
    print dbNum
    srcDir  = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"
    idbPath = srcDir + "/" + "db_%05d.bin"%(dbNum)
    if not os.path.exists(idbPath):
        print "-"*50
        print "NO File"
        print idbPath
        continue


    unitsize  = 1334  # byte
    #nrec = 132384
    nrec = os.path.getsize(idbPath)/unitsize  # byte
    ldb  = f_read_db.read_db_pc_pr(idbPath, nrec)
    a2pc = ldb[0].T   # shape= (nobs,11) 
    a1pr = ldb[1]

    #----------------------
    # m = a2pc,  x = pr 
    #----------------------
    a2m = a2pc
    a1x = a1pr

    #----------------------
    # remove a1x (a1pr) ==0
    #----------------------
    a1idx = arange(len(a1x)).astype(int32)
    a1idx = ma.masked_where(a1x==0, a1idx).compressed() 

    a2m  = a2m[a1idx,:]
    a1x  = a1x[a1idx]
    #----------------------
    # A = inv(Cm) (Cm,x) inv(Cx) (Cm,x).T  # npc x npc matrix
    #----------------------
    nobs = (a2m.shape)[0]
    npc  = (a2m.shape)[1]
    if nobs==1:
        print "nobs=1, skip"
        continue

    Cmx = array([ sum( (a2m[:,i] - mean(a2m[:,i]))*(a1x-mean(a1x)))/(nobs-1) for i in range(npc)]).reshape(npc,-1)   # (npc x 1)

    Cm    = np.cov(a2m, rowvar=0, bias=0)     # (npc x npc)
    Cx  = np.var(a1x, ddof=1)   # scalar


    #-- check if Cm is invertable --
    if not linalg.cond(Cm) < 1/sys.float_info.epsilon:
        print "Cm is not invertable: skip"
        continue 
    #-------------------------------

    A   = np.dot(linalg.inv(Cm), Cmx)/Cx
    A   = np.dot(A, Cmx.T)
    leval, levect= linalg.eig(A)

    corrmax = 0
    for i in range(npc):
        ktmp  = levect[:,i]
        if np.iscomplex(ktmp).sum()>0: continue

        km = np.dot(ktmp.real, a2m.T)
        corr= np.corrcoef( km, a1x)[0,1]
        if abs(corrmax) < abs(corr):
            corrmax = corr 
            k       = ktmp
            ikmax   = argmax(abs(k))
    print "%.2f"%(corrmax),nobs
    ##-- check if selected corr is max
    #lans = []
    #for ik in range(len(k)):
    #    ktmp = k.copy()
    #    ktmp[ik] = k[ik]*0.99
    #    kmtmp    = np.dot(ktmp.real, a2m.T)
    #    corr0    = np.corrcoef(kmtmp, a1x)[0,1]

    #    ktmp[ik] = k[ik]*1.01
    #    kmtmp    = np.dot(ktmp.real, a2m.T)
    #    corr1    = np.corrcoef(kmtmp, a1x)[0,1]

    #    ans= max(map(abs,[corr0,corrmax,corr1]))==abs(corrmax)
    #    #if ans !=True:
    #    #    print corr0,corrmax,corr1
    #    lans.append(ans)
    #print sum(lans)

    lid.append(dbNum)
    lnobs.append(nobs)
    lcorr.append(corrmax)
    limax.append(ikmax)

#lid   = list(lid)
#lnobs = list(lnobs)
#lcorr = list(lcorr)
#limax = list(limax)

slabel= "id,Nobs,Corr,imax" 
lout  = map(list, zip(lid, lnobs, lcorr, limax))
print lout
sout  = util.list2csv(lout)
oDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
oPath = oDir + "/PC.all.csv"
f= open(oPath,"w"); f.write(sout); f.close()
print oPath



