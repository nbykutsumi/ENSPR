import numpy as np
import myfunc.util as util
from   numpy import *
from f_read_db import *
import os, sys
from numpy import linalg


#lldbNum =[range(k,k+1000) for k in range(0,9999,1000)]
lldbNum = [[7998]]
#lldbNum = [range(7992,7992+10)]
#lldbNum = [range(8062,8062+1000)]
for ldbNum in lldbNum:
    lnum = []
    
    for dbNum in ldbNum:
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
        # A = inv(Cm) (Cm,x) inv(Cx) (Cm,x).T  # npc x npc matrix
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

        print "%.2f"%(corrmax),nobs
        #-- check if selected corr is max
        lans = []
        for ik in range(len(k)):
            ktmp = k.copy()
            ktmp[ik] = k[ik]*0.99
            kmtmp    = np.dot(ktmp.real, a2m.T)
            corr0    = np.corrcoef(kmtmp, a1x)[0,1]
 
            ktmp[ik] = k[ik]*1.01
            kmtmp    = np.dot(ktmp.real, a2m.T)
            corr1    = np.corrcoef(kmtmp, a1x)[0,1]

            ans= max(map(abs,[corr0,corrmax,corr1]))==abs(corrmax)
            #if ans !=True:
            #    print corr0,corrmax,corr1
            lans.append(ans)
        print sum(lans)




