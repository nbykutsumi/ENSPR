import matplotlib
matplotlib.use("Agg")
from f_read_db import *
from f_eof11 import *
import numpy as np
from sklearn.decomposition import PCA
from numpy import *
import os
import myfunc.util as util

def ret_pc(a2in):
    a2cov = np.cov(a2in-a2in.mean(axis=1).reshape(a2in.shape[0],-1))
    l = f_eof11.eof11(a2cov.T)
    a1egval = l[0]
    a1cont  = l[1]
    a2egvct = l[2]
    a2score = np.dot(a2egvct.T, a2in)
    return a2score, a1cont 


    ##--
    #pca = PCA()
    #pca.fit(a2pc.T)
    #transformed = pca.fit_transform(a2pc.T)
    #print transformed.shape
    #print transformed[:,0]
 

#lldbNum =[range(k,k+100) for k in range(1,9999,100)]
lldbNum =[range(k,k+1000) for k in range(0,9999,1000)]
#lldbNum =[range(k,k+1000) for k in range(2000,9999,1000)]
for ldbNum in lldbNum:
    lnum  = []
    lcor_pc=[]
    lcor_prcp=[]
    lcont  =[]
    lprcp_std  = []
    lprcp_mean = []
    lnrec      = []

    for dbNum in ldbNum:
        srcDir  = "/data4/utsumi/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"
        idbPath = srcDir + "/" + "db_%05d.bin"%(dbNum)
        if not os.path.exists(idbPath):
            print "-"*50
            print "NO File"
            print idbPath
            continue


        unitsize  = 1334  # byte
        #nrec = 132384
        nrec = os.path.getsize(idbPath)/unitsize  # byte
        ldb  = f_read_db.read_db(idbPath, nrec)
        a2pc = ldb[0]
        a1prcp_NS = ldb[1]
        
        #--- screen zero precip --
        thpr = 0.
        if thpr>=0:
            a1idx = ma.masked_where(a1prcp_NS <=thpr, arange(nrec)).compressed()
            a1prcp_tmp = a1prcp_NS[a1idx]
            a2pc_tmp   = a2pc[:,a1idx]
        else:
            a1prcp_tmp = a1prcp_NS
            a2pc_tmp   = a2pc  
    
        #--- find 1st PC index ---
        a2score, a1cont = ret_pc(a2pc_tmp) 
        imax_= a1cont.argmax()
        pc1_ = a2score[imax_]
        pc1  = a2pc_tmp[0]
        cor_pc= np.corrcoef(pc1,pc1_)[0][1]
        cor_prcp= np.corrcoef(pc1_, a1prcp_tmp)[0][1]
        cont  = a1cont[imax_]
        prcp_std = np.std (a1prcp_tmp)
        prcp_mean= np.mean(a1prcp_tmp)
    
        lnum.append(dbNum)
        lcor_pc.append(cor_pc)
        lcor_prcp.append(cor_prcp)
        lcont.append(cont)
        lprcp_std.append(prcp_std)
        lprcp_mean.append(prcp_mean)
        lnrec.append(len(a1prcp_tmp))
    #----- write ----------
    a2out= array([lnum,lnrec, lcor_pc, lcor_prcp, lcont, lprcp_mean, lprcp_std]).T
    slabel= "num,nrec,cor(PC1 PC1_),cor(PC1_ Prcp),contri(PC1_),prcp_mean,precp_std\n"
    sout = slabel + util.array2csv(a2out)
    oDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    oPath= oDir + "/pc_prcp.thpr.%d.%05d-%05d.csv"%(thpr,ldbNum[0],ldbNum[-1])
    f=open(oPath,"w");f.write(sout); f.close()
    print oPath 

        




