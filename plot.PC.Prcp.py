import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
import os
from numpy import *
from collections import deque
from f_read_db import *
from f_eof11   import *
from sklearn.decomposition import PCA

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

def ret_pc2(a2in):

    pca.fit(a2in.T)
    transformed = pca.fit_transform(a2in.T)
    a1cont  = pca.explained_variance_ratio_
    a2score =transformed.T
    return a2score, a1cont


#-- with PCA library --
pca = PCA()
#----------------------
#thpr = -9999. 
thpr = 0 
#lldbNum =[range(k,k+1000) for k in range(0,9999,1000)]
#lldbNum =[range(k,k+1000) for k in range(0,9000+1,1000)]
lldbNum =[range(k,k+1000) for k in [9000]]
lcorrange = [[0.5,0.51],[0.8,1.0]]
#lcorrange = [[0.05,0.06]]

for corrange in lcorrange:
    cormin, cormax = corrange
    
    for ldbNum in lldbNum:
        csvDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
        srcPath = csvDir + "/pc_prcp.thpr.%d.%05d-%05d.csv"%(thpr,ldbNum[0], ldbNum[-1])
        print "*"*50
        print srcPath
        print "*"*50
        f=open(srcPath,"r")
        lines = f.readlines()
        f.close()
    
        for line in lines[1:]:
            line = map(float, line.split(","))
            num,nrec,cor_pc, cor_prcp, cont, prcp_mean, prcp_std = line
    
            if (np.isnan(cor_pc)or(np.isnan(cor_prcp))or(np.isnan(cont))):
                continue
    
            if not ((cormin<abs(cor_prcp))&(abs(cor_prcp)<cormax)): continue
    

            if nrec < 10: continue
            #if (cont<95)or(96<cont):continue
    
            #-- read DB --------------
            srcDir  = "/data4/utsumi/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"
            idbPath = srcDir + "/" + "db_%05d.bin"%(num)
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
            if thpr>=0:
                a1idx = ma.masked_where(a1prcp_NS <=thpr, arange(nrec)).compressed()
                a1prcp_tmp = a1prcp_NS[a1idx]
                a2pc_tmp   = a2pc[:,a1idx]
            else:
                a1prcp_tmp = a1prcp_NS
                a2pc_tmp   = a2pc
    
            #--- PCA -----------------
            a2score, a1cont = ret_pc(a2pc_tmp)
            #a2score, a1cont = ret_pc2(a2pc_tmp)
            imax_= a1cont.argmax()
            pc1_ = a2score[imax_]
    
            #-- draw plot ------------
            fig   = plt.figure(figsize=(4,4))
            ax    = fig.add_axes([0.2,0.2,0.7,0.7])
    
            x     = pc1_
            y     = a1prcp_tmp        
            ax.plot(x,y,"o",markersize="3",color="k")
            
            # axis label
            ax.set_xlabel("PC'1",fontsize=15)
            ax.set_ylabel("Prcp_NS",fontsize=15)
    
            ax.set_title("%05d (corr=%4.2f, thpr=%d )"%(num, cor_prcp, thpr))
     
            outDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
            util.mk_dir(outDir)
            outPath= outDir + "/plot.PC1.vs.Prcp.cont.%04.2f-%04.2f.thpr.%d.%05d.png"%(cormin, cormax, thpr,num)
            fig.savefig(outPath)
            print outPath
            plt.close()
            
            
