import matplotlib
matplotlib.use("Agg")
from numpy import *
from sklearn.decomposition import PCA
import JPLDB as JPLDB
import matplotlib.pyplot as plt
import numpy as np
import sys, os


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
        if ((0.27<corr) & (corr<0.3) & (num<1000)&(num>300)):
        #if ((corr >0.8) & (num>300)):
        #if ((corr >0.8) & (100>num)&(num>30)):
            ldbID.append(dbID)
    return ldbID
#-------------------------------------------------------



#ldbID  = [6363]
ldbID   = ret_ldbID()
varType = "EMIS"
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

jpldb  = JPLDB.JPLDB()
jpldb()
for dbID in ldbID:
    db = jpldb.loadDBsorted(dbID)

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
    a1k, a1km, ikmax, corrmax = findMaxCorr(a2m, a1x)
    if type(a1k) == "NoneType":continue

    #-- Figure -
    fig = plt.figure(figsize=(3,3))
    ax  = fig.add_axes([0.25,0.2, 0.7, 0.7])

    ax.plot(a1km, a1x, ".",color="k",)

    #-- num -
    plt.text(0.3,0.91, "Num=%d"%(len(a1x)),fontsize=12, transform=ax.transAxes)
    #-- title -
    stitle = "%s  ID:%04d (Corr.=%.2f)"%(varType, dbID, corrmax)
    plt.title(stitle)

    #-- axis-label -
    plt.xlabel(varType, fontsize=12)
    plt.ylabel("Precipitation [mm/h]", fontsize=12)

    figDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/pict"
    figPath= figDir + "/plot.%s.vs.Prcp.%04d.png"%(varType,dbID)
    plt.savefig(figPath)
    print figPath


     


