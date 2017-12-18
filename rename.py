import os
import shutil

srcDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
lldbNum =[range(k,k+1000) for k in range(0,9999,1000)]
for ldbNum in lldbNum:
    srcPath = srcDir + "/pc_prcp.%05d-%05d.csv"%(ldbNum[0],ldbNum[-1])
    outPath = srcDir + "/pc_prcp.thpr.-9999.%05d-%05d.csv"%(ldbNum[0],ldbNum[-1])
    shutil.copy(srcPath,outPath)
    print outPath
