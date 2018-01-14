from f_read_db import *
from f_eof11 import *
import numpy as np
from sklearn.decomposition import PCA
from numpy import *
import os, sys

dbNum   = 834
srcDir  = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"
idbPath = srcDir + "/" + "db_%05d.bin"%(dbNum)

unitsize  = 1334  # byte
nrec = os.path.getsize(idbPath)/unitsize  # byte

ldb  = f_read_db.read_db_multi(idbPath, nrec)
print len(ldb)

