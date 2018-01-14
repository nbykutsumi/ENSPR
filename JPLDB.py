import os, sys
from numpy import *
from f_read_db import *

class JPLDB(object):
    def __init__(self):
        self.rootDir = "/home/utsumi/mnt/wellshare/data/JPLDB"

    def __call__(self, dbName="GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"):
        self.srcDir   = self.rootDir + "/" + dbName
        self.unitsize = 1334  # byte

    def loadDBsorted(self, dbID):
        srcPath = self.srcDir + "/db_%05d.bin"%(dbID)

        if not os.path.exists(srcPath):
            print "NO File (in JPLDB.py)"
            sys.exit()
 
        nrec    = os.path.getsize(srcPath)/self.unitsize  # byte
        lout  = f_read_db.read_db_multi(srcPath, nrec)
        self.a1Year  = lout[0]
        self.a1Mon   = lout[1]
        self.a1Day   = lout[2]
        self.a1Hour  = lout[3]
        self.a1Minute= lout[4]
        self.a1Sec   = lout[5]
        self.a1SC_orientation = lout[6]
        self.a1i_S1       = lout[7]
        self.a1j_S1       = lout[8]
        self.a1i_NS       = lout[9]
        self.a1j_NS       = lout[10]

        self.a1sfc_class  = lout[11]
        self.a1precip_NS  = lout[12]
        self.a1precip_max_NS  = lout[13]
        self.a1precip2_NS     = lout[14]
        self.a1precip2_max_NS = lout[15]
        self.a1precip_MS      = lout[16]
        self.a1precip_max_MS  = lout[17]
        self.a1precip2_MS     = lout[18]
        self.a1precip2_max_MS = lout[19]
        self.a1precip_NS_cmb  = lout[20]
        self.a1precip_max_NS_cmb = lout[21]
        self.a1precip_MS_cmb     = lout[22]
        self.a1precip_max_MS_cmb = lout[23]
        self.a1precip_GPROF      = lout[24]
        self.a1frozen_precip_GPROF = lout[25]
        self.a1prob_precip_GPROF   = lout[26]


        self.a1glat1      = lout[27]
        self.a1glon1      = lout[28]
        self.a2tb         = lout[29].T
        self.a2pc_emis    = lout[30].T
        self.a2emis       = lout[31].T
        self.a1elev       = lout[32]
        self.a1ts         = lout[33]
        self.a1t2m        = lout[34]
        self.a1tqv        = lout[35]
        self.a1hs         = lout[36]
        self.a1ps         = lout[37]
        self.a1p_prof     = lout[38]
        self.a2h_prof     = lout[39].T
        self.a2t_prof     = lout[40].T
        self.a2qv_prof    = lout[41].T
        self.a2z_ku       = lout[42].T


        return self
