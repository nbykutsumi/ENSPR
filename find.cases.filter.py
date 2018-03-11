from numpy import *
import myfunc.util as util
import glob
import JPLDB
import sys
from bisect import bisect_right
from operator import itemgetter

Year = 2014
Mon  = 9
srcDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
srcPath= srcDir + "/cases.%04d%02d.csv"%(Year,Mon)
#srcPath= srcDir + "/tmp.cases.%04d%02d.csv"%(Year,Mon)
print srcPath
#nrec_max = 200000
nrec_max = 100000000
ntop   = 20  # top Nth
jplDir = "/home/utsumi/mnt/wellshare/data/JPLDB"
jpl_gra = JPLDB.JPLDB()   # for granule file
jpl_db  = JPLDB.JPLDB()   # for indexed file


def calc_EPC_single(tb, a2coef):
    """
    # tb: (NTBREG,)
    # a2coef: (NREG, NEM)
    # Default dimensions:
    # NTBREG=9, NREG=57, NEM=11
    """
    NTBREG = 9
    NREG   = a2coef.shape[0]-1
    NEM    = a2coef.shape[1]
    """
    print "NEM=",NEM
    print "NTBREG=",NTBREG
    print "NREG=",NREG
    """

    tbcomb = zeros(NREG+1).astype(float32)
    tbcomb[0]= 1.0
    kt = 0 
    for k in range(0,NTBREG):
        tbcomb[kt+1]= tb[k]
        kt = kt+1
        for k1 in range(k,NTBREG):
            tbcomb[kt+1]=tb[k]*tb[k1]
            kt = kt+1


    tbcomb[kt+1]=(tb[0]-tb[1])/(tb[0]+tb[1])
    kt = kt+1
    tbcomb[kt+1]=(tb[2]-tb[3])/(tb[2]+tb[3])
    kt = kt+1
    tbcomb[kt+1]=(tb[5]-tb[6])/(tb[5]+tb[6])

    a1epc = ( a2coef * tbcomb.reshape(NREG+1,1) ).sum(axis=0)
    return a1epc

def load_EPC_coef():
    # PC coefficients
    coefPath = "/home/utsumi/bin/ENSPR/coef_pc.txt"
    f=open(coefPath,"r"); lines = f.readlines(); f.close()
    NREG    = len(lines) -1
    NEM     = len(lines[0].split()) -1
    a2coef  = empty([NREG+1,NEM]).astype(float32)
    for ireg, line in enumerate(lines):
        line = map(float, line.split()[1:])
        a2coef[ireg,:] = line
    return a2coef

def load_ave_std_pc():
    srcPath = "/home/utsumi/bin/ENSPR/ave_pc.txt"
    f=open(srcPath,"r"); lines=f.readlines(); f.close()
    a1ave_pc = empty(len(lines))
    a1std_pc = empty(len(lines))

    for i,line in enumerate(lines):
        line = line.split()
        a1ave_pc[i] = float(line[1])
        a1std_pc[i] = float(line[2])
    return a1ave_pc, a1std_pc


def load_minmax():
    #minmaxPath =  "/home/utsumi/bin/ENSPR/PC_MIN_MAX_10_2pc_overlap.txt"
    minmaxPath =  "/home/utsumi/bin/ENSPR/PC_MIN_MAX_10_no_overlap.txt"
    f=open(minmaxPath, "r"); lines=f.readlines(); f.close()
    
    NEM     = 11
    NPCHIST = 10
    a3minmax = zeros([NEM,NPCHIST,2])
    for iline, line in enumerate(lines):
        line = map(float, line.split()[3:])
        a2line= array(line).reshape(-1,2)
        a3minmax[iline] = a2line
    return a3minmax

def find_indices(a3minmax, a1epc, NEM_USE=4):
    a1idx = [-9999]*NEM_USE
    for i in range(NEM_USE):
        a1idx[i] = bisect_right(a3minmax[i,:,0], a1epc[i]) -1

    return a1idx

def mask_negative(l):
    lout = []
    for x in l:
        if x<0: lout.append("")
        else: lout.append(x)
    return lout 

def mask_value(l, miss):
    lout = []
    for x in l:
        if x==miss: lout.append("")
        else: lout.append(x)
    return lout 


# Filtering --
f = open(srcPath, "r"); lines=f.readlines(); f.close()
ncount_min = 10
#ncount_min = 0

lcase = []
ltmp = []
for il in range(1,len(lines[1:])):
    line0 = lines[il].strip().split(",")
    line1 = lines[il+1].strip().split(",")
    year0,mon0,day0,hour0,mnt0,sec0 = map(int, line0[2:2+6])
    year1,mon1,day1,hour1,mnt1,sec1 = map(int, line1[2:2+6])
    ltmp.append(line0)
    if [year1,mon1,day1,hour1,mnt1] !=[year0,mon0,day0,hour0,mnt0]:
        if len(ltmp)> ncount_min:
            lcase.append(ltmp[ len(ltmp)/2])
        ltmp = []

#lcase=ltmp   # test
sout = util.list2csv(lcase)

#-- save --
outDir = srcDir
outPath= outDir + "/cases.filt.%04d%02d.csv"%(Year,Mon)
f =open(outPath, "w"); f.write(sout); f.close()
print outPath

#----------------------------------------
# Find entries from indexed database
#----------------------------------------
NEM    = 11
NTBREG = 9
NREG   = 57
NEM_USE = 4

a2coef   = load_EPC_coef()
a3minmax = load_minmax()
a1ave_pc_emis_all, a1std_pc_emis_all = load_ave_std_pc()

dbDir0 = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5"
#dbDir0 = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"
dbDir1 = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"

#for case in lcase[3:4]:
for case in lcase:
    granule, idx_gra = map(int, case[:2])
    dbPath0 = glob.glob(dbDir0 + "/*_%06d.strm"%(granule))[0]
    #dbPath0 = dbDir0 + "/db_05554.bin"
    jpl_gra.set_file(srcPath=dbPath0)
    a2tb = jpl_gra.get_var("tb")
    #idx_gra = a2tb.shape[0]/2    # test
    tb  = a2tb[idx_gra]
    

    a1epc  = calc_EPC_single(tb, a2coef)
    a1epc0 = jpl_gra.get_var("pc_emis")[idx_gra]

    IDX    = find_indices(a3minmax, a1epc, NEM_USE=4)

    
    #-- find similar epcs from subset database
    idx_db  = 0
    for i in range(NEM_USE):
        idx_db = idx_db + IDX[i]*(10**i)

    dbPath1 = dbDir1 + "/db_%05d.bin"%(idx_db)
    print ""
    print "dbPath1=",dbPath1
    jpl_db.set_file(srcPath= dbPath1)
    if jpl_db.nchunks < nrec_max:
        nrec = jpl_db.nchunks
    else:
        nrec = nrec_max
    print "read nrec=",nrec
    print "nchunks=",jpl_db.nchunks

    a2pc_emis = jpl_db.get_var("pc_emis", nrec=nrec)
    a1idx_sub = arange( nrec ).astype(int32)

    a2dif   = ((a2pc_emis - a1epc.reshape(-1,NEM))/a1std_pc_emis_all)**2
    a1dist  = sqrt(a2dif.sum(axis=1))
    print "dist min",a1dist.min()
    #-- reduce the number of candidates ----
    dist_max  = a1dist.mean()*0.5
    a1dist_sorted = ma.masked_greater(a1dist, dist_max)
    a1idx_sorted  = ma.masked_where( a1dist_sorted.mask, a1idx_sub)

    a1dist_sorted = a1dist_sorted.compressed()
    a1idx_sorted  = a1idx_sorted.compressed()

    #-- sort ---
    a2sorted= sorted( zip(a1dist_sorted, a1idx_sorted) , key=itemgetter(0))
    a1dist_sorted, a1idx_sorted = zip(*a2sorted)

    #--- check variables ------
    lTBlabels = ["10V", "10H", "19V", "19H", "23V", "37V", "37H", "89V", "89H", "165V", "165H", "183+/-3V", "183+/-8V"]
    NLEV_NS = 88

    lines = []
    line  = ["idx","precip_NS","precip_GPROF","dist","Year","Mon","Day","Hour","Mnt","SfcClass","Lat","Lon","Elev","Geoid"]
    line  = line + ["PC_emis(normed) %d"%(i) for i in range(1,11+1)]
    line  = line + ["Eemis %d"%(i) for i in range(1,11+1)]
    line  = line + ["TB %s"%(s) for s in lTBlabels]
    line  = line + ["Ts","T2m","qtot","Ps"]
    line  = line + ["Ta %d"%p for p in jpl_db.get_var("p_prof")[0][0::4]]
    line  = line + ["Qv %d"%p for p in jpl_db.get_var("p_prof")[0][0::4]]
    line  = line + ["Z_ku %.1f"%(0.25*(NLEV_NS-z)) for z in range(1,NLEV_NS+1)[::1]]
    lines.append(line)

    a1idx_sorted = [-9999] + list(a1idx_sorted)
    for i, idx in enumerate(a1idx_sorted[:ntop]):
        if idx == -9999:
            jpl = jpl_gra
            idx = idx_gra
            a1pc_emis  = a1epc
            dist = 0
        else:
            jpl = jpl_db            
            a1pc_emis  = jpl.get_var("pc_emis", origin=idx, nrec=1)[0]
            dist       = a1dist[idx]

        precip_NS    = jpl.get_var("precip_NS", origin=idx, nrec=1)[0]
        precip_GPROF = jpl.get_var("precip_GPROF", origin=idx, nrec=1)[0]
        Year         = jpl.get_var("yyyy", origin=idx, nrec=1)[0]
        Mon          = jpl.get_var("mm", origin=idx, nrec=1)[0]
        Day          = jpl.get_var("dd", origin=idx, nrec=1)[0]
        Hour         = jpl.get_var("hh", origin=idx, nrec=1)[0]
        Mnt          = jpl.get_var("mn", origin=idx, nrec=1)[0]
        sfc_class    = jpl.get_var("sfc_class", origin=idx, nrec=1)[0]
        glat1        = jpl.get_var("glat1", origin=idx, nrec=1)[0]
        glon1        = jpl.get_var("glon1", origin=idx, nrec=1)[0]
        elev         = jpl.get_var("elev", origin=idx, nrec=1)[0]
        hs           = jpl.get_var("hs", origin=idx, nrec=1)[0]
        a1emis       = jpl.get_var("emis", origin=idx, nrec=1)[0]
        a1tb         = jpl.get_var("tb", origin=idx, nrec=1)[0]
        ts           = jpl.get_var("ts", origin=idx, nrec=1)[0]
        t2m          = jpl.get_var("t2m", origin=idx, nrec=1)[0]
        tqv          = jpl.get_var("tqv", origin=idx, nrec=1)[0]
        ps           = jpl.get_var("ps", origin=idx, nrec=1)[0]
        a1t_prof     = jpl.get_var("t_prof", origin=idx, nrec=1)[0]
        a1qv_prof    = jpl.get_var("qv_prof", origin=idx, nrec=1)[0]
        a1z_ku       = jpl.get_var("z_ku", origin=idx, nrec=1)[0]

        

        line = []
        if ((i==0)&(a1idx_sorted[0]==-9999)):
            line = line+ [-9999]
        else:
            line = line+ [idx]

        line = line+ [precip_NS]
        line = line+ [precip_GPROF]
        line = line+ [dist]
        line = line+ [Year]
        line = line+ [Mon]
        line = line+ [Day]
        line = line+ [Hour]
        line = line+ [Mnt]
        line = line+ [sfc_class]
        line = line+ [glat1]
        line = line+ [glon1]
        line = line+ [elev]
        line = line+ [hs]
        line = line+ ((a1pc_emis-a1ave_pc_emis_all) / a1std_pc_emis_all).tolist()
        line = line+ mask_negative((a1emis).tolist())
        line = line+ (a1tb).tolist()

        line = line+ [ts]
        line = line+ [t2m]
        line = line+ [tqv]
        line = line+ [ps]
        line = line+ mask_negative(a1t_prof[0::4].tolist())
        line = line+ mask_negative(a1qv_prof[0::4].tolist())
        line = line+ mask_value(a1z_ku[0::1].tolist() , -9999)
        lines.append(line)

    sout = util.list2csv(lines)

    outDir  = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    outPath = outDir + "/cmp.case.db.%05d.g.%06d.idx.%05d.csv"%(idx_db, granule, idx_gra)
    #outPath = outDir + "/cmp.Simil.Emis.%05d.csv"%(dbID)

    f=open(outPath, "w"); f.write(sout); f.close()
    print outPath






