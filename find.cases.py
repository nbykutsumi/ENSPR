from numpy import *
import JPLDB
import glob
import myfunc.util as util

srcDir = "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5"
jpldb = JPLDB.JPLDB()

iYM = [2014,10]
eYM = [2015,8]
#eYM = [2014,9]
lYM = util.ret_lYM(iYM, eYM)

for [Year,Mon] in lYM:
    lsrcPath = sorted(glob.glob(srcDir + "/%04d%02d*.strm"%(Year,Mon)))
  
    a1granule_sub = [] 
    a1idx_sub= [] 
    a1ku_sub = []
    a1nsp_sub= []
    a1gp_sub = []
    a1sfc_sub= []
    a1lat_sub= []
    a1lon_sub= []
    a1year_sub=[]
    a1mon_sub =[]
    a1day_sub =[]
    a1hour_sub=[]
    a1min_sub =[]
    a1sec_sub =[]
    a1tdif_sub=[]
    
    sout = ""
    print srcDir
    for srcPath in lsrcPath:
        #lfileName = ["20150216_015019_GPM_005499.strm"]
        fileName = srcPath.split("/")[-1]
        print fileName
        granule  = fileName.split("_")[3]
        granule  = int(granule.split(".")[0])
        
        
        jpldb.set_file(srcPath)
        a1ku  = jpldb.get_var("precip_NS")
        a1nsp = jpldb.get_var("precip2_NS")  # Near Surface
        a1gp  = jpldb.get_var("precip_GPROF")
        a1sfc = jpldb.get_var("sfc_class")
        a1lat = jpldb.get_var("glat1")
        a1lon = jpldb.get_var("glon1")
        a1year= jpldb.get_var("yyyy")
        a1mon = jpldb.get_var("mm")
        a1day = jpldb.get_var("dd")
        a1hour= jpldb.get_var("hh")
        a1min = jpldb.get_var("mn")
        a1sec = jpldb.get_var("ss")
        a1tdif= jpldb.get_var("timediff")
        
        
        prmin = 30
        #prmin = 3
        a1idx = range(jpldb.nchunks)
        a1idx = ma.masked_where(a1ku < prmin, a1idx)
        a1idx = ma.masked_where(a1gp > a1ku*0.5, a1idx)
    
        a1sfc_mask = ma.masked_outside(a1sfc, 3,11)
        a1idx = ma.masked_where(a1sfc_mask.mask , a1idx)
        
        a1idx = a1idx.compressed()
    
    
    
        if len(a1idx) ==0:
            print "no data to satisfy the conditions"
            continue
        a1granule_sub.extend( ones(len(a1idx)).astype(int32)*granule )
        a1idx_sub  .extend( a1idx )
        a1ku_sub   .extend( a1ku[a1idx] )
        a1nsp_sub  .extend( a1nsp[a1idx])
        a1gp_sub   .extend( a1gp[a1idx] )
        a1sfc_sub  .extend( a1sfc[a1idx])
        a1lat_sub  .extend( a1lat[a1idx])
        a1lon_sub  .extend( a1lon[a1idx])
        a1year_sub .extend( a1year[a1idx])
        a1mon_sub  .extend( a1mon[a1idx])
        a1day_sub  .extend( a1day[a1idx])
        a1hour_sub .extend( a1hour[a1idx])
        a1min_sub  .extend( a1min[a1idx])
        a1sec_sub  .extend( a1sec[a1idx])
        a1tdif_sub .extend( a1tdif[a1idx])


    #-- write --
    lout = [
            list( a1granule_sub )
            ,list( a1idx_sub )
            ,list( a1year_sub )
            ,list( a1mon_sub  )
            ,list( a1day_sub  )
            ,list( a1hour_sub )
            ,list( a1min_sub  )
            ,list( a1sec_sub  )
            ,list( a1lat_sub  )
            ,list( a1lon_sub  )
            ,list( a1sfc_sub  )
            ,list( a1ku_sub   )
            ,list( a1nsp_sub  )
            ,list( a1gp_sub   )
            ]
    
      
    llabel = ["granule","idx","year","mon","day","hour","minute","sec","lat","lon","sfc_class","KuNS_eSurf","KuNS_NearSurf","GPROF"]
    lout = map(list, zip(*lout) )
    lout.insert(0, llabel)
    sout = util.list2csv(lout)
    
    outDir = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv"
    outPath= outDir + "/cases.%04d%02d.csv"%(Year,Mon)
    f=open(outPath, "w"); f.write(sout); f.close()
    print outPath
    
