import os, sys
from numpy import *
from f_read_db import *
from collections import OrderedDict as odict
import struct
import numpy as np

'''
Reads JPL observational database format
September 2017

Example input file syntax   20150725_204633_GPM_007985.strm  is provided
In this example, the first record is at 2015/07/25 20:46:33 UTC    7985= GPM orbit revolution

The procedure and of the DPR-GMI matching and calculation of coeffficients to estimate the emissivity principal components (EPC)
from the GMI TB is documented in:

Turk, F.J., Haddad, Z.S. & You, Y., 2016, Estimating Non-Raining Surface Parameters to Assist
GPM Constellation Radiometer Precipitation Algorithms, Journal of Atmospheric and Oceanic Technology, 33(2016), pp. 1333-53

Use "-fpack-struct" gcc compiler option to pack all structure members together without holes
gcc -fpack-struct -o FILE FILE.c  (where FILE.c is the filename of this program)

To test, run as:  FILE 20150725_204633_GPM_007985.strm

sfc_class value:
1= ocean
2= sea ice
3-7 = decreasing vegetation
8-11 = decreasing snow cover
12 = inland water
13 = coast
14 = ocean/sea-ice boundary

NS refers to the retrievals done across the full Ku-band radar swath (49 beam positions), no consideration of Ka-band
MS refers to the retrievals done within NS beam positions 12-36, where both Ku and Ka band radars are available
HS refers to the retrievals done with the Ka-band High Sensivity radar (interleaved with Ka-band)

DESCRIPTION OF EACH RECORD

satid= identifier for the satellite with the radiometer  0=GPM 16-F16 17=F17 18=F18 19=F19 100=ATMS
rev= orbit revolution number for the satellite with the radar (always GPM)
rev2= orbit rev number for satellite with the radiometer (obviously rev2=rev for GMI, but rev2 != rev for SSMIS and ATMS)
SC_orientation= GPM spacecraft orientation, 0 or 180 degrees
i_S1,j_S1= scan and pixel number from radiometer.  For GMI, i_S1 ranges from 0-3000 or so. j_S1 ranges from 0-220.
i_NS,j_NS= scan and pixel number from radar.  For DPR, i_NS ranges from 0-9500 or so. j_NS ranges from 0-48.
sfc_class= TELSEM surface class index, listed above
yyyy,mm,dd,hh,mn,ss= date
timediff = time offset between radar and radiometer, in seconds
slat, slon= spacecraft coordinates of the satellite with the radar, here GPM
glat1, glon1= GMI coordinates
slat2, slon2= spacecraft coordinates of the satellite with the radiometer (slat=slat2 and slon=slon2 for GMI)
tb= GMI TB, from 1B.GPM.GMI, in Kelvin
pc_emis= emissivity principal components (empty, computed below), length 11 array
emis = emissivity that is reconstructed from pc_emis (empty, computed below), length 11 array
sfc_min, sfc_max= land surface type min and max values encountered from the DPR 3x3 profiles
elev = elevation in meters

nku10, nka10, nka10_HS= Number of bins where Z > 10 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku15, nka15, nka15_HS= Number of bins where Z > 15 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku20, nka20, nka20_HS= Number of bins where Z > 20 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku25, nka25, nka25_HS= Number of bins where Z > 25 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
zen_NS= DPR zenith angle, from 2A.GPM.DPR, in degrees
pia_NS, pia_MS= path integrated attenuation in dB for DPR Ku-band only (NS) and Ku+Ka (MS) retrievals
pia_NS_cmb= path integrated attenuation in dB for Ku-band only (NS) combined (DPR+GMI) retrievals
pia_MS_cmb[2]= path integrated attenuation in dB for Ku+Ka band (MS) combined (DPR+GMI) retrievals

precip_NS= Estimated surface precip rate from DPR-only NS retrieval, precipRateESurface in 2A.GPM.DPR
precip_NS_max= max rain rate encountered within the 3x3 profiles

precip2_NS= Near surface precip rate from DPR-only NS retrieval, precipRateNearSurface in 2A.GPM.DPR
precip2_NS_max= max precip rate encountered within the 3x3 profiles

precip_MS= Estimated surface precip rate from DPR-only MS retrieval, precipRateESurface in 2A.GPM.DPR
precip_MS_max= max rain rate encountered within the 3x3 profiles

precip2_MS= Near surface precip rate from DPR-only MS retrieval, precipRateNearSurface in 2A.GPM.DPR
precip2_MS_max= max precip rate encountered within the 3x3 profiles

precip_NS_cmb= Surface precip rate from DPR+GMI combined NS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
precip_NS_cmb_max= max precip rate encountered within the 3x3 profiles

precip_MS_cmb= Surface precip rate from DPR+GMI combined MS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
precip_MS_cmb_max= max precip rate encountered within the 3x3 profiles

precip_GPROF=  GPROF precip, surfacePrecipitation, copied over from 2A.GPM.GMI.GPROF
frozen_precip_GPROF= same as above but frozenPrecipitation, copied over from 2A.GPM.GMI.GPROF
prob_precip_GPROF= probability of precipitation, probabilityOfPrecip copied over from 2A.GPM.GMI.GPROF

ts, t2m= surface temperature (Kelvin), 2-m air temp (Kelvin), interpolated from MERRA2 1-hour reanalysis
tqv= total precip column water (mm), interpolated from MERRA2 3-hourly reanalysis
hs,ps= surface geopotential height (meters) and surface pressure (hPa), interpolated from MERRA2 3-hourly reanalysis
p_prof= pressure levels of MERRA2 (42 levels), interpolated from MERRA2 3-hourly reanalysis, scaled by 10
h_prof= geopotential height profile, in meters
t_prof= temperature profile, in Kelvin
qv_prof= specific humidity profile, in g/g

z_ku= Ku-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
z_ka= Ka-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
'''


class JPLDB(object):
    """
    Based on jpldb.py by hjkim @ IIS.2018-01-18
    """
    def __init__(self): 
        #self.srcDir = "/home/utsumi/mnt/wellshare/data/JPLDB"

        NCHAN   = 13        # 13 channels of GMI Tb
        NEM     = 11        # 9 emissivity [10V 10H 19V 19H 23V 37V 37H 89V 89H], Ts, Column vapor
        NLEV    = 42        # 42 levels of MERRA2
        NLEV_NS = 88        # 88 vertical bins of DPR (250m each)

        dictvars    = odict((
                            ('satid',               'i'),           #
                            ('rev',                 'i'),           #
                            ('rev2',                'i'),           #
                            ('SC_orientation',      'h'),           #
                            ('i_S1',                'h'),           #
                            ('j_S1',                'h'),           #
                            ('i_NS',                'h'),           #
                            ('j_NS',                'h'),           #
                            ('sfc_class',           'h'),           #
                            ('yyyy',                'h'),           #
                            ('mm',                  'h'),           #
                            ('dd',                  'h'),           #
                            ('hh',                  'h'),           #
                            ('mn',                  'h'),           #
                            ('ss',                  'h'),           #
                            ('timediff',            'i'),           #
                            ('slat',                'f'),           #
                            ('slon',                'f'),           #
                            ('glat1',               'f'),           #
                            ('glon1',               'f'),           #
                            ('slat2',               'f'),           #
                            ('slon2',               'f'),           #
                            ('tb',                '%if'%NCHAN),     #
                            ('pc_emis',           '%if'%NEM),       #
                            ('emis',              '%if'%NEM),       #
                            ('emis_cmb',          '%if'%NCHAN),     #
                            ('sfc_min',             'h'),           #
                            ('sfc_max',             'h'),           #
                            ('elev',                'h'),           #
                            ('nku10',               'h'),           #
                            ('nka10',               'h'),           #
                            ('nka10_HS',            'h'),           #
                            ('nku15',               'h'),           #
                            ('nka15',               'h'),           #
                            ('nka15_HS',            'h'),           #
                            ('nku20',               'h'),           #
                            ('nka20',               'h'),           #
                            ('nka20_HS',            'h'),           #
                            ('nku25',               'h'),           #
                            ('nka25',               'h'),           #
                            ('nka25_HS',            'h'),           #
                            ('zen_NS',              'f'),           #
                            ('pia_NS',              'f'),           #
                            ('pia_MS',              'f'),           #
                            ('pia_NS_cmb',          'f'),           #
                            ('pia_MS_cmb',         '2f'),           #
                            ('precip_NS',           'f'),           #
                            ('precip_max_NS',       'f'),           #
                            ('precip2_NS',          'f'),           #
                            ('precip2_max_NS',      'f'),           #
                            ('precip_MS',           'f'),           #
                            ('precip_max_MS',       'f'),           #
                            ('precip2_MS',          'f'),           #
                            ('precip2_max_MS',      'f'),           #
                            ('precip_NS_cmb',       'f'),           #
                            ('precip_max_NS_cmb',   'f'),           #
                            ('precip_MS_cmb',       'f'),           #
                            ('precip_max_MS_cmb',   'f'),           #
                            ('precip_GPROF',        'f'),           #
                            ('prob_precip_GPROF',   'f'),           #
                            ('frozen_precip_GPROF', 'f'),           #
                            ('qual_GPROF',          'i'),           #
                            ('ts',                  'f'),           #
                            ('t2m',                 'f'),           #
                            ('tqv',                 'f'),           #
                            ('hs',                  'f'),           #
                            ('ps',                  'f'),           #
                            ('p_prof',            '%ih'%NLEV),      #
                            ('h_prof',            '%if'%NLEV),      #
                            ('t_prof',            '%if'%NLEV),      #
                            ('qv_prof',           '%if'%NLEV),      #
                            ('z_ku',              '%ih'%NLEV_NS),   #
                            ('z_ka',              '%ih'%NLEV_NS),   #
                           ))

        self.vars,      \
        self.fmts       = zip(*dictvars.items())
        self.dictvars   = dictvars.copy()

        #self.fmtsizes           = [struct.calcsize( fmt ) for fmt in self.fmts]
        self.fmtsize    = struct.calcsize( '<'+''.join(self.fmts) )



    #def __call__(self, dbName="GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"):
    #    self.srcDir   = self.rootDir + "/" + dbName


    def __getattr__(self, k):

        if k in self.dictvars:
            return self.get_var( k )

        else:
            raise KeyError

    def set_file(self, srcPath):
        self.filesize   = os.stat( srcPath ).st_size
        self.nchunks    = self.filesize / self.fmtsize

        self.dbmmap     = np.memmap( srcPath, dtype='S1', mode='r',
                                     shape=(self.nchunks, self.fmtsize) )

        #self.curr       = 0



    def get_var( self, vname, nrec=None, origin=0 ):

        nrec    = self.nchunks if nrec == None       \
             else nrec

        vidx    = self.vars.index( vname )
        fmt     = self.fmts[ vidx ]
        size    = struct.calcsize( fmt )
        seek    = struct.calcsize( '<'+''.join(self.fmts[:vidx]) )

        data    = np.array( self.dbmmap[origin : origin+nrec,
                                        seek   : seek+size] )

        data.dtype  = fmt[-1]

        if data.shape[-1] == 1:
            data.shape  = data.shape[:-1]

        return data


    

class JPLDB_old(object):
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
        self.a1timediff  = lout[6]
        self.a1SC_orientation = lout[7]
        self.a1i_S1       = lout[8]
        self.a1j_S1       = lout[9]
        self.a1i_NS       = lout[10]
        self.a1j_NS       = lout[11]

        self.a1sfc_class  = lout[12]
        self.a1precip_NS  = lout[13]
        self.a1precip_max_NS  = lout[14]
        self.a1precip2_NS     = lout[15]
        self.a1precip2_max_NS = lout[16]
        self.a1precip_MS      = lout[17]
        self.a1precip_max_MS  = lout[18]
        self.a1precip2_MS     = lout[19]
        self.a1precip2_max_MS = lout[20]
        self.a1precip_NS_cmb  = lout[21]
        self.a1precip_max_NS_cmb = lout[22]
        self.a1precip_MS_cmb     = lout[23]
        self.a1precip_max_MS_cmb = lout[24]
        self.a1precip_GPROF      = lout[25]
        self.a1frozen_precip_GPROF = lout[26]
        self.a1prob_precip_GPROF   = lout[27]


        self.a1glat1      = lout[28]
        self.a1glon1      = lout[29]
        self.a2tb         = lout[30].T
        self.a2pc_emis    = lout[31].T
        self.a2emis       = lout[32].T
        self.a1elev       = lout[33]
        self.a1ts         = lout[34]
        self.a1t2m        = lout[35]
        self.a1tqv        = lout[36]
        self.a1hs         = lout[37]
        self.a1ps         = lout[38]
        self.a1p_prof     = lout[39]
        self.a2h_prof     = lout[40].T
        self.a2t_prof     = lout[41].T
        self.a2qv_prof    = lout[42].T
        self.a2z_ku       = lout[43].T


        return self
