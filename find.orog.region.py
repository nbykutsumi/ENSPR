from numpy import *


orogDir  = "/data1/hjkim/GTOPO30"
orogPath = orogDir + "/E060N40.DEM"

"""
BYTEORDER      M
LAYOUT       BIL
NROWS         6000
NCOLS         4800
NBANDS        1
NBITS         16
BANDROWBYTES         9600
TOTALROWBYTES        9600
BANDGAPBYTES         0
NODATA        -9999
ULXMAP        60.00416666666667
ULYMAP        39.99583333333333
XDIM          0.00833333333333
YDIM          0.00833333333333
"""

ny, nx = 6000, 4800
a = fromfile(orogPath, "int16").byteswap().reshape(ny,nx)
print a


