from f_read_db import *
from f_eof11 import *
import numpy as np
from sklearn.decomposition import PCA
from numpy import *

idbname = "/home/utsumi/bin/ENSPR/20150830_204112_GPM_008545.strm"

nrec = 132384
ldb  = f_read_db.read_db(idbname, nrec)
a2pc = ldb[0]
print a2pc.shape
a1prcp_NS = ldb[1]
print  a1prcp_NS.shape
print  a1prcp_NS


#--- screen zero precip --
a1idx = ma.masked_where(a1prcp_NS <=20., arange(nrec)).compressed()
a1prcp_NS_nz = a1prcp_NS[a1idx]
a2pc_nz      = a2pc[:,a1idx]

#-------------------------
#a2cov = np.cov(a2pc)
a2cov = np.cov(a2pc_nz)
print a2cov.shape
print "*"*30

l = f_eof11.eof11(a2cov.T)
a1egval = l[0]
a1cont  = l[1]
#a2egvct = l[2].T
a2egvct = l[2]

print a1cont
print sum(a1cont)
print a2egvct[:,0:3]

#a2score = np.dot(a2egvct.T, a2pc)
a2score = np.dot(a2egvct.T, a2pc_nz)
print a2score.shape

#--- find 1st PC index ---
imax = a1cont.argmax()



##--
#pca = PCA()
#pca.fit(a2pc.T)
#transformed = pca.fit_transform(a2pc.T)
#print transformed.shape
#print transformed[:,0]
