from f_read_db import *
from f_eof11 import *
import numpy as np
from sklearn.decomposition import PCA
from numpy import *
import os, sys
import JPLDB2

jpldb = JPLDB2.JPLDB()
jpldb("GMI_dbase_PC_prof2_V5")
gn = jpldb.loadDBgranule(2014,9,1,2885)
