import os
import subprocess
import glob

srcDir = "/home/utsumi/bin/JPLCODE/"


outDir = "/home/utsumi/bin/JPLCODE/"
outHost=  "rainbow.iis.u-tokyo.ac.jp"
scmd = "rsync %s -avr utsumi@%s:%s"%(srcDir, outHost, outDir)
print scmd
