import subprocess
import sys
srcName = sys.argv[1]
if srcName[-2:]==".c":
    binName = srcName[:-2]
else:
    print "check file name",srcName
#-- compile -------
cmd = "gcc -fpack-struct %s -L/usr/local/lib -lnetcdf -lm -I/usr/local/include -o %s"%(srcName, binName)
subprocess.call(cmd, shell=True)
print "Compiled"
print binName

