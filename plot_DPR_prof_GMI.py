#! /opt/local/bin/python

from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
#import h5py
import sys
import os
import matplotlib.dates as mdate
import datetime 


ncfile = sys.argv[1]
beampos = int(sys.argv[2])
outfile = "dpr_prof.png"
if ( len(sys.argv) == 4 ):
  outfile = sys.argv[3]

nc = NetCDFFile(ncfile, 'r')

lats = nc.variables['latitude'][:]
lons = nc.variables['longitude'][:]
epoch = nc.variables['time'][:]
tb= nc.variables['Tb'][:]
tb1= 1.0*np.array(nc['MS/Tb1'])
dpr1= 0.01*np.array(nc['MS/z_ku'])
dpr2= 0.01*np.array(nc['MS/z_ka'])
frz2= 0.001*np.array(nc['MS/h273'])
frz= 0.001*np.array(nc['MERRA2/h273'])


clat = getattr(nc,'center_lat')
clon = getattr(nc,'center_lon')
cdate = getattr(nc,'center_date')
#nsms = getattr(nc,'NS_or_MS')
nsms = 'MS'
sat = getattr(nc,'satname')
if ( sat == 'GPM' ): sat = "GMI"

cdate2 = cdate.split()[0]
ctime2 = cdate.split()[1].split(':')
cdate = cdate2 + ' ' + ctime2[0] + ctime2[1] + ' UTC' 

clat = 21
clon = -95

nx = dpr1.shape[0]
ny = dpr1.shape[1]
nz = dpr1.shape[2]
print cdate2, ctime2[0], ctime2[1]
print cdate, clat, clon
print nx,ny,nz

nc.close()


dpr1[dpr1 < 15.0] = np.nan;
dpr2[dpr2 < 15.0] = np.nan;
tb[tb < 0.0] = np.nan;
tb1[tb1 < 0.0] = np.nan;

iscan1 = 0
iscan2 = nx-1
ibin1 = 24  # only up to 16-km
ibin2 = ny-1
#ibin1 = 0
#ibin2 = ny-1


dpr1= dpr1[iscan1:iscan2,beampos,ibin1:ibin2]
dpr2= dpr2[iscan1:iscan2,beampos,ibin1:ibin2]
tb= tb[iscan1:iscan2,beampos,:]
tb1= tb1[iscan1:iscan2,beampos,:]
epoch= epoch[iscan1:iscan2]
frz= frz[iscan1:iscan2,beampos]
frz2= frz2[iscan1:iscan2,beampos]
lats= lats[iscan1:iscan2,beampos]
lons= lons[iscan1:iscan2,beampos]

nx = dpr1.shape[0]
nbin = dpr1.shape[1]
nchan = tb.shape[1]

print 'DPR dims ', dpr1.shape[0], dpr1.shape[1]
print 'TB dims ', tb.shape[0], tb.shape[1]
print 'start coords= ', lats[0], lons[0]
print 'end coords=   ', lats[nx-1], lons[nx-1]

for i in range(0,nx):
  frz[i]= 64*(1-frz[i]/16)
  frz2[i]= 64*(1-frz2[i]/16)


fig = plt.figure(figsize=(19,14),dpi=80)


#aspect = 0.20*nx/nbin
aspect = 0.16*nx/nbin
print nx, nbin, aspect

x = np.arange(0,nx,1)
#y = np.arange(0,15,0.250)
#y = np.arange(0,bin*nbin,bin)
#print "X Y dims=", x.shape,y.shape

plt.set_cmap('jet')
#tick_locs = [87, 79, 71, 63, 55, 47, 39, 31, 23, 15, 7]
#tick_lbls = [' 0', ' 2', ' 4', ' 6', ' 8', '10', '12', '14', '16', '18', '20']
tick_locs = [63, 55, 47, 39, 31, 23, 15, 7, 0]
tick_lbls = [' 0', ' 2', ' 4', ' 6', ' 8', '10', '12', '14', '16']

print 'X dims=', x.shape

ax=plt.subplot(4, 1, 1)
plt.xlim([0,nx])
plt.ylim([nbin-1,0])
title = cdate + '   Ret=' + nsms + '   Top-Ranked DPR Ku-band  ' + sat + '   BeamPos=' + str(beampos)
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
clevs = np.arange(15,45,1)
plt.imshow(dpr1.T, interpolation = 'none', aspect=aspect, vmin=15, vmax=45)
plt.plot(x,frz,'b--',color='blue',linewidth=2)
plt.plot(x,frz2,'b--',color='green',linewidth=2)
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
plt.xlabel('Along-Track Scan Index', labelpad=3, family='Arial', fontsize=14, fontweight='bold')
plt.ylabel('Height Above Ref Surface', labelpad=5, family='Arial', fontsize=14, fontweight='bold')
plt.yticks(tick_locs, tick_lbls, fontsize=14)
plt.xticks(np.arange(0, nx+1, 10), family='Arial', fontsize=16, fontweight='bold')
ax.grid()


ax=plt.subplot(4, 1, 2)
title = cdate + '   Ret=' + nsms + '   Top-Ranked DPR Ka-band  ' + sat + '   BeamPos=' + str(beampos)
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
plt.xlim([0,nx])
plt.ylim([nbin-1,0])
clevs = np.arange(15,45,1)
plt.imshow(dpr2.T, interpolation = 'none', aspect=aspect, vmin=15, vmax=45)
plt.plot(x,frz,'b--',color='blue',linewidth=2)
plt.xlabel('Along-Track Scan Index', labelpad=3, family='Arial', fontsize=14, fontweight='bold')
plt.ylabel('Height Above Ref Surface', labelpad=5, family='Arial', fontsize=14, fontweight='bold')
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
plt.xticks(np.arange(0, nx+1, 10), family='Arial', fontsize=16, fontweight='bold')
ax.grid()


#aspect = 0.20*nx/nchan
aspect = 0.15*nx/nchan

ax=plt.subplot(4, 1, 3)
title = cdate + '   Ret=' + nsms + '   Top-Ranked ' + sat + '   BeamPos=' + str(beampos)
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
plt.xlim([0,nx])
plt.ylim([-0.5,nchan-0.5])
clevs = np.arange(180,290,1)
plt.imshow(tb1.T, interpolation = 'none', aspect=aspect, vmin=180, vmax=290)
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('(Kelvin)', family='Arial', fontsize=16)
plt.xlabel('Along-Track Scan Index', labelpad=3, family='Arial', fontsize=14, fontweight='bold')
plt.ylabel('Channel', family='Arial', labelpad=5, fontsize=14, fontweight='bold')
#plt.yticks(fontsize=14)
plt.yticks(np.arange(0, nchan, 2), family='Arial', fontsize=16, fontweight='normal')
plt.xticks(np.arange(0, nx+1, 10), family='Arial', fontsize=16, fontweight='bold')

plt.text(1, 0, r'10V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 1, r'10H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 2, r'18V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 3, r'18H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 4, r'21V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 5, r'36V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 6, r'36H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 7, r'89V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 8, r'89H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 9, r'166V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 10, r'166H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 11, r'183.31$\pm$3', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 12, r'183.31$\pm$7', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')

ax=plt.subplot(4, 1, 4)
title = 'Observed ' + sat + ' (K)   BeamPos=' + str(beampos)
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
plt.xlim([0,nx])
plt.ylim([-0.5,nchan-0.5])
clevs = np.arange(180,290,1)
plt.imshow(tb.T, interpolation = 'none', aspect=aspect, vmin=180, vmax=290)
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('(Kelvin)', family='Arial', fontsize=16)
plt.xlabel('Along-Track Scan Index', labelpad=3, family='Arial', fontsize=14, fontweight='bold')
plt.ylabel('Channel', family='Arial', labelpad=5, fontsize=14, fontweight='bold')
#plt.yticks(fontsize=14)
plt.yticks(np.arange(0, nchan, 2), family='Arial', fontsize=16, fontweight='normal')
plt.xticks(np.arange(0, nx+1, 10), family='Arial', fontsize=16, fontweight='bold')

plt.text(1, 0, r'10V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 1, r'10H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 2, r'18V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 3, r'18H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 4, r'21V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 5, r'36V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 6, r'36H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 7, r'89V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 8, r'89H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 9, r'166V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 10, r'166H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 11, r'183.31$\pm$3', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(1, 12, r'183.31$\pm$7', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')


#plt.show()
plt.savefig(outfile, transparent='True', bbox_inches='tight', pad_inches=0.25)
sys.exit()


#secs = mdate.epoch2num(epoch)
#print secs
#ax.plot_date(secs, x)
#date_fmt = '%H:%M:%S'
#date_formatter = mdate.DateFormatter(date_fmt)
#ax.xaxis.set_major_formatter(date_formatter)
#fig.autofmt_xdate()

#tick_locs = [0,100,200,300,400,500,600]
#tick_lbls = [ht[0],ht[100],ht[200],ht[300],ht[400],ht[500],ht[600]]
tick_locs = [0,50,100,150,200,250,300,350,400]
tick_lbls = [ht[0],ht[50],ht[100],ht[150],ht[200],ht[250],ht[300],ht[350],ht[400]]

xtick_locs = [0, nx/8, nx/4, 3*nx/8, nx/2, 5*nx/8, 3*nx/4, 7*nx/8, nx-1]
x0= epoch[0]
x2= epoch[nx/4]
x4= epoch[nx/2]
x6= epoch[3*nx/4]
x8= epoch[nx-1]
xtick_lbls = [datetime.datetime.utcfromtimestamp(x0).strftime('%H:%M:%S'),
              datetime.datetime.utcfromtimestamp(x2).strftime('%H:%M:%S'),
              datetime.datetime.utcfromtimestamp(x4).strftime('%H:%M:%S'),
              datetime.datetime.utcfromtimestamp(x6).strftime('%H:%M:%S'),
              datetime.datetime.utcfromtimestamp(x8).strftime('%H:%M:%S')]

#plt.set_cmap('rainbow')

print 'X dims=', x.shape

ax=plt.subplot(3, 1, 1)
plt.xlim([0,nx])
plt.ylim([0,nbin])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'DPR Ku-band  '

secs1 = epoch[0]
secs2 = epoch[nx-1]
print secs1, secs2
date1 = datetime.datetime.utcfromtimestamp(secs1).strftime('%Y/%m/%d %H:%M:%S')
date2 = datetime.datetime.utcfromtimestamp(secs2).strftime('%Y/%m/%d %H:%M:%S')
#print date1,date2
#date1 = datetime.datetime.strptime(sdate1, "%Y/%m/%d %H:%M:%S")
#date2 = datetime.datetime.strptime(sdate2, "%Y/%m/%d %H:%M:%S")
print 'Subsetted date range=',date1, date2


htfrz[htfrz < 0] = np.nan;
isfc[isfc< 0] = np.nan;

ihtfrz= nbin*(htfrz+1500)/12000;

print htfrz
print ihtfrz

#aspect = 0.2*nx/ny
aspect = 0.14*nx/nbin
aspect2 = 5*aspect
print nx, nbin, aspect

x = np.arange(0,nx,1)
#y = np.arange(0,15,0.250)
#y = np.arange(0,bin*nbin,bin)
#print "X Y dims=", x.shape,y.shape

#secs = mdate.epoch2num(epoch)
#print secs
#ax.plot_date(secs, x)
#date_fmt = '%H:%M:%S'
#date_formatter = mdate.DateFormatter(date_fmt)
#ax.xaxis.set_major_formatter(date_formatter)
#fig.autofmt_xdate()

#tick_locs = [0,100,200,300,400,500,600]
#tick_lbls = [ht[0],ht[100],ht[200],ht[300],ht[400],ht[500],ht[600]]
tick_locs = [0,50,100,150,200,250,300,350,400]
tick_lbls = [ht[0],ht[50],ht[100],ht[150],ht[200],ht[250],ht[300],ht[350],ht[400]]


plt.set_cmap('jet')

print 'X dims=', x.shape

ax=plt.subplot(6, 1, 1)
plt.xlim([0,nx])
plt.ylim([0,nbin])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'APR-3 Ku-band  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(-20,40,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(dpr1.T, interpolation = 'none', aspect=aspect, vmin=-20, vmax=40)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
plt.plot(x,isfc,color='black',linewidth=1)
plt.plot(x,ihtfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Time', labelpad=-2, family='Arial', fontsize=14)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'Ku-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')


ax=plt.subplot(6, 1, 2)
plt.xlim([0,nx])

#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'APR-3 Ka-band  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(-20,40,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(dpr2.T, interpolation = 'none', aspect=aspect, vmin=-20, vmax=40)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
plt.plot(x,isfc,color='black',linewidth=1)
plt.plot(x,ihtfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Ray Index', labelpad=-2, family='Arial', fontsize=16)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'Ka-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')




ax=plt.subplot(6, 1, 3)
plt.xlim([0,nx])
plt.ylim([0,nbin])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'APR-3 W-band  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(-20,20,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(dpr3.T, interpolation = 'none', aspect=aspect, vmin=-20, vmax=20)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
plt.plot(x,isfc,color='black',linewidth=1)
plt.plot(x,ihtfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Ray Index', labelpad=-2, family='Arial', fontsize=16)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'W-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')

aspect = 0.14*nx/ntbm
tick_locs = [0,1,2,3,4,5,6,7]
tick_lbls = [1,2,3,4,5,6,7,8]

ax=plt.subplot(6, 1, 4)
#plt.xlim([0,nx])
#plt.ylim([0,ntb])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'MASC  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(200,280,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(tbm.T, interpolation = 'none', aspect=aspect, vmin=200, vmax=280)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
#plt.plot(x,4*dem,color='black',linewidth=2)
#plt.plot(x,4*htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('Kelvin', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Ray Index', labelpad=-2, family='Arial', fontsize=16)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Channel ', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'W-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')
#plt.text(10, 3, 'GMI', family='Arial', fontsize=18, fontweight='bold')
xoff = (iscan2 - iscan1)/14

plt.text(xoff, 0, r'183.31$\pm$1', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 1, r'183.31$\pm$3', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 2, r'183.31$\pm$7', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 3, r'183.31$\pm$8', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 4, r'118.75$\pm$1', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 5, r'118.75$\pm$2', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 6, r'118.75$\pm$7', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(xoff, 7, r'118.75$\pm$8', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')





aspect = 0.14*nx/ntbc
tick_locs = [0,1,2,3,4,5,6,7,8]
tick_lbls = [1,2,3,4,5,6,7,8,9]

ax=plt.subplot(6, 1, 5)
#plt.xlim([0,nx])
#plt.ylim([0,ntb])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'COSMIR  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(200,280,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(tbc.T, interpolation = 'none', aspect=aspect, vmin=200, vmax=280)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
#plt.plot(x,4*dem,color='black',linewidth=2)
#plt.plot(x,4*htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('Kelvin', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Ray Index', labelpad=-2, family='Arial', fontsize=16)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Channel ', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'W-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')
#plt.text(10, 3, 'GMI', family='Arial', fontsize=18, fontweight='bold')
xoff = (iscan2 - iscan1)/60

#50.3, 52.8, 89V, 89H, 165V, 165H, 183+/-1, 183+/-3, 183+/-7

plt.text(xoff, 0, r'50.3', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 1, r'52.8', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 2, r'89V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 3, r'89H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 4, r'165V', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 5, r'165H', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 6, r'183.31$\pm$1', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 7, r'183.31$\pm$3', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 8, r'183.31$\pm$7', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')



aspect = 0.14*nx/ntba
tick_locs = [0,1,2,3,4,5,6,7]
tick_lbls = [1,2,3,4,5,6,7,8]

ax=plt.subplot(6, 1, 6)
#plt.xlim([0,nx])
#plt.ylim([0,ntb])
#title = date1 + ' - ' + date2[11:] + '   (' + clat + ', ' + clon + ')' + '  deltaT= ' + tdiff
title = 'AMPR  ' + date1 + ' - ' + date2
plt.title(title, fontsize=18, family='Arial', fontweight='bold')
#plt.title(title, fontsize=18, fontweight='bold')
clevs = np.arange(200,280,1)
#ax.set_yticks([0,bin*nbin])
plt.imshow(tba.T, interpolation = 'none', aspect=aspect, vmin=200, vmax=280)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,csatT,clevs)
#plt.plot(x,4*dem,color='black',linewidth=2)
#plt.plot(x,4*htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('Kelvin', family='Arial', fontsize=16)
#cbar.set_label('dB (uncorrected)', fontsize=16)
#plt.xlabel('APR3 Ray Index', labelpad=-2, family='Arial', fontsize=16)
#plt.xlabel('CloudSat Ray Index', labelpad=-2, fontsize=16)
plt.ylabel('Channel ', family='Arial', fontsize=16)
#plt.ylabel('Height (km)', fontsize=16)
plt.yticks(tick_locs, tick_lbls, fontsize=14)
#plt.xticks(np.arange(0, nx+1, 500))
plt.xticks(xtick_locs, xtick_lbls, fontsize=14)
ax.grid()
#plt.text(20, 350, 'W-band', family='Arial', fontsize=18, fontweight='bold')
#plt.text(10, 50, 'CloudSat', fontsize=18, fontweight='bold')
#plt.text(10, 3, 'GMI', family='Arial', fontsize=18, fontweight='bold')
xoff = (iscan2 - iscan1)/60


plt.text(xoff, 0, r'10A', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 1, r'10B', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 2, r'19A', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 3, r'19B', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 4, r'37A', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 5, r'37B', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 6, r'85A', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')
plt.text(xoff, 7, r'85B', family='Arial', va='center', ha='left', fontsize=14, fontweight='bold')







#plt.show()
#plt.savefig(ncfile + '.profile.TB2.png', transparent='True', bbox_inches='tight', pad_inches=0.25)
plt.savefig(ncfile + '.profile.TEST.png', transparent='True', bbox_inches='tight', pad_inches=0.25)
sys.exit()


ax=plt.subplot(5, 1, 2)
plt.xlim([0,nx])
plt.ylim([0,nbin])
ax.set_yticks([0,15])
plt.imshow(dpr1T, interpolation = 'none', aspect=aspect, vmin=10, vmax=45)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,dpr1T,clevs)
plt.plot(x,4*elev,color='black',linewidth=2)
plt.plot(x,4*dpr_htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
plt.xlabel('CloudSat Ray Index', labelpad=-2, family='Arial', fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
tick_locs = [0,20,40,60]
tick_lbls = [0,5,10,15]
plt.yticks(tick_locs, tick_lbls)
plt.xticks(np.arange(0, nx+1, 50))
ax.grid()
plt.text(10, 40, 'DPR Ku-band\nNS', family='Arial', fontsize=18, fontweight='bold')
if ( dpr_htfrz[10] > 0 ):
  plt.text(10, 2+4*dpr_htfrz[10], 'T=273 K', family='Arial', fontsize=14, fontstyle='italic', fontweight='bold')

ax=plt.subplot(5, 1, 3)
plt.xlim([0,nx])
plt.ylim([0,60])
ax.set_yticks([0,15])
plt.imshow(dpr2T, interpolation = 'none', aspect=aspect, vmin=10, vmax=45)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,dpr2T,clevs)
plt.plot(x,4*elev,color='black',linewidth=2)
plt.plot(x,4*dpr_htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
plt.xlabel('CloudSat Ray Index', labelpad=-2, family='Arial', fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
tick_locs = [0,20,40,60]
tick_lbls = [0,5,10,15]
plt.yticks(tick_locs, tick_lbls)
plt.xticks(np.arange(0, nx+1, 50))
ax.grid()
plt.text(10, 40, 'DPR Ka-band\nMS', family='Arial', fontsize=18, fontweight='bold')
if ( dpr_htfrz[10] > 0 ):
  plt.text(10, 2+4*dpr_htfrz[10], 'T=273 K', family='Arial', fontsize=14, fontstyle='italic', fontweight='bold')


ax=plt.subplot(5, 1, 4)
plt.xlim([0,nx])
plt.ylim([0,60])
ax.set_yticks([0,15])
plt.imshow(dpr3T, interpolation = 'none', aspect=aspect, vmin=10, vmax=45)
cbar = plt.colorbar(orientation='vertical')
#cs = plt.contourf(x,y,dpr3T,clevs)
plt.plot(x,4*elev,color='black',linewidth=2)
plt.plot(x,4*dpr_htfrz,'b--',color='blue',linewidth=2)
#cbar = plt.colorbar(cs)
cbar.set_label('dB (uncorrected)', family='Arial', fontsize=16)
plt.xlabel('CloudSat Ray Index', labelpad=-2, family='Arial', fontsize=16)
plt.ylabel('Height (km)', family='Arial', fontsize=16)
tick_locs = [0,20,40,60]
tick_lbls = [0,5,10,15]
plt.yticks(tick_locs, tick_lbls)
plt.xticks(np.arange(0, nx+1, 50))
ax.grid()
plt.text(10, 40, 'DPR Ka-band\nHS', family='Arial', fontsize=18, fontweight='bold')
if ( dpr_htfrz[10] > 0 ):
  plt.text(10, 2+4*dpr_htfrz[10], 'T=273 K', family='Arial', fontsize=14, fontstyle='italic', fontweight='bold')


ax=plt.subplot(5, 1, 5)
clevs = np.arange(140,300,2)
#plt.ylim([0,20])
#plt.ylim([0,12])
#cs = plt.contourf(x,y2,tbT,clevs)
plt.imshow(tbT, interpolation = 'none', aspect=aspect2, vmin=140, vmax=300)
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('TB (K)', family='Arial', fontsize=16)
plt.xlabel('CloudSat Ray Index', labelpad=-2, family='Arial', fontsize=16)
#ax.set_yticklabels(['10V','10H','18V','18H','23V','36V','36H','89V','89H','166V','166H','183-3','183-8'])
ax.set_yticklabels([' ']);
plt.xticks(np.arange(0, nx+1, 50))
#ax.set_ylabel('Channel (GMI)', family='Arial', fontsize=16)
#plt.grid()
ax.xaxis.grid(True, which='major')
plt.text(10, 3, 'GMI', family='Arial', fontsize=18, fontweight='bold')
plt.text(-5, 0, 'GMI 89V', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 1, 'GMI 89H', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 2, 'MHS 91', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
#plt.text(-5, 3, '18H', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 4, 'GMI 166V', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 5, 'GMI 166H', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 6, 'MHS 157', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
#plt.text(-5, 7, '89V', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 8, r'GMI 183$\pm$3', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 9, r'MHS 183$\pm$3', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 11, r'GMI 183$\pm$8', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')
plt.text(-5, 12, r'MHS 183$\pm$7', family='Arial', va='center', ha='right', fontsize=14, fontweight='bold')


plt.show()
#plt.savefig(ncfile + '.profile.png', transparent='True', bbox_inches='tight', pad_inches=0.25)
#plt.savefig('profile.png', transparent='True', bbox_inches='tight', pad_inches=0.25)
#print imgname
