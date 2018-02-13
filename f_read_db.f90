MODULE f_read_db

CONTAINS

!sfc_class value:
!1= ocean
!2= sea ice
!3-7 = decreasing vegetation
!8-11 = decreasing snow cover
!12 = inland water
!13 = coast
!14 = ocean/sea-ice boundary
!
!NS refers to the retrievals done across the full Ku-band radar swath (49 beam positions), no consideration of Ka-band
!MS refers to the retrievals done within NS beam positions 12-36, where both Ku and Ka band radars are available
!HS refers to the retrievals done with the Ka-band High Sensivity radar (interleaved with Ka-band)
!
!
!DESCRIPTION OF EACH RECORD
!
!satid= identifier for the satellite with the radiometer  0=GPM 16-F16 17=F17 18=F18 19=F19 100=ATMS
!rev= orbit revolution number for the satellite with the radar (always GPM)
!rev2= orbit rev number for satellite with the radiometer (obviously rev2=rev for GMI, but rev2 != rev for SSMIS and ATMS)
!SC_orientation= GPM spacecraft orientation, 0 or 180 degrees
!i_S1,j_S1= scan and pixel number from radiometer.  For GMI, i_S1 ranges from 0-3000 or so. j_S1 ranges from 0-220.
!i_NS,j_NS= scan and pixel number from radar.  For DPR, i_NS ranges from 0-9500 or so. j_NS ranges from 0-48.
!sfc_class= TELSEM surface class index, listed above
!yyyy,mm,dd,hh,mn,ss= date
!timediff = time offset between radar and radiometer, in seconds
!slat, slon= spacecraft coordinates of the satellite with the radar, here GPM
!glat1, glon1= GMI coordinates
!slat2, slon2= spacecraft coordinates of the satellite with the radiometer (slat=slat2 and slon=slon2 for GMI)
!tb= GMI TB, from 1B.GPM.GMI, in Kelvin
!pc_emis= emissivity principal components (empty, computed below), length 11 array
!emis = emissivity that is reconstructed from pc_emis (empty, computed below), length 11 array
!sfc_min, sfc_max= land surface type min and max values encountered from the DPR 3x3 profiles
!elev = elevation in meters
!
!nku10, nka10, nka10_HS= Number of bins where Z > 10 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku15, nka15, nka15_HS= Number of bins where Z > 15 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku20, nka20, nka20_HS= Number of bins where Z > 20 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku25, nka25, nka25_HS= Number of bins where Z > 25 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!zen_NS= DPR zenith angle, from 2A.GPM.DPR, in degrees
!pia_NS, pia_MS= path integrated attenuation in dB for DPR Ku-band only (NS) and Ku+Ka (MS) retrievals
!pia_NS_cmb= path integrated attenuation in dB for Ku-band only (NS) combined (DPR+GMI) retrievals
!pia_MS_cmb[2]= path integrated attenuation in dB for Ku+Ka band (MS) combined (DPR+GMI) retrievals
!
!precip_NS= Estimated surface precip rate from DPR-only NS retrieval, precipRateESurface in 2A.GPM.DPR
!precip_NS_max= max rain rate encountered within the 3x3 profiles
!
!precip2_NS= Near surface precip rate from DPR-only NS retrieval, precipRateNearSurface in 2A.GPM.DPR
!precip2_NS_max= max precip rate encountered within the 3x3 profiles
!
!precip_MS= Estimated surface precip rate from DPR-only MS retrieval, precipRateESurface in 2A.GPM.DPR
!precip_MS_max= max rain rate encountered within the 3x3 profiles
!
!precip2_MS= Near surface precip rate from DPR-only MS retrieval, precipRateNearSurface in 2A.GPM.DPR
!precip2_MS_max= max precip rate encountered within the 3x3 profiles
!
!precip_NS_cmb= Surface precip rate from DPR+GMI combined NS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
!precip_NS_cmb_max= max precip rate encountered within the 3x3 profiles
!
!precip_MS_cmb= Surface precip rate from DPR+GMI combined MS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
!precip_MS_cmb_max= max precip rate encountered within the 3x3 profiles
!
!precip_GPROF=  GPROF precip, surfacePrecipitation, copied over from 2A.GPM.GMI.GPROF
!frozen_precip_GPROF= same as above but frozenPrecipitation, copied over from 2A.GPM.GMI.GPROF
!prob_precip_GPROF= probability of precipitation, probabilityOfPrecip copied over from 2A.GPM.GMI.GPROF
!
!ts, t2m= surface temperature (Kelvin), 2-m air temp (Kelvin), interpolated from MERRA2 1-hour reanalysis
!tqv= total precip column water (mm), interpolated from MERRA2 3-hourly reanalysis
!hs,ps= surface geopotential height (meters) and surface pressure (hPa), interpolated from MERRA2 3-hourly reanalysis
!p_prof= pressure levels of MERRA2 (42 levels), interpolated from MERRA2 3-hourly reanalysis, scaled by 10
!h_prof= geopotential height profile, in meters
!t_prof= temperature profile, in Kelvin
!qv_prof= specific humidity profile, in g/g
!
!z_ku= Ku-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
!z_ka= Ka-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100



!*****************************************************
!Format of PC_MIN_MAX file (11 lines)
!EPC     min        max       bin0_min   bin0_max    bin1_min    bin1_max      then same for bins 2-9
!
!  0 -199.35802  320.60181    -60.01204  -22.19836  -22.19836  -15.78338  -15.78338  -11.00258  -11.00258   -5.36850   -5.36850    0.60104    0.60104    6.83276    6.83276   13.54176   13.54176   20.40424   20.40424   29.00614   29.00614   45.48376
!  1 -182.16125  145.55750    -30.19480   -7.04496   -7.04496   -4.60922   -4.60922   -2.87756   -2.87756   -1.49092   -1.49092   -0.20506   -0.20506    1.25172    1.25172    3.15748    3.15748    6.05210    6.05210   10.28600   10.28600   24.38654
!  2   -1.76337    1.15656     -1.40294   -0.96568   -0.96568   -0.66566   -0.66566   -0.24926   -0.24926    0.19294    0.19294    0.24864    0.24864    0.28968    0.28968    0.34658    0.34658    0.40812    0.40812    0.47638    0.47638    0.65992
!  3   -0.95563    1.32334     -0.40594   -0.12074   -0.12074   -0.06542   -0.06542   -0.04040   -0.04040   -0.02248   -0.02248   -0.00452   -0.00452    0.01480    0.01480    0.03384    0.03384    0.05490    0.05490    0.08490    0.08490    0.51144
!  4   -0.55545    0.41323     -0.17628   -0.02770   -0.02770   -0.01578   -0.01578   -0.00886   -0.00886   -0.00346   -0.00346    0.00146    0.00146    0.00620    0.00620    0.01114    0.01114    0.01650    0.01650    0.02282    0.02282    0.17140
!  5   -0.32075    0.38067     -0.08374   -0.01780   -0.01780   -0.01238   -0.01238   -0.00770   -0.00770   -0.00346   -0.00346    0.00026    0.00026    0.00384    0.00384    0.00794    0.00794    0.01378    0.01378    0.02476    0.02476    0.14406
!  6   -0.48772    0.87471     -0.19236   -0.01520   -0.01520   -0.00834   -0.00834   -0.00482   -0.00482   -0.00204   -0.00204    0.00056    0.00056    0.00332    0.00332    0.00668    0.00668    0.01142    0.01142    0.02162    0.02162    0.46860
!  7   -0.47112    0.42215     -0.06606   -0.00920   -0.00920   -0.00606   -0.00606   -0.00412   -0.00412   -0.00250   -0.00250   -0.00088   -0.00088    0.00094    0.00094    0.00332    0.00332    0.00692    0.00692    0.01330    0.01330    0.09990
!  8   -0.35452    0.40279     -0.09296   -0.00660   -0.00660   -0.00374   -0.00374   -0.00234   -0.00234   -0.00128   -0.00128   -0.00030   -0.00030    0.00070    0.00070    0.00186    0.00186    0.00330    0.00330    0.00562    0.00562    0.07870
!  9   -0.99769    0.31942     -0.09440   -0.00514   -0.00514   -0.00288   -0.00288   -0.00172   -0.00172   -0.00090   -0.00090   -0.00018   -0.00018    0.00056    0.00056    0.00138    0.00138    0.00248    0.00248    0.00476    0.00476    0.06016
! 10   -0.18611    0.44873     -0.03368   -0.00364   -0.00364   -0.00224   -0.00224   -0.00134   -0.00134   -0.00060   -0.00060    0.00002    0.00002    0.00066    0.00066    0.00138    0.00138    0.00220    0.00220    0.00344    0.00344    0.04300
!*****************************************************





!--------------------------------------------------
SUBROUTINE read_db_multi(idbname, nrec&
                        ,a1yyyy,a1mm,a1dd,a1hh,a1mn,a1ss &
                        ,a1timediff &
                        ,a1SC_orientation, a1i_S1, a1j_S1, a1i_NS, a1j_NS &
                        ,a1sfc_class &
                        ,a1precip_NS &
                        ,a1precip_max_NS &
                        ,a1precip2_NS &
                        ,a1precip2_max_NS &
                        ,a1precip_MS &
                        ,a1precip_max_MS &
                        ,a1precip2_MS &
                        ,a1precip2_max_MS &
                        ,a1precip_NS_cmb & 
                        ,a1precip_max_NS_cmb & 
                        ,a1precip_MS_cmb & 
                        ,a1precip_max_MS_cmb & 
                        ,a1precip_GPROF &
                        ,a1frozen_precip_GPROF &
                        ,a1prob_precip_GPROF &
                        ,a1glat1, a1glon1 &
                        ,a2tb &
                        ,a2pc_emis &
                        ,a2emis  &
                        ,a1elev &
                        ,a1ts, a1t2m, a1tqv  &
                        ,a1hs, a1ps  &
                        ,a1p_prof  &
                        ,a2h_prof, a2t_prof, a2qv_prof &
                        ,a2z_ku  &
                        )

implicit none

!-- parameters ---------
integer,parameter :: NLEV_NS = 88  ! DPR 250-m bins
integer,parameter :: NLEV = 42  ! from MERRA2  42 levels
integer,parameter :: NCHAN =13  ! GMI TB
integer,parameter :: NEM = 11  ! 11= first nine emis then Ts then total column vapor
integer,parameter :: NTBREG = 9
integer,parameter :: NREG = 57
integer,parameter :: NEM_USE = 4
integer,parameter :: NPCHIST = 10
integer,parameter :: BYTE = 1334  ! UTSUMI


!-- input --------------
character(len=512)  idbname
!f2py intent(in)    idbname

integer             nrec
!f2py intent(in)    nrec

!-- out -----------------
integer*2,dimension(nrec) :: a1yyyy,a1mm,a1dd,a1hh,a1mn,a1ss, a1timediff
!f2py intent(out)            a1yyyy,a1mm,a1dd,a1hh,a1mn,a1ss, a1timediff

integer*2,dimension(nrec) :: a1SC_orientation, a1i_S1, a1j_S1, a1i_NS, a1j_NS
!f2py intent(out)            a1SC_orientation, a1i_S1, a1j_S1, a1i_NS, a1j_NS

integer*2,dimension(nrec) :: a1sfc_class
!f2py intent(out)            a1sfc_class

real,dimension(nrec)    :: a1precip_NS, a1precip_max_NS
!f2py intent(out)          a1precip_NS, a1precip_max_NS
real,dimension(nrec)    :: a1precip2_NS, a1precip2_max_NS
!f2py intent(out)          a1precip2_NS, a1precip2_max_NS

real,dimension(nrec)    :: a1precip_MS, a1precip_max_MS
!f2py intent(out)          a1precip_MS, a1precip_max_MS
real,dimension(nrec)    :: a1precip2_MS, a1precip2_max_MS
!f2py intent(out)          a1precip2_MS, a1precip2_max_MS


real,dimension(nrec)    :: a1precip_NS_cmb, a1precip_max_NS_cmb
!f2py intent(out)          a1precip_NS_cmb, a1precip_max_NS_cmb
real,dimension(nrec)    :: a1precip_MS_cmb, a1precip_max_MS_cmb
!f2py intent(out)          a1precip_MS_cmb, a1precip_max_MS_cmb


real,dimension(nrec)    :: a1precip_GPROF, a1frozen_precip_GPROF
!f2py intent(out)          a1precip_GPROF, a1frozen_precip_GPROF
real,dimension(nrec)    :: a1prob_precip_GPROF
!f2py intent(out)          a1prob_precip_GPROF

real,dimension(nrec)    :: a1glat1, a1glon1
!f2py intent(out)          a1glat1, a1glon1

real,dimension(nchan,nrec) :: a2tb
!f2py intent(out)             a2tb

real,dimension(nem,nrec) :: a2emis
!f2py intent(out)           a2emis

real,dimension(nem,nrec) :: a2pc_emis
!f2py intent(out)           a2pc_emis

real,dimension(nrec)    :: a1ts, a1t2m, a1tqv
!f2py intent(out)          a1ts, a1t2m, a1tqv

real,dimension(nrec)    :: a1hs, a1ps
!f2py intent(out)          a1hs, a1ps

integer*2,dimension(nlev) :: a1p_prof
!f2py intent(out)            a1p_prof

real,dimension(nlev,nrec)      :: a2h_prof, a2t_prof, a2qv_prof
!f2py intent(out)                 a2h_prof, a2t_prof, a2qv_prof

real,dimension(nlev_NS,nrec) :: a2z_ku
!f2py intent(out)               a2z_ku

integer*2,dimension(nrec) :: a1elev
!f2py intent(out)            a1elev
!------------------------



! record structure (1334 bytes)

type strm_record
    integer satid, rev, rev2
    integer*2 SC_orientation,i_S1,j_S1,i_NS,j_NS,sfc_class
    integer*2 yyyy,mm,dd,hh,mn,ss
    integer timediff
    real slat,slon,glat1,glon1, slat2, slon2
    real tb(NCHAN)
    real pc_emis(NEM), emis(NEM)
    real emis_cmb(NCHAN)
    integer*2 sfc_min,sfc_max,elev
    integer*2 nku10, nka10, nka10_HS
    integer*2 nku15, nka15, nka15_HS
    integer*2 nku20, nka20, nka20_HS
    integer*2 nku25, nka25, nka25_HS
    real zen_NS
    real pia_NS, pia_MS
    real pia_NS_cmb, pia_MS_cmb(2)
    real precip_NS, precip_max_NS
    real precip2_NS, precip2_max_NS
    real precip_MS, precip_max_MS
    real precip2_MS, precip2_max_MS
    real precip_NS_cmb, precip_max_NS_cmb
    real precip_MS_cmb, precip_max_MS_cmb
    real precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF
    integer qual_GPROF
    real ts, t2m, tqv
    real hs, ps
    integer*2 p_prof(NLEV)
    real h_prof(NLEV), t_prof(NLEV), qv_prof(NLEV)
    integer*2 z_ku(88)
    integer*2 z_ka(88)
end type strm_record

type(strm_record) prof

! Added by UTSUMI

integer             ios

!--------------

integer i, j, k, m, nbad, bytes, irec
integer :: nrec_all=0, ndbf=0, nf=0, nrain=0
integer i1, iclutter
real thick, qv, delp, kPa, abstemp, densair, rmsd, ht, z_ku, z_ka, hs1
character(len=32) pdate

real ave_emis(NEM), std_emis(NEM)
real ave_pc(NEM), std_pc(NEM)
!real b_all(NREG+1)(NEM)  ! regression coeffs for all NTBREG
real b_all(NEM,NREG+1)   ! regression coeffs for all NTBREG
real tmp(NREG+1), psum
!real pc_e(NEM), e(NEM), u(NEM)(NEM)
real pc_e(NEM), e(NEM), u(NEM,NEM)
character(len=512) svar
integer id, k1, k2, kt
real coeff_disc(NEM), disc, disc_mid
real x1,x2,y1,y2,z1,z2

real pc_e_min(NEM), pc_e_max(NEM)
!real pcmin, pcmax, pc_range(NEM)(NPCHIST)(2), pc_lo, pc_hi
real pcmin, pcmax, pc_range(2,NPCHIST,NEM), pc_lo, pc_hi
integer idx, idxn(NEM_USE), found, idx2

integer ndb_files
integer p, pow2(NEM_USE), ndb_files_used

integer,allocatable,dimension(:) :: nrec_dbfile
real gmi_freq(13) 
data gmi_freq/10.7,10.7,18.6,18.6,21.65,36.5,36.5,89.0,89.0,166.0,166.0,186.31,190.31/


ndb_files= NPCHIST**NEM_USE
allocate( nrec_dbfile(ndb_files) )
nrec_dbfile = 0

!bytes= sizeof(strm_record)
bytes= 1334

!print *,"record length=",bytes
!printf("N of EPC used=%d  N of EPC histogram bins=%d   N of DB files= %d\n", NEM_USE, NPCHIST, ndb_files);
!printf("N of EPC used=%d  N of EPC histogram bins=%d   N of DB files= %d\n", NEM_USE, NPCHIST, ndb_files);

!nrec_dbfile= (int *) malloc(ndb_files*sizeof(int));
!for (k= 0; k< ndb_files; k++) nrec_dbfile[k]= 0;

do k= 0,NEM_USE-1
    pow2(k)= NPCHIST**k
end do

!!  Read EPC 10-bin histogram file */
!!  This is used to compute the database indices from the EPC structure */
!
!svar="PC_MIN_MAX_10_no_overlap.txt"
!!svar="PC_MIN_MAX_10_2pc_overlap.txt"
!open(10,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
!pc_range = -9999.
!do i = 1,NEM
!    read(10,*) k, pcmin, pcmax, (pc_range(:,j,i), j=1,NPCHIST)
!    pc_range(1,1,i)       = pcmin-100.0
!    pc_range(2,NPCHIST,i) = pcmax+100.0
!    !print *,j,pc_range(1,:,i),pc_range(2,:,i)
!
!end do
!close(10)


!  Read EPC coefficients file 
!  This calculates EPC terms from the TB combinations 

svar="coef_pc.txt"
open(11,file=svar,status="old")
!print '("Opened ",A)', trim(svar)
b_all = -9999.
do i = 1,NREG+1
    read(11,*) k,(b_all(j,i), j=1,NEM)
    !print *,i,b_all(:,i)
end do
close(11)


!  Read eigenvalues matrix file
!  This reconstructs the emissivity vector from the EPC

svar="wmatrix.txt"
open(12,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
u = -9999.
do i = 1,NEM
    read(12,*) (u(j,i), j=1,NEM)
    !print *,i,u(:,i)
end do
close(12)


!  Read ave emissivity file

svar="ave_emis.txt"
open(13,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
ave_emis = -9999.
do i = 1,NEM
    read(13,*) k,ave_emis(i), std_emis(i)
    !print *,i,ave_emis(i)
end do
close(13)


!!  Read ave EPC file
!svar="ave_pc.txt"
!open(14,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
!ave_pc = -9999.
!do i = 1,NEM
!    read(14,*) k,ave_pc(i), std_pc(i)
!    !print *,i,ave_pc(i), std_pc(i)
!end do
!close(14)



!---------------------------------------------------------------------------
!print *,trim(idbname)
open(15, file=trim(idbname), form="unformatted", status="old", access="direct", recl=byte)
nrain= 0

ios = 0    ! For gfortran
irec=0
do irec =1,nrec
    !irec=irec+1
    read(15,rec=irec, iostat=ios) prof
    !if (ios.ne.0) exit

    !---------------------------------------------------
    ! quality control on each DB record 

    hs1= -99.99
    do i=2,NLEV
        if (prof%h_prof(i).le.0.0) then
            hs1 = prof%h_prof(i-1)
            exit   ! necessary? 
        end if     
    end do

    if ( hs1.lt.0.0) hs1= prof%h_prof(NLEV)
    if ( prof%hs .gt. hs1 )then
        print *,"hs larger than lowest geopotential ht"
        print *,"hs=",prof%hs
        print *,"h_lowest=", hs1
        stop        
    end if
 
    nbad= 0

    do i=1,9
        if ((prof%tb(i).lt.20.0).or.(prof%tb(i).gt.350)) nbad=nbad+1
    end do
    if ( nbad .gt.0) cycle
    if ( (prof%precip_NS .lt.0.0).and.(prof%precip_NS_cmb .lt.0.0) ) cycle
    if ( prof%sfc_class  .lt. 0 ) cycle
    if (( abs(prof%glat1) .gt. 90).and.(abs(prof%glon1) .gt. 180 )) cycle


    if (( (prof%yyyy .lt. 2014) .or. (prof%yyyy .gt. 2020 )) .or. &
      ( (prof%mm .lt. 1) .or. (prof%mm .gt. 12 )).or. &
      ( (prof%dd .lt. 1) .or. (prof%dd .gt. 31 )).or. &
      ( (prof%hh .lt. 0) .or. (prof%hh .gt. 23 )).or. &
      ( (prof%mn .lt. 0) .or. (prof%mn .gt. 59 )).or. &
      ( (prof%ss .lt. 0) .or. (prof%ss .gt. 59 )))then

      print *,"Bad date or time", prof%yyyy, prof%mm, prof%dd, prof%hh, prof%mn, prof%ss
      cycle

    end if

    !---------------------------------------------------*/

    ! precip
    a1precip_NS(irec)     = prof%precip_NS
    a1precip_max_NS(irec) = prof%precip_max_NS
    a1precip2_NS(irec)    = prof%precip2_NS
    a1precip2_max_NS(irec)= prof%precip2_max_NS

    a1precip_MS(irec)     = prof%precip_MS
    a1precip_max_MS(irec) = prof%precip_max_MS
    a1precip2_MS(irec)    = prof%precip2_MS
    a1precip2_max_MS(irec)= prof%precip2_max_MS


    a1precip_NS_cmb(irec)     = prof%precip_NS_cmb
    a1precip_max_NS_cmb(irec) = prof%precip_max_NS_cmb
    a1precip_MS_cmb(irec)     = prof%precip_MS_cmb
    a1precip_max_MS_cmb(irec) = prof%precip_max_MS_cmb

    a1precip_GPROF(irec)        = prof%precip_GPROF
    a1frozen_precip_GPROF(irec) = prof%frozen_precip_GPROF
    a1prob_precip_GPROF(irec)   = prof%prob_precip_GPROF


    ! set up regression variables
    tmp(1)= 1.0
    kt= 1
    do k= 1,NTBREG
        tmp(kt+1)= prof%tb(k)
        kt=kt+1
        do k1= k,NTBREG
            tmp(kt+1)= prof%tb(k)*prof%tb(k1)
            kt=kt+1
        end do
    end do
    tmp(kt+1)= (prof%tb(1)-prof%tb(2))/(prof%tb(1)+prof%tb(2))
    kt = kt+1
    tmp(kt+1)= (prof%tb(3)-prof%tb(4))/(prof%tb(3)+prof%tb(4))
    kt = kt+1
    tmp(kt+1)= (prof%tb(6)-prof%tb(7))/(prof%tb(6)+prof%tb(7))
    kt = kt+1



    ! emissivity PC
    !write(*,fmt="(A)",advance="no") "EPC= "
    do k= 1,NEM
        psum= 0.0
        do k2= 1,NREG+1
            psum = psum+ b_all(k,k2)*tmp(k2)
            !print *, k,k2,b_all(k,k2),tmp(k2)
        end do
        pc_e(k)= psum
        prof%pc_emis(k)= pc_e(k)
        !write(*, fmt='(f8.4 )',advance="no") pc_e(k)
    end do
    !print *,""

    ! reconstruct emissivites from PC
    do k= 1,NEM
        psum= 0.0

        do k2= 1,NEM
            psum = psum + u(k2,k) * pc_e(k2)
        end do
        e(k)= psum + ave_emis(k)
        prof%emis(k)= e(k)
        !if ( k .eq. 1 ) write(*,fmt="(A)",advance="no") "emis= "
        !if ( k .eq. NEM-1 ) write(*,fmt="(A)",advance="no") " Ts= "
        !if ( k .eq. NEM ) write(*,fmt="(A)",advance="no") "  Vap= "
        !write(*,"(f8.3 )",advance="no"), e(k)
    end do
    !print *,""


    ! Profile of geopot height, pressure, temperature, specific humidity

    ! DPR radar reflectivity profiles
    do j= 1,NLEV_NS
        ht= 0.25*(NLEV_NS-j)
        !write(*,fmt='(I3," ",f6.3)',advance="no") j, ht
        z_ku= 0.01*prof%z_ku(j)
        iclutter= 0
        if (( z_ku .lt. 0.0) .and. (z_ku .gt. -99.0 )) then
            iclutter= 1
            z_ku = z_ku*(-1.0)
        end if
            continue
        a2z_ku(j,irec) = z_ku

        !write(*,fmt='("  Ku=",I1," ",f6.2)',advance="no") iclutter, z_ku
        !if (( prof%j_NS .lt. 12).or.(prof%j_NS .gt. 36 )) then
        !    print *,""
        !    cycle
        !end if
        !z_ka= 0.01*prof%z_ka(j)
        !iclutter= 0
        !if ( (z_ka .lt. 0.0) .and. (z_ka .gt. -99.0 ))then
        !    iclutter= 1
        !    z_ka = z_ka*(-1.0)
        !end if
        !write(*,fmt='("  Ka= ",I1," ",f6.2)') iclutter, z_ka


    end do
    !if ( prof%precip_NS_cmb .gt. 0.0 ) nrain=nrain+1

    ! Profile of pressure level: only one time
    if (a1p_prof(1).eq.0)then
        a1p_prof(:)  = prof%p_prof
    end if
    !--------------------------
    a1yyyy(irec)      = prof%yyyy
    a1mm(irec)        = prof%mm
    a1dd(irec)        = prof%dd
    a1hh(irec)        = prof%hh
    a1mn(irec)        = prof%mn
    a1ss(irec)        = prof%ss
    a1timediff(irec)  = prof%timediff
    a1SC_orientation(irec) = prof%SC_orientation
    a1i_S1(irec)      = prof%i_S1
    a1j_S1(irec)      = prof%j_S1
    a1i_NS(irec)      = prof%i_NS
    a1j_NS(irec)      = prof%j_NS
    a1sfc_class(irec) = prof%sfc_class
    a1glat1(irec)     = prof%glat1
    a1glon1(irec)     = prof%glon1
    a2tb(:,irec)      = prof%tb
    a2pc_emis(:,irec) = prof%pc_emis
    a2emis(:,irec)    = prof%emis
    a1elev(irec)      = prof%elev
    a1ts(irec)        = prof%ts
    a1t2m(irec)       = prof%t2m
    a1tqv(irec)       = prof%tqv
    a1hs(irec)        = prof%hs
    a1ps(irec)        = prof%ps
    a2h_prof(:,irec)  = prof%h_prof
    a2t_prof(:,irec)  = prof%t_prof
    a2qv_prof(:,irec) = prof%qv_prof
    !a2z_ku(:,irec)    = prof%z_ku   ! multiplied by 100

end do
close(15)

!ndbf = sum(nrec_dbfile)
!!write(*,'(A,"  N=",I7,"  Nraining=",I7, "  N_dbfiles=",I7)') trim(idbname),irec, nrain, ndbf 
return 
END SUBROUTINE read_db_multi




!--------------------------------------------------
SUBROUTINE read_db_pc_pr(idbname, nrec, a2pc, a1precip_NS)
implicit none

!-- parameters ---------
integer,parameter :: NLEV_NS = 88  ! DPR 250-m bins
integer,parameter :: NLEV = 42  ! from MERRA2  42 levels
integer,parameter :: NCHAN =13  ! GMI TB
integer,parameter :: NEM = 11  ! 11= first nine emis then Ts then total column vapor
integer,parameter :: NTBREG = 9
integer,parameter :: NREG = 57
integer,parameter :: NEM_USE = 4
integer,parameter :: NPCHIST = 10
integer,parameter :: BYTE = 1334  ! UTSUMI


!-- input --------------
character(len=512)  idbname
!f2py intent(in)    idbname

integer             nrec
!f2py intent(in)    nrec

!-- out -----------------
real,dimension(11,nrec) :: a2pc
!f2py intent(out)          a2pc

real,dimension(nrec)    :: a1precip_NS
!f2py intent(out)          a1precip_NS
!------------------------
!sfc_class value:
!1= ocean
!2= sea ice
!3-7 = decreasing vegetation
!8-11 = decreasing snow cover
!12 = inland water
!13 = coast
!14 = ocean/sea-ice boundary
!
!NS refers to the retrievals done across the full Ku-band radar swath (49 beam positions), no consideration of Ka-band
!MS refers to the retrievals done within NS beam positions 12-36, where both Ku and Ka band radars are available
!HS refers to the retrievals done with the Ka-band High Sensivity radar (interleaved with Ka-band)
!
!
!DESCRIPTION OF EACH RECORD
!
!satid= identifier for the satellite with the radiometer  0=GPM 16-F16 17=F17 18=F18 19=F19 100=ATMS
!rev= orbit revolution number for the satellite with the radar (always GPM)
!rev2= orbit rev number for satellite with the radiometer (obviously rev2=rev for GMI, but rev2 != rev for SSMIS and ATMS)
!SC_orientation= GPM spacecraft orientation, 0 or 180 degrees
!i_S1,j_S1= scan and pixel number from radiometer.  For GMI, i_S1 ranges from 0-3000 or so. j_S1 ranges from 0-220.
!i_NS,j_NS= scan and pixel number from radar.  For DPR, i_NS ranges from 0-9500 or so. j_NS ranges from 0-48.
!sfc_class= TELSEM surface class index, listed above
!yyyy,mm,dd,hh,mn,ss= date
!timediff = time offset between radar and radiometer, in seconds
!slat, slon= spacecraft coordinates of the satellite with the radar, here GPM
!glat1, glon1= GMI coordinates
!slat2, slon2= spacecraft coordinates of the satellite with the radiometer (slat=slat2 and slon=slon2 for GMI)
!tb= GMI TB, from 1B.GPM.GMI, in Kelvin
!pc_emis= emissivity principal components (empty, computed below), length 11 array
!emis = emissivity that is reconstructed from pc_emis (empty, computed below), length 11 array
!sfc_min, sfc_max= land surface type min and max values encountered from the DPR 3x3 profiles
!elev = elevation in meters
!
!nku10, nka10, nka10_HS= Number of bins where Z > 10 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku15, nka15, nka15_HS= Number of bins where Z > 15 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku20, nka20, nka20_HS= Number of bins where Z > 20 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!nku25, nka25, nka25_HS= Number of bins where Z > 25 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
!zen_NS= DPR zenith angle, from 2A.GPM.DPR, in degrees
!pia_NS, pia_MS= path integrated attenuation in dB for DPR Ku-band only (NS) and Ku+Ka (MS) retrievals
!pia_NS_cmb= path integrated attenuation in dB for Ku-band only (NS) combined (DPR+GMI) retrievals
!pia_MS_cmb[2]= path integrated attenuation in dB for Ku+Ka band (MS) combined (DPR+GMI) retrievals
!
!precip_NS= Estimated surface precip rate from DPR-only NS retrieval, precipRateESurface in 2A.GPM.DPR
!precip_NS_max= max rain rate encountered within the 3x3 profiles
!
!precip2_NS= Near surface precip rate from DPR-only NS retrieval, precipRateNearSurface in 2A.GPM.DPR
!precip2_NS_max= max precip rate encountered within the 3x3 profiles
!
!precip_MS= Estimated surface precip rate from DPR-only MS retrieval, precipRateESurface in 2A.GPM.DPR
!precip_MS_max= max rain rate encountered within the 3x3 profiles
!
!precip2_MS= Near surface precip rate from DPR-only MS retrieval, precipRateNearSurface in 2A.GPM.DPR
!precip2_MS_max= max precip rate encountered within the 3x3 profiles
!
!precip_NS_cmb= Surface precip rate from DPR+GMI combined NS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
!precip_NS_cmb_max= max precip rate encountered within the 3x3 profiles
!
!precip_MS_cmb= Surface precip rate from DPR+GMI combined MS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
!precip_MS_cmb_max= max precip rate encountered within the 3x3 profiles
!
!precip_GPROF=  GPROF precip, surfacePrecipitation, copied over from 2A.GPM.GMI.GPROF
!frozen_precip_GPROF= same as above but frozenPrecipitation, copied over from 2A.GPM.GMI.GPROF
!prob_precip_GPROF= probability of precipitation, probabilityOfPrecip copied over from 2A.GPM.GMI.GPROF
!
!ts, t2m= surface temperature (Kelvin), 2-m air temp (Kelvin), interpolated from MERRA2 1-hour reanalysis
!tqv= total precip column water (mm), interpolated from MERRA2 3-hourly reanalysis
!hs,ps= surface geopotential height (meters) and surface pressure (hPa), interpolated from MERRA2 3-hourly reanalysis
!p_prof= pressure levels of MERRA2 (42 levels), interpolated from MERRA2 3-hourly reanalysis, scaled by 10
!h_prof= geopotential height profile, in meters
!t_prof= temperature profile, in Kelvin
!qv_prof= specific humidity profile, in g/g
!
!z_ku= Ku-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
!z_ka= Ka-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100






! record structure (1334 bytes)

type strm_record
    integer satid, rev, rev2
    integer*2 SC_orientation,i_S1,j_S1,i_NS,j_NS,sfc_class
    integer*2 yyyy,mm,dd,hh,mn,ss
    integer timediff
    real slat,slon,glat1,glon1, slat2, slon2
    real tb(NCHAN)
    real pc_emis(NEM), emis(NEM)
    real emis_cmb(NCHAN)
    integer*2 sfc_min,sfc_max,elev
    integer*2 nku10, nka10, nka10_HS
    integer*2 nku15, nka15, nka15_HS
    integer*2 nku20, nka20, nka20_HS
    integer*2 nku25, nka25, nka25_HS
    real zen_NS
    real pia_NS, pia_MS
    real pia_NS_cmb, pia_MS_cmb(2)
    real precip_NS, precip_max_NS
    real precip2_NS, precip2_max_NS
    real precip_MS, precip_max_MS
    real precip2_MS, precip2_max_MS
    real precip_NS_cmb, precip_max_NS_cmb
    real precip_MS_cmb, precip_max_MS_cmb
    real precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF
    integer qual_GPROF
    real ts, t2m, tqv
    real hs, ps
    integer*2 p_prof(NLEV)
    real h_prof(NLEV), t_prof(NLEV), qv_prof(NLEV)
    integer*2 z_ku(88)
    integer*2 z_ka(88)
end type strm_record

type(strm_record) prof

! Added by UTSUMI

integer             ios

!--------------

integer i, j, k, m, nbad, bytes, irec
integer :: nrec_all=0, ndbf=0, nf=0, nrain=0
integer i1, iclutter
real thick, qv, delp, kPa, abstemp, densair, rmsd, ht, z_ku, z_ka, hs1
character(len=32) pdate

real ave_emis(NEM), std_emis(NEM)
real ave_pc(NEM), std_pc(NEM)
!real b_all(NREG+1)(NEM)  ! regression coeffs for all NTBREG
real b_all(NEM,NREG+1)   ! regression coeffs for all NTBREG
real tmp(NREG+1), psum
!real pc_e(NEM), e(NEM), u(NEM)(NEM)
real pc_e(NEM), e(NEM), u(NEM,NEM)
character(len=512) svar
integer id, k1, k2, kt
real coeff_disc(NEM), disc, disc_mid
real x1,x2,y1,y2,z1,z2

real pc_e_min(NEM), pc_e_max(NEM)
!real pcmin, pcmax, pc_range(NEM)(NPCHIST)(2), pc_lo, pc_hi
real pcmin, pcmax, pc_range(2,NPCHIST,NEM), pc_lo, pc_hi
integer idx, idxn(NEM_USE), found, idx2

integer ndb_files
integer p, pow2(NEM_USE), ndb_files_used

integer,allocatable,dimension(:) :: nrec_dbfile
real gmi_freq(13) 
data gmi_freq/10.7,10.7,18.6,18.6,21.65,36.5,36.5,89.0,89.0,166.0,166.0,186.31,190.31/


ndb_files= NPCHIST**NEM_USE
allocate( nrec_dbfile(ndb_files) )
nrec_dbfile = 0

!bytes= sizeof(strm_record)
bytes= 1334

!print *,"record length=",bytes
!printf("N of EPC used=%d  N of EPC histogram bins=%d   N of DB files= %d\n", NEM_USE, NPCHIST, ndb_files);
!printf("N of EPC used=%d  N of EPC histogram bins=%d   N of DB files= %d\n", NEM_USE, NPCHIST, ndb_files);

!nrec_dbfile= (int *) malloc(ndb_files*sizeof(int));
!for (k= 0; k< ndb_files; k++) nrec_dbfile[k]= 0;

do k= 0,NEM_USE-1
    pow2(k)= NPCHIST**k
end do


!Format of PC_MIN_MAX file (11 lines)
!EPC     min        max       bin0_min   bin0_max    bin1_min    bin1_max      then same for bins 2-9
!
!  0 -199.35802  320.60181    -60.01204  -22.19836  -22.19836  -15.78338  -15.78338  -11.00258  -11.00258   -5.36850   -5.36850    0.60104    0.60104    6.83276    6.83276   13.54176   13.54176   20.40424   20.40424   29.00614   29.00614   45.48376
!  1 -182.16125  145.55750    -30.19480   -7.04496   -7.04496   -4.60922   -4.60922   -2.87756   -2.87756   -1.49092   -1.49092   -0.20506   -0.20506    1.25172    1.25172    3.15748    3.15748    6.05210    6.05210   10.28600   10.28600   24.38654
!  2   -1.76337    1.15656     -1.40294   -0.96568   -0.96568   -0.66566   -0.66566   -0.24926   -0.24926    0.19294    0.19294    0.24864    0.24864    0.28968    0.28968    0.34658    0.34658    0.40812    0.40812    0.47638    0.47638    0.65992
!  3   -0.95563    1.32334     -0.40594   -0.12074   -0.12074   -0.06542   -0.06542   -0.04040   -0.04040   -0.02248   -0.02248   -0.00452   -0.00452    0.01480    0.01480    0.03384    0.03384    0.05490    0.05490    0.08490    0.08490    0.51144
!  4   -0.55545    0.41323     -0.17628   -0.02770   -0.02770   -0.01578   -0.01578   -0.00886   -0.00886   -0.00346   -0.00346    0.00146    0.00146    0.00620    0.00620    0.01114    0.01114    0.01650    0.01650    0.02282    0.02282    0.17140
!  5   -0.32075    0.38067     -0.08374   -0.01780   -0.01780   -0.01238   -0.01238   -0.00770   -0.00770   -0.00346   -0.00346    0.00026    0.00026    0.00384    0.00384    0.00794    0.00794    0.01378    0.01378    0.02476    0.02476    0.14406
!  6   -0.48772    0.87471     -0.19236   -0.01520   -0.01520   -0.00834   -0.00834   -0.00482   -0.00482   -0.00204   -0.00204    0.00056    0.00056    0.00332    0.00332    0.00668    0.00668    0.01142    0.01142    0.02162    0.02162    0.46860
!  7   -0.47112    0.42215     -0.06606   -0.00920   -0.00920   -0.00606   -0.00606   -0.00412   -0.00412   -0.00250   -0.00250   -0.00088   -0.00088    0.00094    0.00094    0.00332    0.00332    0.00692    0.00692    0.01330    0.01330    0.09990
!  8   -0.35452    0.40279     -0.09296   -0.00660   -0.00660   -0.00374   -0.00374   -0.00234   -0.00234   -0.00128   -0.00128   -0.00030   -0.00030    0.00070    0.00070    0.00186    0.00186    0.00330    0.00330    0.00562    0.00562    0.07870
!  9   -0.99769    0.31942     -0.09440   -0.00514   -0.00514   -0.00288   -0.00288   -0.00172   -0.00172   -0.00090   -0.00090   -0.00018   -0.00018    0.00056    0.00056    0.00138    0.00138    0.00248    0.00248    0.00476    0.00476    0.06016
! 10   -0.18611    0.44873     -0.03368   -0.00364   -0.00364   -0.00224   -0.00224   -0.00134   -0.00134   -0.00060   -0.00060    0.00002    0.00002    0.00066    0.00066    0.00138    0.00138    0.00220    0.00220    0.00344    0.00344    0.04300

!!  Read EPC 10-bin histogram file */
!!  This is used to compute the database indices from the EPC structure */
!
!svar="PC_MIN_MAX_10_no_overlap.txt"
!!svar="PC_MIN_MAX_10_2pc_overlap.txt"
!open(10,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
!pc_range = -9999.
!do i = 1,NEM
!    read(10,*) k, pcmin, pcmax, (pc_range(:,j,i), j=1,NPCHIST)
!    pc_range(1,1,i)       = pcmin-100.0
!    pc_range(2,NPCHIST,i) = pcmax+100.0
!    !print *,j,pc_range(1,:,i),pc_range(2,:,i)
!
!end do
!close(10)


!  Read EPC coefficients file 
!  This calculates EPC terms from the TB combinations 

svar="coef_pc.txt"
open(11,file=svar,status="old")
!print '("Opened ",A)', trim(svar)
b_all = -9999.
do i = 1,NREG+1
    read(11,*) k,(b_all(j,i), j=1,NEM)
    !print *,i,b_all(:,i)
end do
close(11)


!  Read eigenvalues matrix file
!  This reconstructs the emissivity vector from the EPC

svar="wmatrix.txt"
open(12,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
u = -9999.
do i = 1,NEM
    read(12,*) (u(j,i), j=1,NEM)
    !print *,i,u(:,i)
end do
close(12)


!  Read ave emissivity file

svar="ave_emis.txt"
open(13,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
ave_emis = -9999.
do i = 1,NEM
    read(13,*) k,ave_emis(i), std_emis(i)
    !print *,i,ave_emis(i)
end do
close(13)


!!  Read ave EPC file
!svar="ave_pc.txt"
!open(14,file=svar,status="old")
!print '("")'
!print '("Opened ",A)', trim(svar)
!ave_pc = -9999.
!do i = 1,NEM
!    read(14,*) k,ave_pc(i), std_pc(i)
!    !print *,i,ave_pc(i), std_pc(i)
!end do
!close(14)



!---------------------------------------------------------------------------
!print *,trim(idbname)
open(15, file=trim(idbname), form="unformatted", status="old", access="direct", recl=byte)
nrain= 0

ios = 0    ! For gfortran
irec=0
do irec =1,nrec
    !irec=irec+1
    read(15,rec=irec, iostat=ios) prof
    !if (ios.ne.0) exit

    !---------------------------------------------------
    ! quality control on each DB record 

    hs1= -99.99
    do i=2,NLEV
        if (prof%h_prof(i).le.0.0) then
            hs1 = prof%h_prof(i-1)
            exit   ! necessary? 
        end if     
    end do

    if ( hs1.lt.0.0) hs1= prof%h_prof(NLEV)
    if ( prof%hs .gt. hs1 )then
        print *,"hs larger than lowest geopotential ht"
        print *,"hs=",prof%hs
        print *,"h_lowest=", hs1
        stop        
    end if
 
    nbad= 0

    do i=1,9
        if ((prof%tb(i).lt.20.0).or.(prof%tb(i).gt.350)) nbad=nbad+1
    end do
    if ( nbad .gt.0) cycle
    if ( (prof%precip_NS .lt.0.0).and.(prof%precip_NS_cmb .lt.0.0) ) cycle
    if ( prof%sfc_class  .lt. 0 ) cycle
    if (( abs(prof%glat1) .gt. 90).and.(abs(prof%glon1) .gt. 180 )) cycle


    if (( (prof%yyyy .lt. 2014) .or. (prof%yyyy .gt. 2020 )) .or. &
      ( (prof%mm .lt. 1) .or. (prof%mm .gt. 12 )).or. &
      ( (prof%dd .lt. 1) .or. (prof%dd .gt. 31 )).or. &
      ( (prof%hh .lt. 0) .or. (prof%hh .gt. 23 )).or. &
      ( (prof%mn .lt. 0) .or. (prof%mn .gt. 59 )).or. &
      ( (prof%ss .lt. 0) .or. (prof%ss .gt. 59 )))then

      print *,"Bad date or time", prof%yyyy, prof%mm, prof%dd, prof%hh, prof%mn, prof%ss
      cycle

    end if

    !---------------------------------------------------*/


    ! print only the records where R > 10 

    !if ( prof%precip_NS_cmb .lt. 10 ) cycle
    !if ( prof%precip_NS_cmb .lt. 0 ) cycle

    !print *,""
    !write(pdate, "(I4,'/',I2.2,'/',I2.2, ' ',I2.2,':',I2.2,':',I2.2)") prof%yyyy, prof%mm, prof%dd, prof%hh, prof%mn, prof%ss
    !print '(I7, I7.6, I5, I5, I5,I5," ", A )', irec, prof%rev, prof%i_S1, prof%j_S1, prof%i_NS, prof%j_NS, pdate
    !print '(f8.3," ", f8.3)' , prof%glat1, prof%glon1
    !print '(I4, I5, I5, I5 )', prof%sfc_class, prof%sfc_min, prof%sfc_max, prof%elev
    !print '(f9.3, f9.3, f9.3 )', prof%ts, prof%t2m, prof%tqv
    !print '(f8.3 )', prof%precip_NS
    !print '(f8.3 )', prof%precip_NS_cmb
    !do i= 1,9
    !    write(*,fmt='(f7.2 )',advance='no') prof%tb(i)
    !end do
    !print *,""


    ! precip
    a1precip_NS(irec) = prof%precip_NS


    ! set up regression variables
    tmp(1)= 1.0
    kt= 1
    do k= 1,NTBREG
        tmp(kt+1)= prof%tb(k)
        kt=kt+1
        do k1= k,NTBREG
            tmp(kt+1)= prof%tb(k)*prof%tb(k1)
            kt=kt+1
        end do
    end do
    tmp(kt+1)= (prof%tb(1)-prof%tb(2))/(prof%tb(1)+prof%tb(2))
    kt = kt+1
    tmp(kt+1)= (prof%tb(3)-prof%tb(4))/(prof%tb(3)+prof%tb(4))
    kt = kt+1
    tmp(kt+1)= (prof%tb(6)-prof%tb(7))/(prof%tb(6)+prof%tb(7))
    kt = kt+1



    ! emissivity PC
    !write(*,fmt="(A)",advance="no") "EPC= "
    do k= 1,NEM
        psum= 0.0
        do k2= 1,NREG+1
            psum = psum+ b_all(k,k2)*tmp(k2)
            !print *, k,k2,b_all(k,k2),tmp(k2)
        end do
        pc_e(k)= psum
        prof%pc_emis(k)= pc_e(k)
        !write(*, fmt='(f8.4 )',advance="no") pc_e(k)
    end do
    !print *,""

    a2pc(:,irec) = pc_e


    ! reconstruct emissivites from PC
    do k= 1,NEM
        psum= 0.0

        do k2= 1,NEM
            psum = psum + u(k2,k) * pc_e(k2)
        end do
        e(k)= psum + ave_emis(k)
        prof%emis(k)= e(k)
        !if ( k .eq. 1 ) write(*,fmt="(A)",advance="no") "emis= "
        !if ( k .eq. NEM-1 ) write(*,fmt="(A)",advance="no") " Ts= "
        !if ( k .eq. NEM ) write(*,fmt="(A)",advance="no") "  Vap= "
        !write(*,"(f8.3 )",advance="no"), e(k)
    end do
    !print *,""


    ! Profile of geopot height, pressure, temperature, specific humidity */
    ! for (j= 0; j< NLEV; j++) {
    !  printf("%2d %8.2f %7.2f %6.2f %11.8f\n",j,prof.h_prof[j],prof.p_prof[j],prof.t_prof[j], prof.qv_prof[j]);
    !}*/

    !! DPR radar reflectivity profiles
    !do j= 1,NLEV_NS
    !    ht= 0.25*(NLEV_NS-j)
    !    write(*,fmt='(I3," ",f6.3)',advance="no") j, ht
    !    z_ku= 0.01*prof%z_ku(j)
    !    iclutter= 0
    !    if (( z_ku .lt. 0.0) .and. (z_ku .gt. -99.0 )) then
    !        iclutter= 1
    !        z_ku = z_ku*(-1.0)
    !    end if
    !    write(*,fmt='("  Ku=",I1," ",f6.2)',advance="no") iclutter, z_ku
    !    if (( prof%j_NS .lt. 12).or.(prof%j_NS .gt. 36 )) then
    !        print *,""
    !        cycle
    !    end if
    !    z_ka= 0.01*prof%z_ka(j)
    !    iclutter= 0
    !    if ( (z_ka .lt. 0.0) .and. (z_ka .gt. -99.0 ))then
    !        iclutter= 1
    !        z_ka = z_ka*(-1.0)
    !    end if
    !    write(*,fmt='("  Ka= ",I1," ",f6.2)') iclutter, z_ka
    !end do


    if ( prof%precip_NS_cmb .gt. 0.0 ) nrain=nrain+1

    ! next record */
    !print *,"irec=",irec
end do
close(15)

!ndbf = sum(nrec_dbfile)
!!write(*,'(A,"  N=",I7,"  Nraining=",I7, "  N_dbfiles=",I7)') trim(idbname),irec, nrain, ndbf 
return 
END SUBROUTINE read_db_pc_pr

END MODULE f_read_db
