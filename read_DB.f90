program read_DB
implicit none


integer,parameter :: NLEV_NS = 88  ! DPR 250-m bins
integer,parameter :: NLEV = 42  ! from MERRA2  42 levels
integer,parameter :: NCHAN =13  ! GMI TB
integer,parameter :: NEM = 11  ! 11= first nine emis then Ts then total column vapor
integer,parameter :: NTBREG = 9
integer,parameter :: NREG = 57
integer,parameter :: NEM_USE = 4
integer,parameter :: NPCHIST = 10
integer,parameter :: BYTE = 1334  ! UTSUMI

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
integer :: nrec=0, nrec_all=0, ndbf=0, nf=0, nrain=0
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


! argc, argv
character(len=512) :: ifilename




ndb_files= NPCHIST**NEM_USE
allocate( nrec_dbfile(ndb_files) )
nrec_dbfile = 0

!bytes= sizeof(strm_record)
bytes= 1334

print *,"record length=",bytes
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

!  Read EPC 10-bin histogram file */
!  This is used to compute the database indices from the EPC structure */

svar="PC_MIN_MAX_10_no_overlap.txt"
!svar="PC_MIN_MAX_10_2pc_overlap.txt"
open(10,file=svar,status="old")
print '("")'
print '("Opened ",A)', trim(svar)
pc_range = -9999.
do i = 1,NEM
    read(10,*) k, pcmin, pcmax, (pc_range(:,j,i), j=1,NPCHIST)
    pc_range(1,1,i)       = pcmin-100.0
    pc_range(2,NPCHIST,i) = pcmax+100.0
    !print *,j,pc_range(1,:,i),pc_range(2,:,i)

end do
close(10)


!  Read EPC coefficients file 
!  This calculates EPC terms from the TB combinations 

svar="coef_pc.txt"
open(11,file=svar,status="old")
print '("")'
print '("Opened ",A)', trim(svar)
b_all = -9999.
do i = 1,NREG+1
    read(11,*) k,(b_all(j,i), j=1,NEM)
    print *,i,b_all(:,i)
end do
close(11)


!  Read eigenvalues matrix file
!  This reconstructs the emissivity vector from the EPC

svar="wmatrix.txt"
open(12,file=svar,status="old")
print '("")'
print '("Opened ",A)', trim(svar)
u = -9999.
do i = 1,NEM
    read(12,*) (u(j,i), j=1,NEM)
    !print *,i,u(:,i)
end do
close(12)


!  Read ave emissivity file

svar="ave_emis.txt"
open(13,file=svar,status="old")
print '("")'
print '("Opened ",A)', trim(svar)
ave_emis = -9999.
do i = 1,NEM
    read(13,*) k,ave_emis(i), std_emis(i)
    !print *,i,ave_emis(i)
end do
close(13)


!  Read ave EPC file
svar="ave_pc.txt"
open(14,file=svar,status="old")
print '("")'
print '("Opened ",A)', trim(svar)
ave_pc = -9999.
do i = 1,NEM
    read(14,*) k,ave_pc(i), std_pc(i)
    !print *,i,ave_pc(i), std_pc(i)
end do
close(14)


!
!---------------------------------------------------------------------------
if (iargc().lt.1) then
    print *,"NO Input File"
    print *,"USAGE:  cmd IFILE"
    print *,"EXIT PROGRAM"
    stop
end if

call getarg(1, ifilename)
print *,trim(ifilename)
open(15, file=trim(ifilename), form="unformatted", status="old", access="direct", recl=byte)
nrec = 0
nrain= 0

ios = 0    ! For gfortran
irec=0
do
    irec=irec+1
    read(15,rec=irec, iostat=ios) prof
    if (ios.ne.0) exit

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

    if ( prof%precip_NS_cmb .lt. 10 ) cycle
    !if ( prof%precip_NS_cmb .lt. 0 ) cycle

    print *,""
    write(pdate, "(I4,'/',I2.2,'/',I2.2, ' ',I2.2,':',I2.2,':',I2.2)") prof%yyyy, prof%mm, prof%dd, prof%hh, prof%mn, prof%ss
    print '(I7, I7.6, I5, I5, I5,I5," ", A )', nrec, prof%rev, prof%i_S1, prof%j_S1, prof%i_NS, prof%j_NS, pdate
    print '(f8.3," ", f8.3)' , prof%glat1, prof%glon1
    print '(I4, I5, I5, I5 )', prof%sfc_class, prof%sfc_min, prof%sfc_max, prof%elev
    print '(f9.3, f9.3, f9.3 )', prof%ts, prof%t2m, prof%tqv
    print '(f8.3 )', prof%precip_NS
    print '(f8.3 )', prof%precip_NS_cmb
    do i= 1,9
        write(*,fmt='(f7.2 )',advance='no') prof%tb(i)
    end do
    print *,""


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
    write(*,fmt="(A)",advance="no") "EPC= "
    do k= 1,NEM
        psum= 0.0
        do k2= 1,NREG+1
            psum = psum+ b_all(k,k2)*tmp(k2)
            !print *, k,k2,b_all(k,k2),tmp(k2)
        end do
        pc_e(k)= psum
        prof%pc_emis(k)= pc_e(k)
        write(*, fmt='(f8.4 )',advance="no") pc_e(k)
    end do
    print *,""

    ! reconstruct emissivites from PC
    do k= 1,NEM
        psum= 0.0

        do k2= 1,NEM
            psum = psum + u(k2,k) * pc_e(k2)
        end do
        e(k)= psum + ave_emis(k)
        prof%emis(k)= e(k)
        if ( k .eq. 1 ) write(*,fmt="(A)",advance="no") "emis= "
        if ( k .eq. NEM-1 ) write(*,fmt="(A)",advance="no") " Ts= "
        if ( k .eq. NEM ) write(*,fmt="(A)",advance="no") "  Vap= "
        write(*,"(f8.3 )",advance="no"), e(k)
    end do
    print *,""


    ! Profile of geopot height, pressure, temperature, specific humidity */
    ! for (j= 0; j< NLEV; j++) {
    !  printf("%2d %8.2f %7.2f %6.2f %11.8f\n",j,prof.h_prof[j],prof.p_prof[j],prof.t_prof[j], prof.qv_prof[j]);
    !}*/

    ! DPR radar reflectivity profiles
    do j= 1,NLEV_NS
        ht= 0.25*(NLEV_NS-j)
        write(*,fmt='(I3," ",f6.3)',advance="no") j, ht
        z_ku= 0.01*prof%z_ku(j)
        iclutter= 0
        if (( z_ku .lt. 0.0) .and. (z_ku .gt. -99.0 )) then
            iclutter= 1
            z_ku = z_ku*(-1.0)
        end if
        write(*,fmt='("  Ku=",I1," ",f6.2)',advance="no") iclutter, z_ku
        if (( prof%j_NS .lt. 12).or.(prof%j_NS .gt. 36 )) then
            print *,""
            cycle
        end if
        z_ka= 0.01*prof%z_ka(j)
        iclutter= 0
        if ( (z_ka .lt. 0.0) .and. (z_ka .gt. -99.0 ))then
            iclutter= 1
            z_ka = z_ka*(-1.0)
        end if
        write(*,fmt='("  Ka= ",I1," ",f6.2)') iclutter, z_ka
    end do


    ! Calculate database bin indices that would be populated with this observation 

    ! first index moving fastest
    do i1= 1, ndb_files
        ! First PC  ---------
        k= 1
        idxn(k)= mod(i1-1, NPCHIST)+1
        found= 0
        do j= 1,NPCHIST
            pc_lo= pc_range(1,j,k)
            pc_hi= pc_range(2,j,k)
            !printf("%2d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_e[k], pc_lo, pc_hi, j , idxn[k] )



            ! test
            !if ( (pc_e(k) .gt. pc_lo) .and. (pc_e(k) .le. pc_hi) .and. (j .eq. idxn(k)) ) then
            !found =1
            !end if
            !write(*,'("i1=",I6," k=",I6, " j=",I3" idxn=",I6 ," found=",I1)') i1, k, j, idxn(k), found
            !if (found.eq.1) stop


            if ( (pc_e(k) .gt. pc_lo) .and. (pc_e(k) .le. pc_hi) .and. (j .eq. idxn(k)) ) then
                found=1
                exit
            end if
        end do
        if ( found .ne. 1 ) cycle
        !write(*,'("AA found!! i1=",I6," k=",I3," j=",I3," idxn=",I3," found=",I2)') i1,k,j,idxn(k),found

        ! 2nd ~ 2nd-last PC  ---------
        do k= 2, NEM_USE-1
          p= pow2(k-1)
          idxn(k)= mod(((i1-1)/p) , NPCHIST)+1

          found= 0
          do j= 1, NPCHIST
            pc_lo= pc_range(1,j,k)
            pc_hi= pc_range(2,j,k)
            !print *,k,"il=",il,"pc=",pc_e(k),"range=",pc_lo,pc_hi,"j=",j,"idx=",idxn(k)

            if (( pc_e(k) .gt. pc_lo).and.(pc_e(k) .le. pc_hi).and.(j .eq. idxn(k) )) then
                found=1
                exit
            end if
          end do
          if ( found .ne. 1 ) goto 10
          !write(*,'("BB found!! i1=",I6," k=",I3," j=",I3," idxn=",I3," found=",I2)') i1,k,j,idxn(k),found
        end do

        ! last PC  ---------
        k= NEM_USE
        idxn(k)= (i1-1)/pow2(k-1) +1
        !print *,"i1=",i1,"k=",k,"idxn(k)=",idxn(k)

        found= 0
        do j= 1, NPCHIST
          pc_lo= pc_range(1,j,k)
          pc_hi= pc_range(2,j,k)
          !write(*,'("CC i1=",I6," k=",I3," j=",I3," idxn=",I3," found=",I2)') i1,k,j,idxn(k),found
          if ( (pc_e(k) .gt. pc_lo) .and. (pc_e(k) .le. pc_hi) .and. (j .eq. idxn(k)) )then
            found=1
            exit
          end if  
        end do
        if ( found .ne. 1 ) cycle
        !write(*,'("CC found!! i1=",I6," k=",I3," j=",I3," idxn=",I3," found=",I2)') i1,k,j,idxn(k),found

        idx= 0
        do k=1,NEM_USE
          idx = idx + (idxn(k)-1)*pow2(k-1)
        end do
        nrec_dbfile(idx)=1

        do k=1, NEM_USE
          write(*,fmt='(I6," ")',advance="no") idxn(k)
        end do

        write(*,fmt='("  index=",I6)') idx
        if (( idx .lt. 1).or.(idx .ge. ndb_files)) then
          write(*,'("Exceeded max files", I6)') idx
          stop
        end if

!        nrec_dbfile[idx]++;
!        /*printf("%5d  N=%d\n", idx, nrec_dbfile[idx]);*/

10      continue
    end do

    nrec = nrec + 1
    if ( prof%precip_NS_cmb .gt. 0.0 ) nrain=nrain+1

    ! next record */
end do
close(15)

ndbf = sum(nrec_dbfile)
write(*,'(A,"  N=",I7,"  Nraining=",I7, "  N_dbfiles=",I7)') trim(ifilename),nrec, nrain, ndbf 

end program
