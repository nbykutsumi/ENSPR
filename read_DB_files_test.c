/*
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

*/
/*------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <time.h>

#define NLEV_NS 88  /* DPR 250-m bins */
#define NLEV 42  /* from MERRA2  42 levels */
#define NCHAN 13  /* GMI TB */

#define NEM 11  /* 11= first nine emis then Ts then total column vapor */
#define NTBREG 9
#define NREG 57

#define NEM_USE 4
#define NPCHIST 10


/* record structure (1334 bytes) */

typedef struct {
     int satid, rev, rev2;
     short SC_orientation,i_S1,j_S1,i_NS,j_NS,sfc_class;
     short yyyy,mm,dd,hh,mn,ss;
     int timediff;
     float slat,slon,glat1,glon1, slat2, slon2;
     float tb[NCHAN];
     float pc_emis[NEM], emis[NEM];
     float emis_cmb[NCHAN];
     short sfc_min,sfc_max,elev;
     short nku10, nka10, nka10_HS;
     short nku15, nka15, nka15_HS;
     short nku20, nka20, nka20_HS;
     short nku25, nka25, nka25_HS;
     float zen_NS;
     float pia_NS, pia_MS;
     float pia_NS_cmb, pia_MS_cmb[2];
     float precip_NS, precip_max_NS;
     float precip2_NS, precip2_max_NS;
     float precip_MS, precip_max_MS;
     float precip2_MS, precip2_max_MS;
     float precip_NS_cmb, precip_max_NS_cmb;
     float precip_MS_cmb, precip_max_MS_cmb;
     float precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF;
     int qual_GPROF;
     float ts, t2m, tqv;
     float hs, ps;
     short p_prof[NLEV];
     float h_prof[NLEV], t_prof[NLEV], qv_prof[NLEV];
     short z_ku[88];
     short z_ka[88];
} strm_record;

strm_record prof;


main(argc,argv)
int argc;
char *argv[];
{

  FILE *fdb;
  int i, j, k, m, nbad, bytes, irec, nrec=0, nrec_all=0, ndbf=0, nf=0, nrain=0, i1, iclutter;
  float thick, qv, delp, kPa, abstemp, densair, rmsd, ht, z_ku, z_ka, hs1;
  char pdate[32];

  FILE *fpc;
  float ave_emis[NEM], std_emis[NEM];
  float ave_pc[NEM], std_pc[NEM];
  float b_all[NREG+1][NEM];  /* regression coeffs for all NTBREG */
  float tmp[NREG+1], psum;
  float pc_e[NEM], e[NEM], u[NEM][NEM];
  char svar[512];
  int id, k1, k2, kt;
  float coeff_disc[NEM], disc, disc_mid;
  float x1,x2,y1,y2,z1,z2;

  float pc_e_min[NEM], pc_e_max[NEM];
  float pcmin, pcmax, pc_range[NEM][NPCHIST][2], pc_lo, pc_hi;
  int idx, idxn[NEM_USE], found, idx2;

  int ndb_files;
  int p, pow2[NEM_USE], ndb_files_used;

  int *nrec_dbfile;
  float gmi_freq[13]= {10.7,10.7,18.6,18.6,21.65,36.5,36.5,89.0,89.0,166.0,166.0,186.31,190.31};


  int itmp;   /* NU */
  /*---------------------------------------------------------------*/

  if (argc != 2 ) exit(-1);

  ndb_files= pow(NPCHIST, NEM_USE);
  bytes= sizeof(strm_record);
  printf("record length= %d\n", bytes);
  printf("N of EPC used=%d  N of EPC histogram bins=%d   N of DB files= %d\n", NEM_USE, NPCHIST, ndb_files);

  nrec_dbfile= (int *) malloc(ndb_files*sizeof(int));
  for (k= 0; k< ndb_files; k++) nrec_dbfile[k]= 0;


  for (k= 0; k< NEM_USE; k++)
    pow2[k]= pow(NPCHIST,k);

/*
Format of PC_MIN_MAX file (11 lines)
EPC     min        max       bin0_min   bin0_max    bin1_min    bin1_max      then same for bins 2-9

  0 -199.35802  320.60181    -60.01204  -22.19836  -22.19836  -15.78338  -15.78338  -11.00258  -11.00258   -5.36850   -5.36850    0.60104    0.60104    6.83276    6.83276   13.54176   13.54176   20.40424   20.40424   29.00614   29.00614   45.48376 
  1 -182.16125  145.55750    -30.19480   -7.04496   -7.04496   -4.60922   -4.60922   -2.87756   -2.87756   -1.49092   -1.49092   -0.20506   -0.20506    1.25172    1.25172    3.15748    3.15748    6.05210    6.05210   10.28600   10.28600   24.38654 
  2   -1.76337    1.15656     -1.40294   -0.96568   -0.96568   -0.66566   -0.66566   -0.24926   -0.24926    0.19294    0.19294    0.24864    0.24864    0.28968    0.28968    0.34658    0.34658    0.40812    0.40812    0.47638    0.47638    0.65992 
  3   -0.95563    1.32334     -0.40594   -0.12074   -0.12074   -0.06542   -0.06542   -0.04040   -0.04040   -0.02248   -0.02248   -0.00452   -0.00452    0.01480    0.01480    0.03384    0.03384    0.05490    0.05490    0.08490    0.08490    0.51144 
  4   -0.55545    0.41323     -0.17628   -0.02770   -0.02770   -0.01578   -0.01578   -0.00886   -0.00886   -0.00346   -0.00346    0.00146    0.00146    0.00620    0.00620    0.01114    0.01114    0.01650    0.01650    0.02282    0.02282    0.17140 
  5   -0.32075    0.38067     -0.08374   -0.01780   -0.01780   -0.01238   -0.01238   -0.00770   -0.00770   -0.00346   -0.00346    0.00026    0.00026    0.00384    0.00384    0.00794    0.00794    0.01378    0.01378    0.02476    0.02476    0.14406 
  6   -0.48772    0.87471     -0.19236   -0.01520   -0.01520   -0.00834   -0.00834   -0.00482   -0.00482   -0.00204   -0.00204    0.00056    0.00056    0.00332    0.00332    0.00668    0.00668    0.01142    0.01142    0.02162    0.02162    0.46860 
  7   -0.47112    0.42215     -0.06606   -0.00920   -0.00920   -0.00606   -0.00606   -0.00412   -0.00412   -0.00250   -0.00250   -0.00088   -0.00088    0.00094    0.00094    0.00332    0.00332    0.00692    0.00692    0.01330    0.01330    0.09990 
  8   -0.35452    0.40279     -0.09296   -0.00660   -0.00660   -0.00374   -0.00374   -0.00234   -0.00234   -0.00128   -0.00128   -0.00030   -0.00030    0.00070    0.00070    0.00186    0.00186    0.00330    0.00330    0.00562    0.00562    0.07870 
  9   -0.99769    0.31942     -0.09440   -0.00514   -0.00514   -0.00288   -0.00288   -0.00172   -0.00172   -0.00090   -0.00090   -0.00018   -0.00018    0.00056    0.00056    0.00138    0.00138    0.00248    0.00248    0.00476    0.00476    0.06016 
 10   -0.18611    0.44873     -0.03368   -0.00364   -0.00364   -0.00224   -0.00224   -0.00134   -0.00134   -0.00060   -0.00060    0.00002    0.00002    0.00066    0.00066    0.00138    0.00138    0.00220    0.00220    0.00344    0.00344    0.04300 
*/

  /*  Read EPC 10-bin histogram file */
  /*  This is used to compute the database indices from the EPC structure */

  sprintf(svar,"%s", "PC_MIN_MAX_10_no_overlap.txt");
  /*sprintf(svar,"%s", "PC_MIN_MAX_10_2pc_overlap.txt");*/
  if ((fpc= fopen(svar,"r")) == 0) {
    printf("Unable to open %s\n", svar);
    exit(-1);
  }
  printf("Opened %s\n", svar);
  for (i=0; i< NEM; i++) {
    fscanf(fpc,"%d %f %f ", &k, &pcmin, &pcmax);
    printf("%2d %12.6f %12.6f ", k, pcmin, pcmax);
    for (j=0; j< NPCHIST; j++) {
      fscanf(fpc,"%f %f ", &pc_range[i][j][0], &pc_range[i][j][1]);
      if ( j == 0 ) pc_range[i][j][0]= pcmin - 100.0;
      if ( j == NPCHIST-1 ) pc_range[i][j][1]= pcmax + 100.0;
      printf("%12.6f %12.6f ", pc_range[i][j][0], pc_range[i][j][1]);
    }
    fscanf(fpc,"\n");
    printf("\n");
  }
  fclose(fpc);

  for (j=0; j< NPCHIST; j++) {
      printf("%d %12.6f %12.6f \n", j,pc_range[0][j][0], pc_range[0][j][1]);
  } 
  printf("--------------\n");




  /*  Read EPC coefficients file */
  /*  This calculates EPC terms from the TB combinations */

  sprintf(svar,"%s", "coef_pc.txt");
  if ((fpc= fopen(svar,"r")) == 0) {
    printf("Unable to open %s\n", svar);
    exit(-1);
  }
  printf("Opened %s\n", svar);
  for (i=0; i<= NREG; i++) {
    fscanf(fpc,"%d ", &k);
    for (j=0; j< NEM; j++) {
      fscanf(fpc,"%f ", &b_all[i][j]);
      printf("%12.6f ", b_all[i][j]);
    }
    fscanf(fpc,"\n");
    printf("\n");
  }
  fclose(fpc);

  /*  Read eigenvalues matrix file */
  /*  This reconstructs the emissivity vector from the EPC */

  sprintf(svar,"%s", "wmatrix.txt");
  if ((fpc= fopen(svar,"r")) == 0) {
    printf("Unable to open %s\n", svar);
    exit(-1);
  }
  printf("Opened %s\n", svar);
  for (i = 0; i < NEM; i++) {
    for (j = 0; j < NEM; j++) fscanf(fpc,"%f ",&u[i][j]);
    fscanf(fpc,"\n");
    for (j = 0; j < NEM; j++) printf("%12.6f ",u[i][j]);
    printf("\n");
  }
  fclose(fpc);

  /*  Read ave emissivity file */

  sprintf(svar,"%s", "ave_emis.txt");
  if ((fpc= fopen(svar,"r")) == 0) {
    printf("Unable to open %s\n", svar);
    exit(-1);
  }
  printf("Opened %s\n", svar);
  for (i = 0; i < NEM; i++) {
    fscanf(fpc,"%d %f %f\n", &j, &ave_emis[i], &std_emis[i]);
    printf("%2d %12.6f\n", i, ave_emis[i]);
  }
  fclose(fpc);

  /*  Read ave EPC file */

  sprintf(svar,"%s", "ave_pc.txt");
  if ((fpc= fopen(svar,"r")) == 0) {
    printf("Unable to open %s\n", svar);
    exit(-1);
  }
  printf("Opened %s\n", svar);
  for (i = 0; i < NEM; i++) {
    fscanf(fpc,"%d %f %f\n", &j, &ave_pc[i], &std_pc[i]);
    printf("%2d %12.6f %12.6f\n", i, ave_pc[i], std_pc[i]);
  }
  fclose(fpc);

  /*---------------------------------------------------------------------------*/

  if ((fdb = fopen(argv[1],"r")) == 0) {
    printf("Unable to open %s\n", argv[1]);
    exit(-1);
  }

  nrec= nrain= 0;
  itmp = -1;  /* NU */
  while ( fread(&prof,bytes,1,fdb) == 1 ) {
    itmp++;
    /*---------------------------------------------------*/
    /* quality control on each DB record */

    hs1= -99.99;
    for (i= 1; i< NLEV; i++) {
      if ( prof.h_prof[i] < 0.0 ) {
        hs1= prof.h_prof[i-1];
        break;
      }
    }
    if ( hs1 < 0.0 ) hs1= prof.h_prof[NLEV-1];
    if ( prof.hs > hs1 ) {
      printf("hs larger than lowest geopotential ht  hs=%f h_lowest=%f", prof.hs, hs1);
      exit(1);
    }

    nbad= 0;
    for (i= 0; i< 9; i++) if ( prof.tb[i] < 20.0 || prof.tb[i] > 350 ) nbad++;
    if ( nbad > 0 ) continue;

    if ( prof.precip_NS < 0.0 && prof.precip_NS_cmb < 0.0 ) continue;
    if ( prof.sfc_class < 0 ) continue;
    if ( fabs(prof.glat1) > 90 || fabs(prof.glon1) > 180 ) continue;

    if (( prof.yyyy < 2014 || prof.yyyy > 2020 ) ||
      ( prof.mm < 1 || prof.mm > 12 ) ||
      ( prof.dd < 1 || prof.dd > 31 ) ||
      ( prof.hh < 0 || prof.hh > 23 ) ||
      ( prof.mn < 0 || prof.mn > 59 ) ||
      ( prof.ss < 0 || prof.ss > 59 )) {
      printf("Bad date or time  %d %d %d %d %d %d\n", prof.yyyy, prof.mm, prof.dd, prof.hh, prof.mn, prof.ss);
      continue;
    }
    /*---------------------------------------------------*/
    if ( itmp !=930) continue;

    /* print only the records where R > 10 */
    /*if ( prof.precip_NS_cmb < 10 ) continue;*/
    /*if ( prof.precip_NS_cmb < 0 ) continue; */


    printf("\n\n");
    sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d", prof.yyyy, prof.mm, prof.dd, prof.hh, prof.mn, prof.ss);
    printf("%6d %06d %4d %4d %4d %4d %s ", nrec, prof.rev, prof.i_S1, prof.j_S1, prof.i_NS, prof.j_NS, pdate);
    printf("%8.3f %8.3f ", prof.glat1, prof.glon1);
    printf("%3d %4d %4d %5d ", prof.sfc_class, prof.sfc_min, prof.sfc_max, prof.elev);
    printf("%8.3f %8.3f %8.3f ", prof.ts, prof.t2m, prof.tqv);
    printf("%8.3f ", prof.precip_NS);
    printf("%8.3f ", prof.precip_NS_cmb);
    for (i= 0; i< 9; i++) printf("%6.2f ", prof.tb[i]);
    printf("\n");

    /* set up regression variables */
    tmp[0]= 1.0;
    kt= 0;
    for (k= 0; k< NTBREG; k++) {
      tmp[kt+1]= prof.tb[k];
      kt++;
      for (k1= k; k1< NTBREG; k1++) {
        tmp[kt+1]= prof.tb[k]*prof.tb[k1];
        kt++;
      }
    }
    printf("last three\n");
    printf("0 %8.4f\n", prof.tb[0]);
    printf("1 %8.4f\n", prof.tb[1]);
    printf("2 %8.4f\n", prof.tb[2]);
    printf("3 %8.4f\n", prof.tb[3]);
    printf("5 %8.4f\n", prof.tb[5]);
    printf("6 %8.4f\n", prof.tb[6]);
    
    tmp[kt+1]= (prof.tb[0]-prof.tb[1])/(prof.tb[0]+prof.tb[1]);
    kt++;
    tmp[kt+1]= (prof.tb[2]-prof.tb[3])/(prof.tb[2]+prof.tb[3]);
    kt++;
    tmp[kt+1]= (prof.tb[5]-prof.tb[6])/(prof.tb[5]+prof.tb[6]);
    kt++;

    /* NU */
    printf("tbcomb\n");
    for (k2= 0; k2< NREG+1; k2++) {
    printf("%d ", k2);
    printf("%8.4f \n", tmp[k2]);
    }
    
    printf("\n");

    /* emissivity PC */
    printf("EPC= ");
    for (k= 0; k< NEM; k++) {
      psum= 0.0;
      for (k2= 0; k2< NREG+1; k2++) psum+= b_all[k2][k]*tmp[k2];
      pc_e[k]= psum;
      prof.pc_emis[k]= pc_e[k];
      printf("%8.4f ", pc_e[k]);
    }
    printf("\n");


    /* reconstruct emissivites from PC */
    for (k= 0; k< NEM; k++) {
      psum= 0.0;
      for (k2= 0; k2< NEM; k2++) psum+= u[k][k2] * pc_e[k2];
      e[k]= psum + ave_emis[k];
      prof.emis[k]= e[k];
      if ( k == 0 ) printf("emis= ");
      if ( k == NEM-2 ) printf("  Ts= ");
      if ( k == NEM-1 ) printf("  Vap= ");
      printf("%8.3f ", e[k]);
    }
    printf("\n");

    /* Profile of geopot height, pressure, temperature, specific humidity */
    /*for (j= 0; j< NLEV; j++) {
      printf("%2d %8.2f %7.2f %6.2f %11.8f\n",j,prof.h_prof[j],prof.p_prof[j],prof.t_prof[j], prof.qv_prof[j]);
    }*/

    /* DPR radar reflectivity profiles */
    for (j= 0; j< NLEV_NS; j++) {
      ht= 0.25*(NLEV_NS-j-1);
      printf("%3d %6.3f ", j, ht); 
      z_ku= 0.01*prof.z_ku[j]; 
      iclutter= 0;
      if ( z_ku < 0.0 && z_ku > -99.0 ) {
        iclutter= 1;
        z_ku*= -1.0;
      }
      printf("Ku=%1d %6.2f ", iclutter, z_ku); 
      if ( prof.j_NS < 12 || prof.j_NS > 36 ) {
        printf("\n");
        continue;
      }
      z_ka= 0.01*prof.z_ka[j]; 
      iclutter= 0;
      if ( z_ka < 0.0 && z_ka > -99.0 ) {
        iclutter= 1;
        z_ka*= -1.0;
      }
      printf("  Ka=%1d %6.2f \n", iclutter, z_ka); 
    }

    /* Calculate database bin indices that would be populated with this observation */

    /* first index moving fastest */
    for (i1= 0; i1< ndb_files; i1++) {
        /* 1st PC   -------*/
        k= 0;
        idxn[k]= i1 % NPCHIST;
        found= 0;
        for (j= 0; j < NPCHIST; j++) {
          pc_lo= pc_range[k][j][0];
          pc_hi= pc_range[k][j][1];

          /*
          printf("AA  i1=%d k=%d j=%d idxn=%d found=%d     pc_e(k)=%8.2f  lo=%8.2f  hi=%8.2f\n", i1, k, j, idxn[k], found, pc_e[k],pc_lo, pc_hi);
          */


          /* test
          if ( pc_e[k] > pc_lo && pc_e[k] <= pc_hi && j == idxn[k] ){
          found=1;
          }
          printf("i1=%d  k=%d  j=%d  idxn=%d  found=%d\n", i1, k, j, idxn[k], found);
          if (found == 1) {exit(0);}
          */



          /*printf("%2d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_e[k], pc_lo, pc_hi, j , idxn[k] );*/
          if ( pc_e[k] > pc_lo && pc_e[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
          /* NU
          Why is the condition (j==idxn[k]) necessary?
          */

        }
        if ( found != 1 ) continue;
        /*
        printf("AA Found!!  i1=%d k=%d j=%d idxn=%d found=%d\n", i1, k, j, idxn[k], found); 
        */ 


        /* 2nd ~ 2nd-last PC   -------*/
        for (k= 1; k< NEM_USE-1; k++) {
          p= pow2[k];
          idxn[k]= (i1/p) % NPCHIST;


          found= 0;
          for (j= 0; j < NPCHIST; j++) {
            pc_lo= pc_range[k][j][0];
            pc_hi= pc_range[k][j][1];
            /*
            printf("BB  i1=%d k=%d j=%d idxn=%d found=%d     pc_e(k)=%8.2f  lo=%8.2f  hi=%8.2f\n", i1, k, j, idxn[k], found, pc_e[k],pc_lo, pc_hi);
            */

            /*printf("%d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_e[k], pc_lo, pc_hi, j , idxn[k] );*/
            if ( pc_e[k] > pc_lo && pc_e[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
          }
          if ( found != 1 ) goto bypass;
          /*
          printf("BB found!!  i1=%d k=%d j=%d idxn=%d found=%d     pc_e(k)=%8.2f  lo=%8.2f  hi=%8.2f\n", i1, k, j, idxn[k], found, pc_e[k],pc_lo, pc_hi);
          */
        }

        /* last PC   -------*/
        k= NEM_USE-1;
        idxn[k]= i1/pow2[k];
        /*printf("i1=%d  k=%d  idxn(k)=%d\n", i1, k, idxn[k]);*/

        found= 0;
        for (j= 0; j < NPCHIST; j++) {
          pc_lo= pc_range[k][j][0];
          pc_hi= pc_range[k][j][1];
          /*
          printf("CC  i1=%d k=%d j=%d idxn=%d found=%d\n", i1, k, j, idxn[k], found);
          */

          /*printf("%d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_e[k], pc_lo, pc_hi, j, idxn[k] );*/
          if ( pc_e[k] > pc_lo && pc_e[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
        }
        if ( found != 1 ) continue;

        /*
        printf("CC found!! i1=%d k=%d j=%d idxn=%d found=%d\n", i1, k, j, idxn[k], found); 
        printf("nrec=%d\n", nrec); 
        */ 

        idx= 0;
        for (k= 0; k< NEM_USE; k++)
          idx+= idxn[k]*pow2[k];

        for (k= 0; k< NEM_USE; k++) {
          printf("%d  ", idxn[k]);
        }
        printf("  index=%d\n", idx);

        if ( idx < 0 || idx >= ndb_files ) {
          printf("Exceeded max files %d\n", idx);
          exit(1);
        }

        nrec_dbfile[idx]++;
        /*printf("nrec_dbfile[idx]=%d\n", nrec_dbfile[idx]);*/ 

        /* NU
        What is nrec_dbfile? What is the difference from nrec?
        */


        /*printf("%5d  N=%d\n", idx, nrec_dbfile[idx]);*/

bypass:;
      }

      nrec++;
      if ( prof.precip_NS_cmb > 0.0 ) nrain++;

      /* NU
      Why does the searchng continue until the end of loop?
      In case the overlapped DB is used?
      */


    }  /* next record */

    fclose(fdb);


    for (k= 0; k< ndb_files; k++) if ( nrec_dbfile[k] > 0 ) ndbf++;
    printf("%s N=%d  Nraining=%d N_dbfiles=%d\n", argv[1], nrec, nrain, ndbf);


  exit(0);

}

