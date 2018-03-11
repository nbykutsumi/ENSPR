/*
Reads PPS 1B GMI production files

Run interp_model on the GMI file to get the model file

Performs EPC-based precip retrievals using first 9 GMI channels
Saves lat, lon, time, precip, TB as netCDF

Compile as:
gcc -fpack-struct -L/opt/local/lib -l netcdf -I /opt/local/include -o  file  file.c 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <netcdf.h>

#define NEM 11
#define NTBREG 9
#define NREG 57

#define RADEARTH 6371
#define RTD 57.29578
#define DTR 0.017453

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(-1);}

#define NCLASS 15
#define NLEV_NS 88

#define NCHAN 13  /* GMI */
#define NEM_USE 4
#define NPCHIST 10

#define PRECIP1 0  /* for precip PDF */
#define PRECIP2 100
#define DPRECIP 0.5
#define NPRECIP 200

#define ZMIN 15   /* For CCFADS */
#define ZMAX 55
#define DZ 0.5
#define NZ 80

/*-------------------------------------------*/
#define DB_MAXREC 2000   /* max records in biggest DB file */
#define N2MAX 500
#define WT2MIN 0.01
#define DB_RAINFRAC 0.1 
#define MAX_T2M_DIFF 15
#define MAX_TQV_DIFF 15

short em_compare[NEM]= {1,1,1,1,1,1,0,0,0,0,0};  /* length=NEM */
/*-------------------------------------------*/


   typedef struct {
     int satid, rev, rev2;
     short SC_orientation,i_S1,j_S1,i_NS,j_NS,sfc_class;
     short yyyy,mm,dd,hh,mn,ss;
     int timediff;
     float slat,slon,glat1,glon1, slat2, slon2;
     float tb[13];
     float pc_emis[11], emis[11];
     float emis_cmb[13];
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
     short p_prof[42];
     float h_prof[42], t_prof[42], qv_prof[42];
     short z_ku[88];
     short z_ka[88];
   } strm_record;


strm_record prof, prof_idx[DB_MAXREC];


int *tmp1;
int *tmp2;
float *tmp3;
float *tmp7;

int max (int *a, int n, int i, int j, int k) {
    int m = i;
    if (j < n && a[j] > a[m]) {
        m = j;
    }
    if (k < n && a[k] > a[m]) {
        m = k;
    }
    return m;
}
 
void downheap (int *a, int *b, int n, int i) {
    while (1) {
        int j = max(a, n, i, 2 * i + 1, 2 * i + 2);
        if (j == i) {
            break;
        }
        int t = a[i];
        int t2 = b[i];
        a[i] = a[j];
        a[j] = t;
        b[i] = b[j];
        b[j] = t2;
        i = j;
    }
}
 
void heapsort2 (int *a, int *b, int n) {
    int i;
    for (i = (n - 2) / 2; i >= 0; i--) {
        downheap(a, b, n, i);
    }
    for (i = 0; i < n; i++) {
        int t = a[n - i - 1];
        int t2 = b[n - i - 1];
        a[n - i - 1] = a[0];
        a[0] = t;
        b[n - i - 1] = b[0];
        b[0] = t2;
        downheap(a, b, n - i - 1, 0);
    }
}
 




main (argc, argv)
int argc;
char *argv[];

{
   int ncid;
   int grpid1;
   int grpid1_scanstatus, grpid1_scantime, grpid1_nav;
   int grpid2;
   int grpid2_scanstatus;
   int grpid_NS, grpid_MS;

   FILE *flist;
   char tbfile[256];

   int i, j, k, m, n, i1,i2,j1, j2, k1, status, isat, jsat, ip;
   int isat2, isat2_min, jsat2, jsat2_min, itest;
   int isat_start=-1, isat_end=-1, nsat=0;
   size_t nl, ns;
   time_t t0=-1.0, t1, t2;
   struct tm tm, tm0, tm1, tm2;
   char tdate[128], pdate[128], cdate[128], outdate[32], mdate[32];
   char tdate1[128], tdate2[128];
   long secs, secs1, secs2;

   /*------------------------------ 1B.GPM.GMI -------------------------------------------*/

   /*--- S1 grid geolocation ---*/
   int lat_S1_id, lon_S1_id;
   float *lat_S1_buf, *lon_S1_buf;
   nc_type lat_S1_type, lon_S1_type;
   int lat_S1_ndims, lon_S1_ndims;
   int lat_S1_dimids[NC_MAX_VAR_DIMS], lon_S1_dimids[NC_MAX_VAR_DIMS];
   int lat_S1_natts, lon_S1_natts;

   /*--- S1 grid TB ---*/
   int tb_S1_id;
   float *tb_S1_buf;
   nc_type tb_S1_type;
   int tb_S1_ndims;
   int tb_S1_dimids[NC_MAX_VAR_DIMS];
   int tb_S1_natts;

   /*--- S2 grid geolocation ---*/
   int lat_S2_id, lon_S2_id;
   float *lat_S2_buf, *lon_S2_buf;
   nc_type lat_S2_type, lon_S2_type;
   int lat_S2_ndims, lon_S2_ndims;
   int lat_S2_dimids[NC_MAX_VAR_DIMS], lon_S2_dimids[NC_MAX_VAR_DIMS];
   int lat_S2_natts, lon_S2_natts;

   /*--- S2 grid TB ---*/
   int tb_S2_id;
   float *tb_S2_buf;
   nc_type tb_S2_type;
   int tb_S2_ndims;
   int tb_S2_dimids[NC_MAX_VAR_DIMS];
   int tb_S2_natts;

   /*--- S1 grid date ---*/
   int yyyy_S1_id, mm_S1_id, dd_S1_id, hh_S1_id, mn_S1_id, ss_S1_id;
   short *yyyy_S1_buf;
   short *mm_S1_buf;
   short *dd_S1_buf;
   short *hh_S1_buf;
   short *mn_S1_buf;
   short *ss_S1_buf;

   /*--- S1 grid sat position ---*/
   int sclat_S1_id, sclon_S1_id;
   float *sclat_S1_buf, *sclon_S1_buf;

   /*--- S1 scan quality ---*/
   int qual_S1_id;
   unsigned char *qual_S1_buf;

   /*--- S2 scan quality ---*/
   int qual_S2_id;
   unsigned char *qual_S2_buf;


   size_t nl_S1, ns_S1, nc_S1;
   size_t nl_S2, ns_S2, nc_S2;

   int FractionalGranuleNumber_S1_id;
   double *FractionalGranuleNumber_S1_buf;

   double *secs_S1_buf;
   int irev, isatname;
   char satname[32], sensor[32];
   double rev;

   int nbad, nbad1;
   float percent;

   int *pixel_index, *n_index, *ncold_index, ncold= -1;
   int npix=0;
   float t2m_MERRA2_min= 1.0E6, t2m_MERRA2_max= -1.0E6;
   float tqv_MERRA2_min= 1.0E6, tqv_MERRA2_max= -1.0E6;
   float epc_min[NEM], epc_max[NEM];

   /*------------------------------ GPROF -------------------------------------------*/

   char gprof_file[256];
   int have_gprof= -1;
   int qualityFlag_SG_id, surfacePrecipitation_SG_id, probabilityOfPrecip_SG_id;
   int sfc_SG_id, totalColumnWaterVaporIndex_SG_id, temp2mIndex_SG_id;
   int frozenPrecipitation_SG_id;
   signed char *qualityFlag_GPROF_buf;
   signed char *sfc_GPROF_buf;
   float *surfacePrecipitation_GPROF_buf, *probabilityOfPrecip_GPROF_buf;
   float *lat_GPROF_buf, *lon_GPROF_buf;
   signed char *totalColumnWaterVaporIndex_GPROF_buf;
   float *frozenPrecipitation_GPROF_buf;
   short *temp2mIndex_GPROF_buf;

   nc_type lat_SG_type, lon_SG_type;
   size_t nl_SG, ns_SG;
   int lat_SG_id, lon_SG_id, lat_SG_ndims, lat_SG_dimids[NC_MAX_VAR_DIMS], lat_SG_natts;

   float precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF;
   int tqv_GPROF, t2m_GPROF;
   int t2m_GPROF_varid, tqv_GPROF_varid;
   int sfc_GPROF_varid;

   /*---------------------------------------------------------------------------------*/

   double km, km_min;
   int ntf=0, qual;
   int iyyyy, imm, idd, ihh, imn, iss;
   float pcent, tb, tb_gmi[20];
   float glat1, glon1, glat2, glon2, slat, slon, hrs, clat, clon;
   int ilat, ilon, idx, idx1, idx2, idx3, s1s2, sfc_class, idx_db, ndb, idx_out;
   float pcent_db;

   /*--- output file ---*/
   int ncid2;
   char outfile[512];
   int dimids5[NC_MAX_VAR_DIMS];
   int dimids4[NC_MAX_VAR_DIMS];
   int dimids3[NC_MAX_VAR_DIMS];
   int dimids2[NC_MAX_VAR_DIMS];
   int scan_dimid, pix_dimid, tb_dimid, em_dimid, lev_dimid;
   char ctmp[256];
   size_t start3[3], count3[3];
   size_t start2[2], count2[2];
   int retval;
   char units[256];
   int lat_varid, lon_varid;
   int secs_varid, tb_varid, emis_varid, pc_emis_varid;
   int precip_NS_varid, precip_prob_NS_varid;
   int precip2_NS_varid, precip2_prob_NS_varid;
   int precip_MS_varid, precip_prob_MS_varid;
   int precip2_MS_varid, precip2_prob_MS_varid;
   int precip_GPROF_varid, precip_prob_GPROF_varid;
   int frozen_precip_GPROF_varid;
   int ts_NS_varid, t2m_NS_varid, tqv_NS_varid, h273_NS_varid, ku_ztop_NS_varid;
   int ts_MS_varid, t2m_MS_varid, tqv_MS_varid, h273_MS_varid, ku_ztop_MS_varid, ka_ztop_MS_varid;
   int ts_MERRA2_varid, t2m_MERRA2_varid, tqv_MERRA2_varid;
   int z_ku_NS_varid;
   int z_ku_MS_varid, z_ka_MS_varid;
   int tb1_NS_varid;
   int tb1_MS_varid;
   int h273_MERRA2_varid;
   int quality_NS_varid;
   int quality_MS_varid;

   int ishift, jshift;
   int ishift_min=1000, jshift_min=1000;
   int ishift_max=-1000, jshift_max=-1000;

   /*tm0.tm_year= 2003-1970;
   tm0.tm_mon= 1-1;
   tm0.tm_mday= 1;
   tm0.tm_hour= 0;
   tm0.tm_min= 0;
   tm0.tm_sec= 0;
   t0= timegm(&tm0);
   start_secs= (int) t0;
   tm= *gmtime(&t0);
   sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
   printf("seconds between 2003 and 1970=%d\n",t0);*/

  FILE *fpc, *fout, *fcc1, *fcc2, *fret;
  float ave_pc[NEM], std_pc[NEM];
  float ave_tb[NCHAN], std_tb[NCHAN];
  float ave_emis[NEM], std_emis[NEM];
  float b_all[NREG+1][NEM];  /* regression coeffs for all NTBREG */
  float tmp5[NREG+1], psum;
  float pc_emis[NEM], emis[NEM], u[NEM+1][NEM+1];
  char svar[512], sdir[512], outdir[512];
  int k2, kt, nlay, ilay;
  int id;

   float *tb_S1S2_buf;
   float *emis_S1S2_buf;
   float *pc_emis_S1S2_buf;

   float *precip_NS_buf;
   float *precip2_NS_buf;
   float *precip_prob_NS_buf;
   float *precip2_prob_NS_buf;
   int *qual_NS_buf;
   float *precip_MS_buf;
   float *precip2_MS_buf;
   float *precip_prob_MS_buf;
   float *precip2_prob_MS_buf;
   int *qual_MS_buf;

   FILE *fdb2, *fdb_nrain;
   char dbdir[256], dbfile[256];
   int ndb_files= pow(NPCHIST, NEM_USE);
   int p, pow2[NEM_USE], ndb_files_used, bytes;
   float pcmin, pcmax, pc_range[NEM][NPCHIST][2], pc_lo, pc_hi;
   int idxn[NEM_USE], found, kuka, idb;
   int nem_compare, nem_compare_varid;

   int ix;
   float x, rmsd, rmsd1, rmsd2, range=0.0, wt, wt0, wt2, std;
   float ht, ku_ztop, ka_ztop;
   int maxrec, nrain_db[6][2], nrain, nrec, nrec1, nrec2, n1, n2, class[NCLASS], iclass;
   float z_ku[NLEV_NS], z_ka[NLEV_NS], zlog, zlog1, zlog2, zlog3, zlin, prob_precip, prob_precip2, precip, precip2, zmax;
   int nz_ku[NLEV_NS], nz_ka[NLEV_NS], hist_precip[NPRECIP], ku_clutter, ka_clutter, iz;
   int nrmsd, nrmsd1, nrmsd2, iprecip;
   int nfound, nfound_all=0;
   float rmsd_min, rmsd_max, rmsd2_min, rmsd2_max;
   short *z_ku_1_NS_buf;
   short *z_ku_1_MS_buf, *z_ka_1_MS_buf;
   float *tb_1_NS_buf, *tb_1_MS_buf;
   short *h273_buf;

   float tlat= -999.0, tlon= -999.0;
   int mrms= 0;  /* 1= look for MRMS data near the target time  */
   int n_mrms=0;
  
   int sum1, sum2, ccfads_ku[NLEV_NS][NZ], ccfads_ka[NLEV_NS][NZ]; 
   float ccfads1[NZ], ccfads2[NZ]; 

   float ts, t2m, tqv;
   float ts_sum, t2m_sum, tqv_sum, h273_sum, precip_sum, precip2_sum;
   float wt_ts_sum, wt_t2m_sum, wt_tqv_sum, wt_h273_sum, wt_precip_sum, wt_precip2_sum;
   int n_precip_sum, n_precip2_sum;
   int n_precip_all_sum, n_precip2_all_sum;
   float ku_ztop_sum, ka_ztop_sum;
   float wt_ku_ztop_sum, wt_ka_ztop_sum;

   int *ibuf;
   short *sbuf;
   char *bbuf;
   float *fbuf, *fbuf3, *fbuf4;
   double *dbuf;
   int nl_out, nc_out;

   FILE *fmrms;
   char mrms_name[512];


   /*---- variables for MERRA2 ----*/

   int have_model= -1;
   char model_file[256];
   char mfile[256], m1dir[256], m2dir[256], m3dir[256];
   float lat, lon, temp, t, qv, h, h1, ps, hs, hs1, p_lowest, h_lowest, slp, ps2;

   float *p_buf;
   float *lat_buf, *lon_buf;
   float *h_buf, *t_buf, *qv_buf;
   float *ps_buf, *hs_buf;
   float *ts_buf, *t2m_buf, *tqv_buf;
   float *ts_epc_NS_buf, *t2m_epc_NS_buf, *tqv_epc_NS_buf;
   short *h273_epc_NS_buf;
   short *ku_ztop_epc_NS_buf;
   float *ts_epc_MS_buf, *t2m_epc_MS_buf, *tqv_epc_MS_buf;
   short *h273_epc_MS_buf;
   short *ku_ztop_epc_MS_buf, *ka_ztop_epc_MS_buf;
   double *time_buf;

   short badshort= -32768;
   float badfloat=-1000.0, weight, total_weight, lat_ratio, lon_ratio, time_ratio;
   int ntotal_weight, ilev_lowest;

   size_t nlat, nlon, nlevel;

   float lat1, lon1, lat2, lon2;
   int icube, isecs, iday, ilat1, ilon1, ilat2, ilon2, itime1, itime2;
   int lat_id, lon_id, time_id, height_id;
   int p_id, h_id, t_id, qv_id, lev_id, ps_id, hs_id, ts_id, t2m_id, slp_id, tqv_id;

   nc_type t_type, ts_type, tqv_type;
   int t_ndims, ts_ndims, tqv_ndims;
   int t_dimids[NC_MAX_VAR_DIMS], ts_dimids[NC_MAX_VAR_DIMS], tqv_dimids[NC_MAX_VAR_DIMS];
   int t_natts, ts_natts, tqv_natts;

   float ts_MERRA2, t2m_MERRA2, tqv_MERRA2;
   short h273_MERRA2;

   int plot_files= 0, elev, nku20, nku25;
   float e0, e1, frac, frac0, frac1, h273;

   double now;
   char cdate_start[64], cdate_end[64];
   /*--------------------------------------------*/

   now= time(NULL);
   t0= now;
   tm= *gmtime(&t0);
   sprintf(cdate_start,"%4d/%02d/%02d %02d:%02d:%02d", 1900+tm.tm_year, 1+tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
   /*
   sprintf(dbdir,"%s", "/Volumes/15TB-RAID/projects/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1");
   */

   sprintf(dbdir,"%s", "/home/utsumi/mnt/wellshare/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1");


   sprintf(gprof_file,"%s", "none");
   sprintf(model_file,"%s", "none");

   if ( argc < 5 ) {
     printf("%s GMI_1B_FILE  target_lat target_lon  OUTDIR  <GPROF-FILE>  <MODEL-FILE>\n", argv[0]);
     exit(1);
   }

   tlat= atof(argv[2]);
   tlon= atof(argv[3]);
   if ( tlat < -90 )
     printf("Processing orbit over continental US region\n");
   else if ( tlat > 90 )
     printf("Processing entire orbit\n");
   else
     printf("Processing orbit near %f %f\n", tlat, tlon);
   sprintf(outdir,"%s", argv[4]);
   printf("outdir= %s\n", outdir);

   if ( argc >= 6 ) {
     sprintf(gprof_file,"%s", argv[5]);
     printf("gprof= %s\n", gprof_file);
   }
   if ( argc == 7 ) {
     sprintf(model_file,"%s", argv[6]);
     printf("model= %s\n", model_file);
   }

   /*sprintf(outdir,"%s", "/Users/jturk/code/GPM/GPM_DB_V5/EPC_DB_ENTRIES_V5_MRMS");*/

   bytes= sizeof(strm_record);
   printf("record length= %d\n", bytes);

   tmp1= (int *) malloc(DB_MAXREC*sizeof(int));
   tmp2= (int *) malloc(DB_MAXREC*sizeof(int));
   tmp3= (float *) malloc(DB_MAXREC*sizeof(float));
   tmp7= (float *) malloc(DB_MAXREC*sizeof(float));

   for (k= 0; k< NEM; k++) {
     epc_min[k]= 1.0E6;
     epc_max[k]= -1.0E6;
   }

   for (k= 0; k< NEM_USE; k++)
     pow2[k]= pow(NPCHIST,k);




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

   /*  Read PC coefficients file */
   /*--- modify ---*/
   /*sprintf(sdir,"%s", "/Users/jturk/code/radtran/CSAT-AMSU-MODELING/EMIS-MODELING/e9_Ts_Vap_NTB9_NREG57_NPC3_NKu20");*/
   /*
   sprintf(sdir,"%s", "/Users/jturk/code/GPM/GPM_DB_V5/e9_Ts_Vap_NTB9_NREG57_NPC3_NKu20");
   */
   sprintf(sdir,"%s", "/home/utsumi/bin/ENSPR");

   sprintf(svar,"%s/%s", sdir, "coef_pc.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i=0; i<= NREG; i++) {
     fscanf(fpc,"%d ", &k);
     for (j=0; j< NEM; j++) {
       fscanf(fpc,"%f ", &b_all[i][j]);
       printf("%f ", b_all[i][j]);
     }
     fscanf(fpc,"\n");
     printf("\n");
   }
   fclose(fpc);


   /*  Read eigenvalues wmatrix file */

   sprintf(svar,"%s/%s", sdir, "wmatrix.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NEM; i++) {
     for (j = 0; j < NEM; j++) fscanf(fpc,"%f ",&u[i+1][j+1]);
     fscanf(fpc,"\n");
     for (j = 0; j < NEM; j++) printf("%f ",u[i+1][j+1]);
     printf("\n");
   }
   fclose(fpc);

   /*  Read ave emissivity file */

   sprintf(svar,"%s/%s", sdir, "ave_emis.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NEM; i++) {
     fscanf(fpc,"%d %f %f\n", &j, &ave_emis[i], &std_emis[i]);
     printf("%f\n", ave_emis[i]);
   }
   fclose(fpc);

   /*  Read ave PC(e) file */

   sprintf(svar,"%s/%s", sdir, "ave_pc.txt");
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

   /*  Read ave TB file */
   /*sprintf(svar,"%s/../../%s", sdir, "ave_tb.txt");*/
   sprintf(svar,"%s/%s", sdir, "ave_tb.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NCHAN; i++) {
     fscanf(fpc,"%d %f %f\n", &j, &ave_tb[i], &std_tb[i]);
     printf("%2d %12.6f %12.6f\n", i, ave_tb[i], std_tb[i]);
   }
   fclose(fpc);

   /*---------------- Read GMI -----------------------------------*/

   isatname= 1;
   sprintf(satname,"%s","GPM");
   sprintf(sensor,"%s","GMI");

   sprintf(tbfile,"%s", argv[1]);

   status = nc_open(tbfile, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n",tbfile);
     exit(-1);
   }
   printf("Opened %s\n",tbfile);

   if ((status= nc_inq_ncid(ncid, "S1", &grpid1))) ERR(status);
   if ((status= nc_inq_ncid(grpid1, "ScanTime", &grpid1_scantime))) ERR(status);
   if ((status= nc_inq_ncid(grpid1, "navigation", &grpid1_nav))) ERR(status);
   if ((status= nc_inq_ncid(grpid1, "scanStatus", &grpid1_scanstatus))) ERR(status);

   if ((status= nc_inq_ncid(ncid, "S2", &grpid2))) ERR(status);
   if ((status= nc_inq_ncid(grpid2, "scanStatus", &grpid2_scanstatus))) ERR(status);

   /*--- S1 latitude and longitude ---*/

   status= nc_inq_varid(grpid1, "Latitude", &lat_S1_id);
   if ( status != 0 )  ERR(status);
   status= nc_inq_varid(grpid1, "Longitude", &lon_S1_id);
   if ( status != 0 )  ERR(status);

   status = nc_inq_var(grpid1, lat_S1_id, 0, &lat_S1_type, &lat_S1_ndims, lat_S1_dimids, &lat_S1_natts);
   if ( status != 0 )  ERR(status);
   printf("latitude ndims= %d\n", lat_S1_ndims);
   status = nc_inq_dimlen(grpid1, lat_S1_dimids[0], &nl_S1); 
   status = nc_inq_dimlen(grpid1, lat_S1_dimids[1], &ns_S1); 
   printf("S1 latitude dims= %d %d\n", (int) nl_S1, (int) ns_S1);

   lat_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
   lon_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
   status = nc_get_var_float(grpid1, lat_S1_id, lat_S1_buf);
   if ( status != 0 ) ERR(status);
   status = nc_get_var_float(grpid1, lon_S1_id, lon_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 TB ---*/

   status= nc_inq_varid(grpid1, "Tb", &tb_S1_id);
   if ( status != 0 )  ERR(status);

   status = nc_inq_var(grpid1, tb_S1_id, 0, &tb_S1_type, &tb_S1_ndims, tb_S1_dimids, &tb_S1_natts);
   if ( status != 0 )  ERR(status);
   printf("S1 TB ndims= %d\n", tb_S1_ndims);
   status = nc_inq_dimlen(grpid1, tb_S1_dimids[0], &nl_S1);                  
   status = nc_inq_dimlen(grpid1, tb_S1_dimids[1], &ns_S1);                  
   status = nc_inq_dimlen(grpid1, tb_S1_dimids[2], &nc_S1);                  
   printf("S1 TB dims= %d %d %d\n", (int) nl_S1, (int) ns_S1, (int) nc_S1);

   tb_S1_buf= (float *) malloc(nl_S1*ns_S1*nc_S1*sizeof(float));
   status = nc_get_var_float(grpid1, tb_S1_id, tb_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S2 latitude and longitude ---*/

   status= nc_inq_varid(grpid2, "Latitude", &lat_S2_id);
   if ( status != 0 )  ERR(status);
   status= nc_inq_varid(grpid2, "Longitude", &lon_S2_id);
   if ( status != 0 )  ERR(status);

   status = nc_inq_var(grpid2, lat_S2_id, 0, &lat_S2_type, &lat_S2_ndims, lat_S2_dimids, &lat_S2_natts);
   if ( status != 0 )  ERR(status);
   printf("latitude ndims= %d\n", lat_S2_ndims);
   status = nc_inq_dimlen(grpid2, lat_S2_dimids[0], &nl_S2);
   status = nc_inq_dimlen(grpid2, lat_S2_dimids[1], &ns_S2);
   printf("S2 latitude dims= %d %d\n", (int) nl_S2, (int) ns_S2);

   lat_S2_buf= (float *) malloc(nl_S2*ns_S2*sizeof(float));
   lon_S2_buf= (float *) malloc(nl_S2*ns_S2*sizeof(float));
   status = nc_get_var_float(grpid2, lat_S2_id, lat_S2_buf);
   if ( status != 0 ) ERR(status);
   status = nc_get_var_float(grpid2, lon_S2_id, lon_S2_buf);
   if ( status != 0 ) ERR(status);

   /*--- S2 TB ---*/

   status= nc_inq_varid(grpid2, "Tb", &tb_S2_id);
   if ( status != 0 )  ERR(status);

   status = nc_inq_var(grpid2, tb_S2_id, 0, &tb_S2_type, &tb_S2_ndims, tb_S2_dimids, &tb_S2_natts);
   if ( status != 0 )  ERR(status);
   printf("S2 TB ndims= %d\n", tb_S2_ndims);
   status = nc_inq_dimlen(grpid2, tb_S2_dimids[0], &nl_S2);
   status = nc_inq_dimlen(grpid2, tb_S2_dimids[1], &ns_S2);
   status = nc_inq_dimlen(grpid2, tb_S2_dimids[2], &nc_S2);
   printf("S2 TB dims= %d %d %d\n", (int) nl_S2, (int) ns_S2, (int) nc_S2);

   tb_S2_buf= (float *) malloc(nl_S2*ns_S2*nc_S2*sizeof(float));
   status = nc_get_var_float(grpid2, tb_S2_id, tb_S2_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 date ---*/

   status= nc_inq_varid(grpid1_scantime, "Year", &yyyy_S1_id);
   if ( status != 0 )  ERR(status);
   yyyy_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, yyyy_S1_id, yyyy_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Month", &mm_S1_id);
   if ( status != 0 )  ERR(status);
   mm_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mm_S1_id, mm_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "DayOfMonth", &dd_S1_id);
   if ( status != 0 )  ERR(status);
   dd_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, dd_S1_id, dd_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Hour", &hh_S1_id);
   if ( status != 0 )  ERR(status);
   hh_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, hh_S1_id, hh_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Minute", &mn_S1_id);
   if ( status != 0 )  ERR(status);
   mn_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mn_S1_id, mn_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Second", &ss_S1_id);
   if ( status != 0 )  ERR(status);
   ss_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, ss_S1_id, ss_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 spacecraft position ---*/

   status= nc_inq_varid(grpid1_nav, "scLat", &sclat_S1_id);
   if ( status != 0 )  ERR(status);
   sclat_S1_buf= (float *) malloc(nl_S1*sizeof(float));
   status = nc_get_var_float(grpid1_nav, sclat_S1_id, sclat_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_nav, "scLon", &sclon_S1_id);
   if ( status != 0 )  ERR(status);
   sclon_S1_buf= (float *) malloc(nl_S1*sizeof(float));
   status = nc_get_var_float(grpid1_nav, sclon_S1_id, sclon_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 granule number ---*/
   status= nc_inq_varid(grpid1_scanstatus, "FractionalGranuleNumber", &FractionalGranuleNumber_S1_id);
   if ( status != 0 )  ERR(status);
   FractionalGranuleNumber_S1_buf= (double *) malloc(nl_S1*sizeof(double));
   status = nc_get_var_double(grpid1_scanstatus, FractionalGranuleNumber_S1_id, FractionalGranuleNumber_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 data quality ---*/

   status= nc_inq_varid(grpid1_scanstatus, "dataQuality", &qual_S1_id);
   if ( status != 0 )  ERR(status);
   qual_S1_buf= (unsigned char *) malloc(nl_S1*sizeof(unsigned char));
   status = nc_get_var_uchar(grpid1_scanstatus, qual_S1_id, qual_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S2 data quality ---*/

   status= nc_inq_varid(grpid2_scanstatus, "dataQuality", &qual_S2_id);
   if ( status != 0 )  ERR(status);
   qual_S2_buf= (unsigned char *) malloc(nl_S2*sizeof(unsigned char));
   status = nc_get_var_uchar(grpid2_scanstatus, qual_S2_id, qual_S2_buf);
   if ( status != 0 ) ERR(status);

   if ((retval = nc_close(ncid))) ERR(retval);
   printf("Closed %s\n", tbfile);


   tb_S1S2_buf= (float *) malloc(nl_S1*ns_S1*(nc_S1+nc_S2)*sizeof(float));
   secs_S1_buf= (double *) malloc(nl_S1*sizeof(double));

   for (i=0; i< nl_S1*ns_S1*(nc_S1+nc_S2); i++) {
     tb_S1S2_buf[i]= -999.0;
   }

   for (i=0; i < nl_S1; i++) {
     secs_S1_buf[i]= -1.0;
     irev= FractionalGranuleNumber_S1_buf[i];
     iyyyy= yyyy_S1_buf[i];
     imm= mm_S1_buf[i];
     idd= dd_S1_buf[i];
     ihh= hh_S1_buf[i];
     imn= mn_S1_buf[i];
     iss= ss_S1_buf[i];
     tm.tm_year= iyyyy-1900;
     tm.tm_mon= imm-1;
     tm.tm_mday= idd;
     tm.tm_hour= ihh;
     tm.tm_min= imn;
     tm.tm_sec= iss;
     t0= timegm(&tm);
     tm= *gmtime(&t0);
     sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
     /*printf("GMI scan=%4d rev=%6d date=%s\n", i, irev, tdate);*/
     if ( i == 0 ) {
       strcpy(tdate1, tdate);
     }
     strcpy(tdate2, tdate);
     secs_S1_buf[i]= ( double ) t0;
   }

   printf("GMI rev=%6d Start=%s  End=%s\n", irev, tdate1, tdate2);

   /*---------------- GPROF -----------------------------------*/

   status = nc_open(gprof_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", gprof_file);
   }
   else {
     have_gprof= 1;
     printf("Opened %s\n", gprof_file);
     if ((status= nc_inq_ncid(ncid, "S1", &grpid1)))
     ERR(status);

     /*--- GPROF latitude and longitude ---*/

     status= nc_inq_varid(grpid1, "Latitude", &lat_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "Longitude", &lon_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "surfaceTypeIndex", &sfc_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "qualityFlag", &qualityFlag_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "surfacePrecipitation", &surfacePrecipitation_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "probabilityOfPrecip", &probabilityOfPrecip_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "temp2mIndex", &temp2mIndex_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "totalColumnWaterVaporIndex", &totalColumnWaterVaporIndex_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "frozenPrecipitation", &frozenPrecipitation_SG_id);
     if ( status != 0 )  ERR(status);

     status = nc_inq_var(grpid1, lat_SG_id, 0, &lat_SG_type, &lat_SG_ndims, lat_SG_dimids, &lat_SG_natts);
     if ( status != 0 )  ERR(status);
     printf("latitude ndims= %d\n", lat_SG_ndims);
     status = nc_inq_dimlen(grpid1, lat_SG_dimids[0], &nl_SG);
     status = nc_inq_dimlen(grpid1, lat_SG_dimids[1], &ns_SG);
     printf("GPROF latitude dims= %d %d\n", (int) nl_SG, (int) ns_SG);
     if ( nl_SG != nl_S1 && ns_SG != ns_S1 ) {
       printf("Dimensions are not parallel  GPROF=%d %d  GMI=%d %d\n",
         (int) nl_SG, (int) ns_SG, (int) nl_S1, (int) ns_S1);
       exit(1);
     }
     lat_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     lon_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     sfc_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     qualityFlag_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     surfacePrecipitation_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     probabilityOfPrecip_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     totalColumnWaterVaporIndex_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     temp2mIndex_GPROF_buf= (short *) malloc(nl_SG*ns_SG*sizeof(short));
     frozenPrecipitation_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));

     status = nc_get_var_float(grpid1, lat_SG_id, lat_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, lon_SG_id, lon_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid1, sfc_SG_id, sfc_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid1, qualityFlag_SG_id, qualityFlag_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, surfacePrecipitation_SG_id, surfacePrecipitation_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, probabilityOfPrecip_SG_id, probabilityOfPrecip_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_short(grpid1, temp2mIndex_SG_id, temp2mIndex_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid1, totalColumnWaterVaporIndex_SG_id, totalColumnWaterVaporIndex_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, frozenPrecipitation_SG_id, frozenPrecipitation_GPROF_buf);
     if ( status != 0 ) ERR(status);

     if ((retval = nc_close(ncid))) ERR(retval);
     printf("Closed %s\n", gprof_file);

     /*for (i=0; i < nl_S1; i++) {
       for (j=0; j < ns_S1; j++) {
         idx= i*ns_S1 + j;
         glat1= lat_GPROF_buf[idx];
         glon1= lon_GPROF_buf[idx];
         sfc_class= sfc_GPROF_buf[idx];
         printf("%4d %4d  %8.3f %8.3f %3d %3d ",i,j,glat1,glon1,qualityFlag_GPROF_buf[idx],sfc_class);
         printf("%9.3f ", surfacePrecipitation_GPROF_buf[idx]);
         printf("%9.3f ", probabilityOfPrecip_GPROF_buf[idx]);
         printf("%9.3f ", frozenPrecipitation_GPROF_buf[idx]);
         printf("%6d ", totalColumnWaterVaporIndex_GPROF_buf[idx]);
         printf("%6d ", temp2mIndex_GPROF_buf[idx]);
         printf("\n");
       }
     }*/
   }

   /*--------------------------------------------*/

   irev= FractionalGranuleNumber_S1_buf[nl_S1/2];

   if ( tlat > -90 && tlat < 90 ) {
     km_min= 1.0E8;
     i= -1;
     for (isat=0; isat < nl_S1; isat++) {
       idx= isat*ns_S1 + (ns_S1/2);
       glat1= lat_S1_buf[idx];
       glon1= lon_S1_buf[idx];
       km= RADEARTH*acos(cos(DTR*glon1-DTR*tlon)*cos(DTR*glat1)*cos(DTR*tlat) + sin(DTR*glat1)*sin(DTR*tlat));
       if ( km < km_min ) {
         km_min= km;
         i= isat;
         clat= glat1;
         clon= glon1;
       }
     }
     if ( i < 0 || km_min > 500 ) {
       printf("Starting line not found  km_min=%f\n", km_min);
       exit(1);
     }
     printf("Starting line %d   km_min=%f\n", i, km_min);
     isat_start= i-30;
     if ( isat_start < 0 ) isat_start= 0;
     isat_end= i+30;
     if ( isat_end >= nl_S1 ) isat_end= nl_S1-1;

   }
   else if ( tlat < -90 ) {
     /* Check for MRMS coverage */
     mrms= 1;
     for (isat=0; isat < nl_S1; isat++) {
       idx= isat*ns_S1 + (ns_S1/2);
       glat1= lat_S1_buf[idx];
       glon1= lon_S1_buf[idx];
       if ( glat1 > 20 && glat1 < 55 && glon1 > -130 && glon1 < -60 ) {
         if (isat_start < 0 ) isat_start= isat;
         isat_end= isat;
         n_mrms++;
       }
     }
     if ( n_mrms < 20 ) {
       mrms= 0;
       printf("Percent= 0 over MRMS domain\n");
       exit(1);
     }
   }
   else {
     isat_start= 0;
     isat_end= nl_S1-1;
   }

   if ( isat_start < 0 || isat_end < 0 ) {
     printf("Starting line not found\n");
     exit(1);
   }

   nl_out= isat_end - isat_start + 1;
   printf("Lines=%d  Range=%d %d\n", nl_out, isat_start, isat_end);

   i= (isat_end + isat_start)/2;
   idx= i*ns_S1 + (ns_S1/2);
   clat= lat_S1_buf[idx];
   clon= lon_S1_buf[idx];
   iyyyy= yyyy_S1_buf[i];
   imm= mm_S1_buf[i];
   idd= dd_S1_buf[i];
   ihh= hh_S1_buf[i];
   imn= mn_S1_buf[i];
   iss= ss_S1_buf[i];
   sprintf(cdate,"%04d/%02d/%02d %02d:%02d:%02d", iyyyy, imm, idd, ihh, imn, iss);
   sprintf(outdate,"%04d%02d%02d_%02d%02d", iyyyy, imm, idd, ihh, imn);
   sprintf(outfile, "%s/%s_EPC_%06d_%s.NS_MS.nc", outdir, satname, irev, outdate);

   printf("Center coords=%8.3f %8.3f date=%s\n", clat, clon, cdate);

   /*--------------------------------------------*/

   /* emis and pc_emis computed from GMI input TB */
   emis_S1S2_buf= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));
   pc_emis_S1S2_buf= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));

   /* estimated precipitation   NS and MS database */
   precip_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip_prob_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_prob_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   qual_NS_buf= (int *) malloc(nl_out*ns_S1*sizeof(int));

   precip_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip_prob_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_prob_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   qual_MS_buf= (int *) malloc(nl_out*ns_S1*sizeof(int));

   /* additional estimated model parameters */
   ts_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   t2m_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   tqv_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   h273_epc_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ku_ztop_epc_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   ts_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   t2m_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   tqv_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   h273_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ku_ztop_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ka_ztop_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   /* top-ranked Z and TB candidates */
   z_ku_1_NS_buf= (short *) malloc(nl_out*ns_S1*NLEV_NS*sizeof(short));
   z_ku_1_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_NS*sizeof(short));
   z_ka_1_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_NS*sizeof(short));
   tb_1_NS_buf= (float *) malloc(nl_out*ns_S1*(nc_S1+nc_S2)*sizeof(float));
   tb_1_MS_buf= (float *) malloc(nl_out*ns_S1*(nc_S1+nc_S2)*sizeof(float));

   for (i=0; i< nl_out*ns_S1*NEM; i++) {
     emis_S1S2_buf[i]= -999.0;
     pc_emis_S1S2_buf[i]= -999.0;
   }
   for (i=0; i< nl_out*ns_S1; i++) {
     precip_NS_buf[i]= -999.0;
     precip2_NS_buf[i]= -999.0;
     precip_prob_NS_buf[i]= -999.0;
     precip2_prob_NS_buf[i]= -999.0;
     ts_epc_NS_buf[i]= -999.0;
     t2m_epc_NS_buf[i]= -999.0;
     h273_epc_NS_buf[i]= -9999;
     tqv_epc_NS_buf[i]= -999.0;
     qual_NS_buf[i]= -999;
     ku_ztop_epc_NS_buf[i]= -9999;

     precip_MS_buf[i]= -999.0;
     precip2_MS_buf[i]= -999.0;
     precip_prob_MS_buf[i]= -999.0;
     precip2_prob_MS_buf[i]= -999.0;
     ts_epc_MS_buf[i]= -999.0;
     t2m_epc_MS_buf[i]= -999.0;
     h273_epc_MS_buf[i]= -9999;
     tqv_epc_MS_buf[i]= -999.0;
     qual_MS_buf[i]= -999;
     ku_ztop_epc_MS_buf[i]= -9999;
     ka_ztop_epc_MS_buf[i]= -9999;
   }
   for (i=0; i < nl_out*ns_S1*NLEV_NS; i++) {
     z_ku_1_NS_buf[i]= -9999;
     z_ku_1_MS_buf[i]= z_ka_1_MS_buf[i]= -9999;
   }
   for (i=0; i< nl_out*ns_S1*(nc_S1+nc_S2); i++) {
     tb_1_NS_buf[i]= -999.0;
     tb_1_MS_buf[i]= -999.0;
   }


   pixel_index = (int *) malloc(nl_S1*ns_S1*sizeof(int));
   for (i=0; i< nl_S1*ns_S1; i++) pixel_index[i]= -1;
   n_index = (int *) malloc(ndb_files*sizeof(int));
   for (i=0; i< ndb_files; i++) n_index[i]= 0;

   ncold_index = (int *) malloc(ndb_files*sizeof(int));
   for (i=0; i< ndb_files; i++) ncold_index[i]= 0;

   
  /*----------------------- optional MRMS check ----------------------------*/

     if ( mrms > 0 ) {
       tm.tm_year= iyyyy-1900;
       tm.tm_mon= imm-1;
       tm.tm_mday= idd;
       tm.tm_hour= ihh;
       tm.tm_min= imn;
       tm.tm_sec= 0;
       t0= timegm(&tm);

       found= -1;
       for (m= -20; m<= 20; m++) {
         t1= t0 + 60*m;
         tm= *gmtime(&t1);
         sprintf(mdate,"%4d%02d%02d.%02d%02d%02d", 1900+tm.tm_year, 1+tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
         sprintf(mrms_name,"/Volumes/12TB-RAID-3/satdata/NMQ/GMI_EPC/GMI-%d-V04A.MRMS-%s.matched-130W_55W_20N_55N.extract.dat", irev, mdate);
         printf("%s %s %s\n", cdate, mdate, mrms_name);

         if ((fmrms= fopen(mrms_name,"r")) == 0) {
           printf("Percent   Not found %s\n", mrms_name);
           continue;
         }
         printf("Percent  Opened %s\n", mrms_name);
         fclose(fmrms);
         found= 1;
         break;
       }
       if ( found < 0 ) {
         printf("Percent  MRMS files not located near %s\n", cdate);
         exit(1);
       }
     }


   /*---------------------------------------------------------------*/

   /* Match S2 to S1 */

   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat=0; jsat <= ns_S1; jsat++) {
       idx1= isat*ns_S1 + jsat;
       glat1= lat_S1_buf[idx1];
       glon1= lon_S1_buf[idx1];

       /*printf ("%4d %4d %8.3f %8.3f  ", isat, jsat, glat1, glon1);*/
       for (k=0; k < nc_S1; k++) {
         /* 10V 10H 18V 18H 21V 36V 36H 89V 89H */
         idx2= isat*ns_S1*nc_S1 + jsat*nc_S1 + k;
         idx3= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;
         tb_S1S2_buf[idx3]= tb_S1_buf[idx2];
         /*printf ("%6.2f ", tb_S1_buf[idx2]);*/
       }

       km_min= 1000.0;
       isat2_min= jsat2_min= -1;
       /* Matching S2 line is 4-12 scans before S1 scan, and S2 sample is within +/-28 samples of S1 sample */
       /* but the spacecraft can re-orient so make it +/-16 lines to be sure */

       /* works but only fore direction     for (isat2=isat; isat2 >= isat-16; isat2--) { */
       for (itest=0; itest<= 32; itest++) {
         if ( km_min < 6.0 ) break;
         isat2= ( itest % 2 == 0 ) ? isat-(itest/2) : isat+(itest/2)+1; 
         qual= qual_S2_buf[isat2];
         if ( qual != 0 ) continue;
         if ( isat2 < 0 || isat2 >= nl_S2 ) continue;
         for (jsat2=jsat-30; jsat2 <= jsat+30; jsat2++) {
           if ( jsat2 < 0 || jsat2 >= ns_S2 ) continue;
           idx2= isat2*ns_S2 + jsat2;
           glat2= lat_S2_buf[idx2];
           if ( fabs(glat1-glat2) > 1.0 ) continue;
           glon2= lon_S2_buf[idx2];
           km= RADEARTH*acos(cos(DTR*glon1-DTR*glon2)*cos(DTR*glat1)*cos(DTR*glat2) + sin(DTR*glat1)*sin(DTR*glat2));
           if ( km < km_min ) {
             km_min= km;
             isat2_min= isat2; jsat2_min= jsat2;
             if ( km_min < 6.0 ) break;
           }
         }
       }
       if ( km_min < 7.0 && isat2_min >= 0 && jsat2_min >= 0 ) { 
         idx2= isat2_min*ns_S2 + jsat2_min;
         glat2= lat_S2_buf[idx2];
         glon2= lon_S2_buf[idx2];
         /*printf ("%6.2f %4d %4d %8.3f %8.3f  ", km_min, isat2_min, jsat2_min, glat2, glon2);*/
 
         for (k=0; k < nc_S2; k++) {
           /* 166V 166H 183/3 183/8 */
           idx2= isat2_min*ns_S2*nc_S2 + jsat2_min*nc_S2 + k;
           idx3= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + (nc_S1+k);
           tb_S1S2_buf[idx3]= tb_S2_buf[idx2];
           /*printf ("%6.2f ", tb_S2_buf[idx2]);*/
         }
         /*printf("\n");*/
         ishift= isat-isat2_min;
         if ( ishift < ishift_min ) ishift_min= ishift;
         if ( ishift > ishift_max ) ishift_max= ishift;
         jshift= jsat-jsat2_min;
         if ( jshift < jshift_min ) jshift_min= jshift;
         if ( jshift > jshift_max ) jshift_max= jshift;
       }
       else { 
         /*printf("\n");*/
       }
     }
   }

   /*---------------------------------------------------------------*/

   /* MERRA2 interpolated model fields */

   status = nc_open(model_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("Model file not found %s\n", model_file);
   }
   else {
     have_model= 1;
     printf("Opened %s\n", model_file);

     if ((status= nc_inq_varid(ncid, "t", &t_id))) ERR(status);
     if ((status = nc_inq_var(ncid, t_id, 0, &t_type, &t_ndims, t_dimids, &t_natts))) ERR(status);
     status = nc_inq_dimlen(ncid, t_dimids[0], &nlat);
     status = nc_inq_dimlen(ncid, t_dimids[1], &nlon);
     status = nc_inq_dimlen(ncid, t_dimids[2], &nlevel);
     t_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, t_id, t_buf))) ERR(status);
     printf("T dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     /* number of levels of the interpolated model fields is the same as the number of levels in the database entries */
     /* For MERRA2, there are 42 levels which are indexed from top to surface */
     if ( nlevel != 42 ) {
       printf("MERRA2 levels=%d are different than the 42-level database entries\n", (int) nlevel);
       exit(1);
     }
     /* interpolated model fields are on the identical grid as the GMI 1B file */
     if ( nlat != nl_S1 || nlon != ns_S1 ) {
       printf("Dimensions are not parallel  MERRA2=%d %d  GMI=%d %d\n", (int) nlat, (int) nlon, (int) nl_S1, (int) ns_S1);
       exit(1);
     }

     if ((status= nc_inq_varid(ncid, "time", &time_id))) ERR(status);
     time_buf= (double *) malloc(nlat*sizeof(double));
     if ((status = nc_get_var_double(ncid, time_id, time_buf))) ERR(status);
     printf("TIME dims= %d\n", (int) nlat);

     if ((status= nc_inq_varid(ncid, "h", &h_id))) ERR(status);
     h_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, h_id, h_buf))) ERR(status);
     printf("H dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     if ((status= nc_inq_varid(ncid, "qv", &qv_id))) ERR(status);
     qv_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, qv_id, qv_buf))) ERR(status);
     printf("QV dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     if ((status= nc_inq_varid(ncid, "ts", &ts_id))) ERR(status);
     ts_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, ts_id, ts_buf))) ERR(status);
     printf("TS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "t2m", &t2m_id))) ERR(status);
     t2m_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, t2m_id, t2m_buf))) ERR(status);
     printf("T2M dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "tqv", &tqv_id))) ERR(status);
     tqv_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, tqv_id, tqv_buf))) ERR(status);
     printf("TQV dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "hs", &hs_id))) ERR(status);
     hs_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, hs_id, hs_buf))) ERR(status);
     printf("HS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "ps", &ps_id))) ERR(status);
     ps_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, ps_id, ps_buf))) ERR(status);
     printf("PS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "latitude", &lat_id))) ERR(status);
     lat_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, lat_id, lat_buf))) ERR(status);
     printf("LAT dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "longitude", &lon_id))) ERR(status);
     lon_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, lon_id, lon_buf))) ERR(status);
     printf("LON dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "levels", &p_id))) ERR(status);
     p_buf= (float *) malloc(nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, p_id, p_buf))) ERR(status);
     printf("P dims= %d\n", (int) nlevel);

     status = nc_close(ncid);
     if ( status != 0 ) ERR(status);

     nbad= 0;
     for (i= 0; i< nlat; i++) {
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         lat1= lat_buf[idx];
         lat2= lat_S1_buf[idx];
         lon1= lon_buf[idx];
         lon2= lon_S1_buf[idx];
         if ( ( fabs(lat1-lat2) > 0.2 ) || ( fabs(lon) < 178 && fabs(lon1-lon2) > 0.2 )) {
           printf("GMI coords=%8.3f %8.3f  model=%8.3f %8.3f\n", lat_S1_buf[idx], lon_S1_buf[idx], lat_buf[idx], lon_buf[idx]);
           nbad++;
         }
       }
     }
     if ( nbad > 0 ) {
       printf("GMI and model coordinates are not aligned  N=%d\n", nbad);
       exit(1);
     }

     /* Freezing level height */

     h273_buf= (short *) malloc(nlat*nlon*sizeof(short));
     for (i= 0; i< nlat; i++) {
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         h273_buf[idx]= -999.0;
         h273= -1.0;
         for (j1= 0; j1 < nlevel; j1++) {
           idx3= i*ns_S1*nlevel + j*nlevel + j1;
           /*printf("%2d ModelH= %9.4f  ModelP=%9.4f  ModelT=%9.4f\n", j1, h_buf[idx3], p_buf[j1], t_buf[idx3]);*/
           if ( j1 == 0 ) continue;
           idx2= i*ns_S1*nlevel + j*nlevel + (j1-1);
           if ( h_buf[idx3] < -2000 || h_buf[idx3] > 10000 ) continue;
           if ( h_buf[idx2] < -2000 || h_buf[idx2] > 10000 ) continue;
           if ( t_buf[idx2] < 273 && t_buf[idx3] >= 273 ) {
             frac= ( 273.0 - t_buf[idx2])/(t_buf[idx3] - t_buf[idx2]);
             h273= h_buf[idx2] + frac*(h_buf[idx3] - h_buf[idx2]);
             /*printf("FREEZLEV  %7.2f %7.2f ModelH=%9.4f %9.4f  ModelT=%9.4f %9.4f Frac=%9.4f H273=%9.4f\n",
               lat_buf[idx], lon_buf[idx], h_buf[idx3], h_buf[idx2], t_buf[idx3], t_buf[idx2], frac, h273);*/
             h273_buf[idx]= h273;
             break;
           }
         }
       }
     }

   }  

/*
   if ( have_model > 0 ) {
     for (i= 0; i< nlat; i++) {
       t1= time_buf[i];
       tm= *gmtime(&t1);
       sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         printf("%4d %4d %8.3f %8.3f %s  %8.3f %8.3f %8.3f  %8.3f %8.3f %6d\n", i, j, lat_buf[idx], lon_buf[idx], tdate, ts_buf[idx], t2m_buf[idx], tqv_buf[idx], hs_buf[idx], ps_buf[idx], h273_buf[idx]);
         for (k= 0; k< nlevel; k++) {
           idx1= i*nlon*nlevel + j*nlevel + k;
           printf(" %4d %8.2f %8.2f %8.2f %12.8f\n", k, p_buf[k], h_buf[idx1], t_buf[idx1], qv_buf[idx1]);
         }
       }
     } 
   }
*/

   /*---------------------------------------------------------------*/

   /* Assign DB file index for each pixel */

   n1= n2= 0;
   printf("Locating number of DB files required \n");

   nsat= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {

     iyyyy= yyyy_S1_buf[isat];
     if ( iyyyy < 2014 ) {
       for (jsat=0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         qual_NS_buf[idx]= -3;
         qual_MS_buf[idx]= -3;
       }
       continue;  /* before GPM launch */
     }
     imm= mm_S1_buf[isat];
     idd= dd_S1_buf[isat];
     ihh= hh_S1_buf[isat];
     imn= mn_S1_buf[isat];
     iss= ss_S1_buf[isat];
     tm.tm_year= iyyyy-1900;
     tm.tm_mon= imm-1;
     tm.tm_mday= idd;
     tm.tm_hour= ihh;
     tm.tm_min= imn;
     tm.tm_sec= iss;
     t1= timegm(&tm);
     if ( t1 < 100 ) continue; /* avoid zero */
     if ( t1 > t0 ) t0= t1;  /* oldest time in the files */
     tm= *gmtime(&t1);
     sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

     isecs= 3600*ihh + 60*imn + iss;
     if ( isecs < 0 || isecs >= 86400 ) {
       printf("ERROR isecs=%d\n", isecs);
       exit(1);
     }

     for (jsat=0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;

       glat1= lat_S1_buf[idx];
       glon1= lon_S1_buf[idx];
       if ( fabs(glat1) > 90 || fabs(glon1) > 180 ) {
         qual_NS_buf[idx]= -4;
         qual_MS_buf[idx]= -4;
         continue;
       }
       slat= sclat_S1_buf[isat];
       slon= sclon_S1_buf[isat];
       if ( fabs(slat) > 90 || fabs(slon) > 180 ) {
         qual_NS_buf[idx]= -5;
         qual_MS_buf[idx]= -5;
         continue;
       }

       nbad= 0;
       for (k= 0; k< nc_S1+nc_S2; k++) {
         idx1= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;
         tb_gmi[k]= tb_S1S2_buf[idx1];
         if ( k < nc_S1 && ( tb_gmi[k] < 20 || tb_gmi[k] > 350 )) nbad++;
       }
       if ( nbad > 0 ) {
         qual_NS_buf[idx]= -1;
         qual_MS_buf[idx]= -1;
         continue;
       }

       /*if ( glat1 < -40 || glat1 > -38 || glon1 < -66 || glon1 > -63 ) continue;
       if ( tb_gmi[6] > 200 && tb_gmi[8] > 180 ) continue;*/

       /* If available, use interpolated MERRA2 fields */

       if ( have_model > 0 ) {
         ts_MERRA2= ts_buf[idx];
         t2m_MERRA2= t2m_buf[idx];
         tqv_MERRA2= tqv_buf[idx];
         h273_MERRA2= h273_buf[idx];
       }


       /*----------------------------------------------------*/

       tmp5[0]= 1.0;
       kt= 0;
       for (k= 0; k< NTBREG; k++) {
         tmp5[kt+1]= tb_gmi[k];
         kt++;
         for (k1= k; k1< NTBREG; k1++) {
           tmp5[kt+1]= tb_gmi[k]*tb_gmi[k1];
           kt++;
         }
       }
       tmp5[kt+1]= (tb_gmi[0]-tb_gmi[1])/(tb_gmi[0]+tb_gmi[1]);
       kt++;
       tmp5[kt+1]= (tb_gmi[2]-tb_gmi[3])/(tb_gmi[2]+tb_gmi[3]);
       kt++;
       tmp5[kt+1]= (tb_gmi[5]-tb_gmi[6])/(tb_gmi[5]+tb_gmi[6]);
       kt++;

       /* emissivity PC */
       for (k= 0; k< NEM; k++) {
         psum= 0.0;
         for (k2= 0; k2< NREG+1; k2++) psum+= b_all[k2][k]*tmp5[k2];
         pc_emis[k]= psum;
         /*printf("%10.6f ", pc_emis[k]);*/
       }
       /*printf("\n");*/

       /* emissivites from PC */
       for (k= 0; k< NEM; k++) {
         psum= 0.0;
         for (k2= 0; k2< NEM; k2++) psum+= u[k+1][k2+1] * pc_emis[k2];
         emis[k]= psum + ave_emis[k];
         /*printf("%10.6f ", emis[k]);*/
       }
       /*printf("\n");*/

       /* Save for output file */
       for (k=0; k < NEM; k++) {
         /*idx2= isat*ns_S1*NEM + jsat*NEM + k;*/
         idx_out= nsat*ns_S1*NEM + jsat*NEM + k;
         emis_S1S2_buf[idx_out]= emis[k];
         pc_emis_S1S2_buf[idx_out]= pc_emis[k];
       }

       /* first index moving fastest */

       for (i1= 0; i1< ndb_files; i1++) {
         k= 0;
         idxn[k]= i1 % NPCHIST;

         found= 0;
         for (j= 0; j < NPCHIST; j++) {
           pc_lo= pc_range[k][j][0];
           pc_hi= pc_range[k][j][1];
           /*printf("%2d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_emis[k], pc_lo, pc_hi, j , idxn[k] );*/
           if ( pc_emis[k] > pc_lo && pc_emis[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
         }
         if ( found != 1 ) continue;

         for (k= 1; k< NEM_USE-1; k++) {
           p= pow2[k];
           idxn[k]= (i1/p) % NPCHIST;

           found= 0;
           for (j= 0; j < NPCHIST; j++) {
             pc_lo= pc_range[k][j][0];
             pc_hi= pc_range[k][j][1];
             /*printf("%d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_emis[k], pc_lo, pc_hi, j , idxn[k] );*/
             if ( pc_emis[k] > pc_lo && pc_emis[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
           }
           if ( found != 1 ) goto bypass;
         }

         k= NEM_USE-1;
         idxn[k]= i1/pow2[k];

         found= 0;
         for (j= 0; j < NPCHIST; j++) {
           pc_lo= pc_range[k][j][0];
           pc_hi= pc_range[k][j][1];
           /*printf("%d  i1=%d  pc=%f  range=%f %f   j=%d idx=%d\n", k, i1, pc_emis[k], pc_lo, pc_hi, j, idxn[k] );*/
           if ( pc_emis[k] > pc_lo && pc_emis[k] <= pc_hi && j == idxn[k] ) {found=1; break;}
         }
         if ( found != 1 ) continue;

         idx_db= 0;
         for (k= 0; k< NEM_USE; k++)
           idx_db+= idxn[k]*pow2[k];

         /*for (k= 0; k< NEM_USE; k++) {
           j= idxn[k];
           printf("%1d ", j);
           printf("%12.6f %12.6f  %12.6f\n", pc_range[k][j][0], pc_range[k][j][1], pc_emis[k]);
         }*/

/*
           for (k= 0; k< NEM_USE; k++) printf("%d ", idxn[k]);
           printf("  index=%d\n", idx_db);
*/

         if ( idx_db < 0 || idx_db >= ndb_files ) {
           printf("Exceeded DB index=%d\n", idx_db);
           exit(1);
         }

         pixel_index[idx]= idx_db;
         if ( n_index[idx_db] == 0 ) n1++;
         n_index[idx_db]++;
         if ( have_model > 0 ) {
           if ( t2m_MERRA2 < t2m_MERRA2_min ) t2m_MERRA2_min= t2m_MERRA2;
           if ( t2m_MERRA2 > t2m_MERRA2_max ) t2m_MERRA2_max= t2m_MERRA2;
           if ( tqv_MERRA2 < tqv_MERRA2_min ) tqv_MERRA2_min= tqv_MERRA2;
           if ( tqv_MERRA2 > tqv_MERRA2_max ) tqv_MERRA2_max= tqv_MERRA2;
           if ( t2m_MERRA2 < 278 ) {
             if ( ncold_index[idx_db] == 0 ) n2++;
             ncold_index[idx_db]++;
           }
         }
         for (k= 0; k< NEM; k++) {
           if ( pc_emis[k] < epc_min[k] ) epc_min[k]= pc_emis[k];
           if ( pc_emis[k] > epc_max[k] ) epc_max[k]= pc_emis[k];
         }

bypass:;
       }


     }  /* next GMI sample */

     nsat++;

   }  /* next GMI line */


   printf("Number of DB files required %d\n", n1);
   if ( have_model > 0 ) {
     ncold= n2;
     printf("Number of DB files with cold T2m pixels %d\n", ncold);
     printf("T2m range= %8.3f %8.3f\n", t2m_MERRA2_min, t2m_MERRA2_max);
     printf("TQV range= %8.3f %8.3f\n", tqv_MERRA2_min, tqv_MERRA2_max);
   }
   for (k= 0; k< NEM; k++) {
     printf("EPC %2d range= %12.6f %12.6f\n", k, epc_min[k], epc_max[k]);
   }

   if ( n1 == 0 ) exit(1);

   /*---------------------------------------------------*/
   /*--- Process all data for each DB file only once ---*/

   for (idx_db= 0; idx_db< ndb_files; idx_db++) {

     /*if ( idx_db < 9272 ) continue;  TEST */

     if ( n_index[idx_db] == 0 ) continue;
     ndb= 0;
     
     
     /*--- Read the number of rain events for cold and warm T2m conditions ---*/
     /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
     /* second six= same for when T2m > 278K */

     sprintf(svar,"%s/db_%05d.bin.nrain.txt", dbdir, idx_db);
     if ((fdb_nrain = fopen(svar,"r")) == 0) {
       printf("Unable to locate database precip CDF file %s\n", svar);
       maxrec= DB_MAXREC;
     }
     else {
       printf("%s  COLD= ", dbfile);
       for (i= 0; i< 6; i++) {
         fscanf(fdb_nrain, "%d ", &j);
         nrain_db[i][0]= j;
         printf("%8d ", j);
       }
       printf("  WARM= ");
       for (i= 0; i< 6; i++) {
         fscanf(fdb_nrain, "%d ", &j);
         nrain_db[i][1]= j;
         printf("%8d ", j);
       }
       fclose(fdb_nrain);
   
       /*--- If the DB file contains very little rain for cold and warm T2m, don't process it */

       frac0= frac1= 0.0;   /* 0=cold T2m  1=warm T2m */
       if ( nrain_db[0][0] > 0 ) frac0= 1.0*nrain_db[1][0]/nrain_db[0][0];
       if ( nrain_db[0][1] > 0 ) frac1= 1.0*nrain_db[1][1]/nrain_db[0][1];

       printf("FracRain= %9.5f %9.5f ", frac0, frac1);
       printf("\n");

       maxrec= DB_MAXREC;
       if ( nrain_db[4][0] == 0 && nrain_db[4][1] == 0 ) 
         maxrec= 0.1*DB_MAXREC;
       else if ( nrain_db[3][0] == 0 && nrain_db[3][1] == 0 ) 
         maxrec= 0.01*DB_MAXREC;
       else if ( nrain_db[2][0] == 0 && nrain_db[2][1] == 0 ) {
         if ( frac0 < DB_RAINFRAC && frac1 < DB_RAINFRAC ) {
           printf("Insufficient fraction of DB entries exceeding 1 mm/hr= %9.5f %9.5f\n", frac0, frac1);
           npix+= n_index[idx_db];
           continue;
         }
         maxrec= 0.001*DB_MAXREC;
       }
     }

     /*--- Open DB file ---*/

     sprintf(dbfile,"%s/db_%05d.bin", dbdir, idx_db);
     /*printf("%s  MAXREC=%d\n", dbfile, maxrec);*/

     if ((fdb2 = fopen(dbfile,"r")) == 0) {
       printf("ERROR  Unable to locate %s\n", dbfile);

       /* Find out where this happened */
       nsat= 0;
       for (isat= isat_start; isat <= isat_end; isat++) {
         for (jsat=0; jsat <= ns_S1; jsat++) {
           idx1= isat*ns_S1 + jsat;
           idx_out= nsat*ns_S1 + jsat;
           if ( pixel_index[idx1] == idx_db ) {
             printf("Occured for %d %d %8.3f %8.3f  TB=", isat, jsat, lat_buf[idx1], lon_buf[idx1]);
             for (k=0; k< nc_S1+nc_S2; k++) {
               idx3= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;
               printf("%6.2f ", tb_S1S2_buf[idx3]);
             }
             printf("\n");
             qual_NS_buf[idx_out]= -2;
             qual_MS_buf[idx_out]= -2;
           }
         }
         nsat++;
       }
       continue;
     }


     /*-----------------------------------------------------------------------*/
     /* Read the index database file one time */

     nrec= nrec1= nrec2= nrain= 0;

     while ( fread(&prof,bytes,1,fdb2) == 1 ) {
       nrec++;

       /*nbad= 0;
       for (k= 0; k< NEM; k++) {
         if ( em_compare[k] != 1 ) continue;
         x= 0.1*(epc_max[k] - epc_min[k]);
         if ( prof.pc_emis[k] < epc_min[k] - x ) nbad++;
         if ( prof.pc_emis[k] > epc_max[k] + x ) nbad++;
       }
       if ( nbad > 0 ) continue;*/
       if ( have_model > 0 ) {
         if ( prof.t2m < t2m_MERRA2_min - MAX_T2M_DIFF ) continue;
         if ( prof.t2m > t2m_MERRA2_max + MAX_T2M_DIFF ) continue;
         /*if ( prof.tqv < tqv_MERRA2_min - MAX_TQV_DIFF ) continue;
         if ( prof.tqv > tqv_MERRA2_max + MAX_TQV_DIFF ) continue;*/
       }

       if ( prof.j_NS < 0 || prof.j_NS >= 49 ) continue;  /* NS positions 0-48 */
       kuka= ( prof.j_NS < 12 || prof.j_NS > 36 ) ? 0 : 1;

       if ( kuka == 1 ) {
         /*if ( prof.j_NS < 12 || prof.j_NS > 36 ) continue;*/  /* NS positions 12-36 correspond to MS positions 0-24 */
       }

       /*--- Reject bad entries ----*/
       nbad= 0;
       for (i= 0; i< NCHAN; i++) {
         if ( i < nc_S1 && ( prof.tb[i] < 20.0 || prof.tb[i] > 350 )) nbad++;
         if ( prof.tb[i] < 20.0 || prof.tb[i] > 350 ) prof.tb[i]= -1.0; /* for proper print format */
       }
       if ( prof.sfc_class < 0 ) nbad++;
       if ( prof.yyyy < 1900 || prof.yyyy > 2100 ) nbad++;
       if ( prof.mm < 1 || prof.mm > 12 ) nbad++;
       if ( prof.dd < 1 || prof.dd > 31 ) nbad++;
       if ( prof.hh < 0 || prof.hh > 23 ) nbad++;
       if ( prof.mn < 0 || prof.mn > 59 ) nbad++;
       if ( prof.ss < 0 || prof.ss > 59 ) nbad++;
       if ( fabs(prof.glat1) > 90 || fabs(prof.glon1) > 180 ) nbad++;
       if ( prof.nku10 == 0 && prof.nku15 == 0 && prof.nku20 == 0 && prof.nku25 == 0 ) nbad++;
       if ( prof.precip_NS < 0.0 || prof.precip_NS_cmb < 0.0 ) nbad++;
       if ( kuka == 1 ) {
         if ( prof.nka10 == 0 && prof.nka15 == 0 && prof.nka20 == 0 && prof.nka25 == 0 ) nbad++;
         if ( prof.precip_MS < 0.0 || prof.precip_MS_cmb < 0.0 ) nbad++;
       }
       if ( prof.ts < 100 || prof.ts > 400 ) nbad++;
       if ( prof.t2m < 100 || prof.t2m > 400 ) nbad++;
       if ( prof.tqv < 0.001 || prof.tqv > 100 ) nbad++;

       if ( nbad > 0 ) {
         /*printf("BAD DPR %6d %5d %3d %8.3f %8.3f ", prof.rev, prof.i_NS, prof.j_NS, prof.glat1, prof.glon1);
         printf("%04d/%02d/%02d %02d:%02d:%02d ", prof.yyyy, prof.mm, prof.dd, prof.hh, prof.mn, prof.ss);
         for (k= 0; k< 3; k++) printf("%8.3f ", prof.pc_emis[k]);
         for (k= 0; k< nc_S1; k++) printf("%6.2f ", prof.tb[k]);
         printf("%3d %6.2f %6.2f %6.2f ", prof.sfc_class, prof.ts, prof.t2m, prof.tqv);
         printf("%4d %4d %4d %4d ", prof.nku10, prof.nku15, prof.nku20, prof.nku25);
         if ( kuka == 1 ) printf("%4d %4d %4d %4d ", prof.nka10, prof.nka15, prof.nka20, prof.nka25);
         printf("%8.2f %8.2f ", prof.precip_NS, prof.precip_NS_cmb);
         if ( kuka == 1 ) printf("%8.2f %8.2f ", prof.precip_MS, prof.precip_MS_cmb);
         printf("\n");*/
         continue;
       }
       /*---------------------------*/

       if ( prof.precip_NS_cmb > 1 ) nrain++;

       prof_idx[nrec1]= prof;
       nrec1++;
       if ( kuka == 1 ) nrec2++;

       if ( nrec1 == maxrec ) {
         printf("Exceeded max records= %d %d\n", nrec1, nrec2);
         break;
       }

     }
     fclose(fdb2);
     printf("%s  NRead=%d  Nrec=%d %d Nrain=%d    Npix using this DB file=%d\n", dbfile, nrec, nrec1, nrec2, nrain, n_index[idx_db]);

     if ( nrec1 == 0 || nrec2 == 0 ) continue;
     /*-----------------------------------------------------------------------*/

     nsat= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {

       t1= secs_S1_buf[isat];
       tm= *gmtime(&t1);
       sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

       for (jsat=0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx_out= nsat*ns_S1 + jsat;
           
         if ( pixel_index[idx] != idx_db ) continue;
         ndb++;
         pcent_db= 100.0*ndb/n_index[idx_db];

         glat1= lat_S1_buf[idx];
         glon1= lon_S1_buf[idx];
         if ( fabs(glat1) > 90 || fabs(glon1) > 180 ) continue;

         if ( have_model > 0 ) {
           ts_MERRA2= ts_buf[idx];
           t2m_MERRA2= t2m_buf[idx];
           tqv_MERRA2= tqv_buf[idx];
           h273_MERRA2= h273_buf[idx];
         }
         else {
           ts_MERRA2= t2m_MERRA2= tqv_MERRA2= -1.0;
           h273_MERRA2= -9999;
         }

         if ( have_gprof > 0 ) {
           sfc_class= sfc_GPROF_buf[idx];
           precip_GPROF= surfacePrecipitation_GPROF_buf[idx];
           prob_precip_GPROF= probabilityOfPrecip_GPROF_buf[idx];
           tqv_GPROF= totalColumnWaterVaporIndex_GPROF_buf[idx];
           t2m_GPROF= temp2mIndex_GPROF_buf[idx];
           frozen_precip_GPROF= frozenPrecipitation_GPROF_buf[idx];
           if ( precip_GPROF < 0.0 ) precip_GPROF= -1.0;
           if ( frozen_precip_GPROF < 0.0 ) frozen_precip_GPROF= -1.0;
           if ( prob_precip_GPROF < 0.0 || prob_precip_GPROF > 100.0 ) prob_precip_GPROF= -1.0;
         }
         else {
           sfc_class= tqv_GPROF= t2m_GPROF= -1;
           precip_GPROF= prob_precip_GPROF= frozen_precip_GPROF= -1.0;
         }


         nbad= 0;
         for (k= 0; k< nc_S1+nc_S2; k++) {
           idx1= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;
           tb_gmi[k]= tb_S1S2_buf[idx1];
           if ( k < nc_S1 && ( tb_gmi[k] < 20 || tb_gmi[k] > 350 )) nbad++;
         }
         if ( nbad != 0 ) {
           printf("bad TB\n");
           continue;
         }

         for (k=0; k < NEM; k++) {
           /*idx3= isat*ns_S1*NEM + jsat*NEM + k;*/
           idx3= nsat*ns_S1*NEM + jsat*NEM + k;
           emis[k]= emis_S1S2_buf[idx3];
           pc_emis[k]= pc_emis_S1S2_buf[idx3];
         }

         /* START LOOKUP */

         printf("\n\n");
         x= -1.0;
         ix= -1;
         printf("OBS %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, isatname, irev, tdate, isat, jsat, glat1, glon1);
         for (k= 0; k< 3; k++) printf("%8.3f ", pc_emis[k]);
         for (k= 0; k< 9; k++) printf("%6.2f ", tb_gmi[k]);
         printf("%3d ", sfc_class);
         printf("%7.2f %7.2f %7.2f %6d ", ts_MERRA2, t2m_MERRA2, tqv_MERRA2, h273_MERRA2);
         printf("%4d %4d ", ix, ix);
         printf("%7.2f %7.2f %7.2f ", x,x,x);
         printf("%6d ", ix);
         e0= ( emis[0] > 0.0 && emis[0] < 2.0 ) ? emis[0] : -1.0;
         e1= ( emis[1] > 0.0 && emis[1] < 2.0 ) ? emis[1] : -1.0;
         printf("%6.3f %6.3f ", e0, e1);
         printf("%6.3f %6.3f ", x, x );
         printf("\n");

         rmsd_min= rmsd2_min= 1.0E9;
         rmsd_max= rmsd2_max= 0.0;
         nrec2= 0;
         for (i= 0; i< nrec1; i++) {

           /* If T2m is warm, reject entries where T2m is cold, and vice-versa */
           if ( have_model > 0 ) {
             if ( fabs(prof_idx[i].t2m - t2m_MERRA2) > MAX_T2M_DIFF ) continue;
             /*if (( fabs(prof_idx[i].t2m - t2m_MERRA2) > MAX_T2M_DIFF ) || ( fabs(prof_idx[i].tqv - tqv_MERRA2) > MAX_TQV_DIFF )) continue;*/
             if (( prof_idx[i].t2m < 278 && t2m_MERRA2 > 278 ) || ( prof_idx[i].t2m > 278 && t2m_MERRA2 < 278 )) { 
               /*printf("OBS= %8.3f %8.3f %8.3f   MODEL= %8.3f %8.3f %8.3f\n", ts_MERRA2, t2m_MERRA2, tqv_MERRA2, prof_idx[i].ts, prof_idx[i].t2m, prof_idx[i].tqv);*/
               continue;
             }
           }

           rmsd1= rmsd2= 0.0;
           nrmsd1= nrmsd2= 0;
           for (k= 0; k< NCHAN; k++) {
             x= (tb_gmi[k] - prof_idx[i].tb[k])/std_tb[k];

             if ( k < nc_S1 ) {
               rmsd2+= x*x;
               nrmsd2++;
             }
             else {
               rmsd1+= x*x;
               nrmsd1++;
             }
           }
           if ( nrmsd2 < nc_S1 ) {
             printf("ERROR in wt assignment - TB out of range N=%d\n", nrmsd2);
             for (k= 0; k< NCHAN; k++) printf("%7.2f ", tb_gmi[k]); printf("\n");
             for (k= 0; k< NCHAN; k++) printf("%7.2f ", prof_idx[i].tb[k]); printf("\n");
             exit(1);
           }
           rmsd2= sqrt(rmsd2/nrmsd2);
           rmsd1= ( nrmsd1 > 0 ) ? sqrt(rmsd1/nrmsd1) : -1.0;

           rmsd= 0.0;
           nem_compare= 0;
           for (k= 0; k< NEM; k++) {
             if ( em_compare[k] != 1 ) continue;
             x= (pc_emis[k] - prof_idx[i].pc_emis[k])/std_pc[k];
             rmsd+= x*x;
             nem_compare++;
           }
           if ( irev == prof_idx[i].rev || rmsd < 1.0E-12 ) continue;
           rmsd= sqrt(rmsd/nem_compare);

           wt= sqrt(0.0*rmsd2*rmsd2 + 1.0*rmsd*rmsd);
           /*
           printf("wt= %f\n",wt);
           */ 

           /*wt= sqrt(1.0*rmsd2*rmsd2 + 0.0*rmsd*rmsd);*/
           /*wt= sqrt(0.7*rmsd1*rmsd1 + 0.3*rmsd*rmsd);*/
           /*wt= sqrt(0.3*rmsd2*rmsd2 + 0.7*rmsd*rmsd);*/

           if ( wt >= INT_MAX/1E6 ) {
             printf("Exceeded INT_MAX=%d  rmsd=%f rmsd2=%f wt=%f\n", INT_MAX, rmsd, rmsd2, wt);
             exit(1);
           }
           tmp1[nrec2]= 1E6*wt;
           tmp2[nrec2]= i;
           tmp3[i]= wt;
           tmp7[nrec2]= wt;
           if ( rmsd > rmsd_max ) rmsd_max= rmsd;
           if ( rmsd < rmsd_min ) rmsd_min= rmsd;
           if ( rmsd2 > rmsd2_max ) rmsd2_max= rmsd2;
           if ( rmsd2 < rmsd2_min ) rmsd2_min= rmsd2;
           nrec2++;

         }
         /*printf("TB rmsd min max= %10.5f %10.5f ", rmsd_min, rmsd_max);
         printf("EPC rmsd min max= %10.5f %10.5f\n", rmsd2_min, rmsd2_max);*/

         heapsort2(tmp1, tmp2, nrec2);
         /*for (i= 0; i< nrec2; i++) {
           k2= tmp2[i];
           if ( i < 10 ) printf("%6d %6d   %d %f   %d %f\n", i, k2, tmp1[i], tmp7[i], tmp1[k2], tmp7[k2]);
         }*/


         /*--- First the NS database then the MS database ---*/

         for (idb= 0; idb< 2; idb++) {

           if ( idb == 0 ) printf("NS database  N=%d %8.3f\n", nrec2, pcent_db);
           else printf("MS database  N=%d %8.3f\n", nrec2, pcent_db);

           zmax= 0.0;
           n1= n2= 0;
           precip_sum= precip2_sum= ts_sum= t2m_sum= tqv_sum= h273_sum= ku_ztop_sum= ka_ztop_sum= 0.0;
           wt_precip_sum= wt_precip2_sum= wt_ts_sum= wt_t2m_sum= wt_tqv_sum= wt_h273_sum= wt_ku_ztop_sum= wt_ka_ztop_sum= 0.0;
           n_precip_sum= n_precip2_sum= n_precip_all_sum= n_precip2_all_sum= 0;
           for (iclass= 0; iclass< NCLASS; iclass++) class[iclass]= 0;
           for (iprecip= 0; iprecip< NPRECIP; iprecip++) hist_precip[iprecip]= 0;
           for (k= 0; k< NLEV_NS; k++) {
             z_ku[k]= z_ka[k]= 0.0;
             nz_ku[k]= nz_ka[k]= 0;
             for (k1= 0; k1< NZ; k1++) {
               ccfads_ku[k][k1]= 0;
               ccfads_ka[k][k1]= 0;
             }
           }

           wt0= -1.0;
           nfound= 0;
           nku20= nku25= 0;

           for (i= 0; i< nrec2; i++) {
             k2= tmp2[i];
             /*wt= tmp7[k2];*/
             wt= tmp3[k2];
             /*printf("i=%d  k2=%d  wt=%f.2\n",i,k,wt);*/
             kuka= ( prof_idx[k2].j_NS < 12 || prof_idx[k2].j_NS > 36 ) ? 0 : 1;
             if ( idb == 1 && kuka == 0 ) continue;

             if ( wt0 < 0.0 ) {
               wt0= wt;
               wt2= 1.0;
             }
             else {
               wt2= exp(-0.5*(wt/wt0)*(wt/wt0));
             }
             /*printf("i=%d  wt %f.3  wt0 %f.3  wt2= %f.3 \n", i, wt, wt0, wt2);*/
             /* Freezing level height */
             h273= -9999;
             for (j1= 0; j1 < 42; j1++) {
               /*printf("%2d ModelH= %9.4f  ModelP=%9.4f  ModelT=%9.4f\n", j1, prof_idx[k2].h_prof[j1], 0.1*prof_idx[k2].p_prof[j1], prof_idx[k2].t_prof[j1]);*/
               if ( j1 == 0 ) continue;
               if ( prof_idx[k2].h_prof[j1] < -2000 || prof_idx[k2].h_prof[j1] > 10000 ) continue;
               if ( prof_idx[k2].h_prof[j1-1] < -2000 || prof_idx[k2].h_prof[j1-1] > 10000 ) continue;
               if ( prof_idx[k2].t_prof[j1-1] < 273 && prof_idx[k2].t_prof[j1] >= 273 ) {
                 frac= ( 273.0 - prof_idx[k2].t_prof[j1-1])/(prof_idx[k2].t_prof[j1] - prof_idx[k2].t_prof[j1-1]);
                 h273= prof_idx[k2].h_prof[j1-1] + frac*(prof_idx[k2].h_prof[j1] - prof_idx[k2].h_prof[j1-1]);
                 /*printf("  FREEZLEV  ModelH=%9.4f %9.4f  ModelT=%9.4f %9.4f Frac=%9.4f H273=%9.4f\n",
                   prof_idx[k2].h_prof[j1], prof_idx[k2].h_prof[j1-1], prof_idx[k2].t_prof[j1], prof_idx[k2].t_prof[j1-1], frac, h273);*/
                 break;
               }
             }
             /*if ( h273 < -2000 || h273 > 9000 ) {
               printf("Freezing level height out of bounds %f\n", h273);
               exit(1);
             }*/
             elev= prof_idx[k2].elev;

             if ( idb == 1 ) {
               precip= prof_idx[k2].precip_MS_cmb;
               precip2= prof_idx[k2].precip_MS;
             }
             else {
               precip= prof_idx[k2].precip_NS_cmb;
               precip2= prof_idx[k2].precip_NS;
             }
             if ( precip < 0.0 || precip2 < 0.0 ) {
               printf("Rainrates out of bounds %f %f\n", precip, precip2);
               exit(1);
             }

             precip_sum+= wt2*precip;
             precip2_sum+= wt2*precip2;
             wt_precip_sum+= wt2;
             if ( precip > 0.0 ) n_precip_sum++;
             if ( precip2 > 0.0 ) n_precip2_sum++;
             n_precip_all_sum++;

             ts= prof_idx[k2].ts;
             ts_sum+= wt2*ts;

             t2m= prof_idx[k2].t2m;
             t2m_sum+= wt2*t2m;

             tqv= prof_idx[k2].tqv;
             tqv_sum+= wt2*tqv;

             if ( h273 > -2000 && h273 < 9000 ) {
               h273_sum+= wt2*h273;
               wt_h273_sum+= wt2;
             }

             iclass= prof_idx[k2].sfc_class;
             if ( iclass >= 0 && iclass < NCLASS ) class[iclass]++;

             /* Store top-ranked Ku/Ka profiles and TB */
             if ( wt2 == 1.0 ) {
               for (k= 0; k< NLEV_NS; k++) {
                 /*idx1= isat*ns_S1*NLEV_NS + jsat*NLEV_NS + k;*/
                 idx1= nsat*ns_S1*NLEV_NS + jsat*NLEV_NS + k;
                 if ( idb == 0 ) {
                   z_ku_1_NS_buf[idx1]= prof_idx[k2].z_ku[k];
                 }
                 else {
                   z_ku_1_MS_buf[idx1]= prof_idx[k2].z_ku[k];
                   z_ka_1_MS_buf[idx1]= prof_idx[k2].z_ka[k];
                 }
               }
               for (k=0; k< nc_S1+nc_S2; k++) {
                 /*idx1= isat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;*/
                 idx1= nsat*ns_S1*(nc_S1+nc_S2) + jsat*(nc_S1+nc_S2) + k;
                 if ( idb == 0 ) 
                   tb_1_NS_buf[idx1]= prof_idx[k2].tb[k];
                 else
                   tb_1_MS_buf[idx1]= prof_idx[k2].tb[k];
               }
             }

             /* CCFADS of Ku and Ka Z profiles */
             ku_clutter= ka_clutter= -1;
             for (k= 0; k< NLEV_NS; k++) {

               zlog= 0.01*prof_idx[k2].z_ku[k];
               if ( zlog < 0 && zlog > -90 ) {
                 zlog*= -1.0;  /* ground clutter */
                 if ( ku_clutter < 0 ) ku_clutter= k;
               }
               if ( zlog > 0 ) {
                 z_ku[k]+= pow(10.0, 0.1*zlog);
               }
               nz_ku[k]++;
               if ( ku_clutter < 0 ) {
                 if ( k > 70 && zlog > zmax ) zmax= zlog;
                 iz= NZ*(zlog-ZMIN)/(ZMAX-ZMIN);
                 if ( iz >= 0 && iz < NZ ) ccfads_ku[k][iz]++;
               }

               zlog= 0.01*prof_idx[k2].z_ka[k];
               if ( zlog < 0 && zlog > -90 ) {
                 zlog*= -1.0;  /* ground clutter */
                 if ( ka_clutter < 0 ) ka_clutter= k;
               }
               if ( zlog > 0 ) {
                 z_ka[k]+= pow(10.0, 0.1*zlog);
               }
               nz_ka[k]++;
               if ( ka_clutter < 0 ) {
                 if ( k > 70 && zlog > zmax ) zmax= zlog;
                 iz= NZ*(zlog-ZMIN)/(ZMAX-ZMIN);
                 if ( iz >= 0 && iz < NZ ) ccfads_ka[k][iz]++;
               }
             }


             /* 20-dB top of Ku and Ka Z profiles */
             if ( ku_clutter > 0 ) {
               ku_ztop= -1.0;
               for (k= 1; k< ku_clutter-3; k++) {
                 ht= 0.25*(NLEV_NS-k-1);
                 zlog1= 0.01*prof_idx[k2].z_ku[k-1];
                 zlog2= 0.01*prof_idx[k2].z_ku[k];
                 zlog3= 0.01*prof_idx[k2].z_ku[k+1];
                 if ( zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) {
                   ku_ztop= ht+0.25;
                   break;
                 }
               }
               if ( ku_ztop > 0 ) ku_ztop_sum+= wt2*ku_ztop;
               wt_ku_ztop_sum+= wt2;
             }
             if ( ka_clutter > 0 ) {
               ka_ztop= -1.0;
               for (k= 1; k< ka_clutter-3; k++) {
                 ht= 0.25*(NLEV_NS-k-1);
                 zlog1= 0.01*prof_idx[k2].z_ka[k-1];
                 zlog2= 0.01*prof_idx[k2].z_ka[k];
                 zlog3= 0.01*prof_idx[k2].z_ka[k+1];
                 if ( zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) {
                   ka_ztop= ht+0.25;
                   break;
                 }
               }
               if ( ka_ztop > 0 ) ka_ztop_sum+= wt2*ka_ztop;
               wt_ka_ztop_sum+= wt2;
             }

             if ( nfound < 10 ) {
               nku20+= prof_idx[k2].nku20;
               nku25+= prof_idx[k2].nku25;
             }

             /*  PDF precip */
             if ( precip < 0.1 ) iprecip= 0;
             else iprecip= NPRECIP*(precip-PRECIP1)/(PRECIP2-PRECIP1) + 1;
             /*printf("precip=%8.3f iprecip=%d\n", precip, iprecip);*/
             if ( iprecip >= 0 && iprecip < NPRECIP ) hist_precip[iprecip]++;

             sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d", prof_idx[k2].yyyy, prof_idx[k2].mm, prof_idx[k2].dd, prof_idx[k2].hh, prof_idx[k2].mn, prof_idx[k2].ss);

             /* Also write to output file for scatter plots */

             if ( nfound < 10 ) {
               printf("DB  %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, prof_idx[k2].satid, prof_idx[k2].rev, pdate, prof_idx[k2].j_S1, prof_idx[k2].j_NS, prof_idx[k2].glat1, prof_idx[k2].glon1);
               for (k= 0; k< 3; k++) printf("%8.3f ", prof_idx[k2].pc_emis[k]);
               for (k= 0; k< 9; k++) printf("%6.2f ", prof_idx[k2].tb[k]);
               printf("%3d ", prof_idx[k2].sfc_class);
               printf("%7.2f %7.2f %7.2f %6d ", prof_idx[k2].ts, prof_idx[k2].t2m, prof_idx[k2].tqv, (int) h273);
               printf("%4d %4d ", prof_idx[k2].nku20, prof_idx[k2].nku25);
               printf("%7.2f ", prof_idx[k2].precip_GPROF);
               printf("%7.2f %7.2f %6d ", precip, precip2, elev);
               e0= e1= -1.0;
               if ( prof_idx[k2].nku20 == 0 && prof_idx[k2].nku25 == 0 ) {
                 e0= ( prof_idx[k2].emis[0] > 0.0 ) ? prof_idx[k2].emis[0] : -1.0;
                 e1= ( prof_idx[k2].emis[1] > 0.0 ) ? prof_idx[k2].emis[1] : -1.0;
               }
               printf("%6.3f %6.3f ", e0, e1);
               e0= ( prof_idx[k2].emis_cmb[0] > 0.0 && prof_idx[k2].emis_cmb[0] < 2.0 ) ? prof_idx[k2].emis_cmb[0] : -1.0;
               e1= ( prof_idx[k2].emis_cmb[1] > 0.0 && prof_idx[k2].emis_cmb[1] < 2.0 ) ? prof_idx[k2].emis_cmb[1] : -1.0;
               printf("%6.3f %6.3f ", e0, e1);
               printf("%7.4f ", wt2);
               printf("%5.2f ", ku_ztop);
               if ( idb == 1 ) printf("%5.2f ", ka_ztop);
               printf("\n");
             }

             nfound++;
             nfound_all++;

             if ( wt2 < WT2MIN || n_precip_all_sum >= N2MAX ) {
               if ( idb == 0 ) qual_NS_buf[idx_out]= nfound;
               else qual_MS_buf[idx_out]= nfound;
               break;
             }

           }  /* nrec2 loop */

           /*--- All done with lookup ---*/

           if ( n_precip_all_sum == 0 ) {
             printf("No database entries located\n");
             continue;
           }
           precip= precip_sum/wt_precip_sum;
           if ( idb == 0 ) precip_NS_buf[idx_out]= precip;
           else precip_MS_buf[idx_out]= precip;

           precip2= precip2_sum/wt_precip_sum;
           if ( idb == 0 ) precip2_NS_buf[idx_out]= precip2;
           else precip2_MS_buf[idx_out]= precip2;

           prob_precip= 100.0*n_precip_sum/n_precip_all_sum;
           if ( idb == 0 ) precip_prob_NS_buf[idx_out]= prob_precip;
           else precip_prob_MS_buf[idx_out]= prob_precip;

           prob_precip2= 100.0*n_precip2_sum/n_precip_all_sum;
           if ( idb == 0 ) precip2_prob_NS_buf[idx_out]= prob_precip2;
           else precip2_prob_MS_buf[idx_out]= prob_precip2;
 
           ts= ts_sum/wt_precip_sum;
           if ( idb == 0 ) ts_epc_NS_buf[idx_out]= ts;
           else ts_epc_MS_buf[idx_out]= ts;

           t2m= t2m_sum/wt_precip_sum;
           if ( idb == 0 ) t2m_epc_NS_buf[idx_out]= t2m;
           else t2m_epc_MS_buf[idx_out]= t2m;

           tqv= tqv_sum/wt_precip_sum;
           if ( idb == 0 ) tqv_epc_NS_buf[idx_out]= tqv;
           else tqv_epc_MS_buf[idx_out]= tqv;

           if ( wt_h273_sum > 0.0 ) 
             h273= h273_sum/wt_h273_sum;
           else
             h273= -9999;
           if ( idb == 0 ) h273_epc_NS_buf[idx_out]= h273;
           else h273_epc_MS_buf[idx_out]= h273;

           if ( wt_ku_ztop_sum > 0.0 ) {
             ku_ztop= ku_ztop_sum/wt_ku_ztop_sum;
             if ( idb == 0 ) ku_ztop_epc_NS_buf[idx_out]= 1000*ku_ztop;
             else ku_ztop_epc_MS_buf[idx_out]= 1000*ku_ztop;
           }
           else
             ku_ztop= -1.0;

           if ( idb == 1 ) {
             if ( wt_ka_ztop_sum > 0.0 ) {
               ka_ztop= ka_ztop_sum/wt_ka_ztop_sum;
               ka_ztop_epc_MS_buf[idx_out]= 1000*ka_ztop;
             }
           }
           else
             ka_ztop= -1.0;


           printf("EST %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, isatname, irev, tdate, isat, jsat, glat1, glon1);
           for (k= 0; k< 3; k++) printf("%8.3f ", pc_emis[k]);
           for (k= 0; k< 9; k++) printf("%6.2f ", tb_gmi[k]);
           printf("%3d ", sfc_class);
           printf("%7.2f %7.2f %7.2f ", ts, t2m, tqv);
           printf("%6d ", (int) h273);
           printf("%4d %4d ", nku20, nku25);
           printf("%7.2f ", precip_GPROF);
           printf("%7.2f %7.2f ", precip, precip2);
           printf("%7.2f %7.2f %7.2f ", prob_precip_GPROF, prob_precip, prob_precip2);
           printf("%5.2f ", ku_ztop);
           if ( idb == 1 ) printf("%5.2f ", ka_ztop);
           printf("\n");

           if ( zmax > 40 && precip > 25 ) {
             for (k= 0; k< NLEV_NS; k++) {
               zlog1= ( z_ku[k] > 0 && nz_ku[k] > 0 ) ? 10.0*log10(z_ku[k]/nz_ku[k]) : -99.99;
               zlog2= ( z_ka[k] > 0 && nz_ka[k] > 0 ) ? 10.0*log10(z_ka[k]/nz_ka[k]) : -99.99;
               printf("%3d %6d %6.2f %6d %6.2f\n", k, nz_ku[k], zlog1, nz_ka[k], zlog2);
             }
             for (iprecip= 0; iprecip< NPRECIP; iprecip++) {
               x= ( iprecip == 0 ) ? 0.0 : PRECIP1 + iprecip*DPRECIP;
               if ( hist_precip[iprecip] > 0 ) printf("%3d %6.2f N=%d\n", iprecip, x, hist_precip[iprecip]);
             }
           }

         }   /* idb loop for NS or MS retrievals */

         npix++;

       }  /* next GMI sample */

       nsat++;  /* increment output buffer line index */

     }  /* next GMI line */

     printf("\n%s KuKa=%2d DBidx=%5d Frac=%6.4f %6.4f Nread=%7d Nrain=%7d Npix=%6d  Percent complete=%7.3f\n", tbfile, idb, idx_db, frac0, frac1, nrec, nrain, ndb, 100.0*npix/((isat_end-isat_start+1)*ns_S1) );

   }  /* next DB file */

   if ( plot_files == 1 && nfound_all > 0 ) {
     fclose(fout);
     fclose(fret);
   }

   /*-------------------------------------------------------------------------*/
   /*  output in netCDF */

   nc_out= nc_S1+nc_S2;

   now= time(NULL);
   t0= now;
   tm= *gmtime(&t0);
   sprintf(cdate_end,"%4d/%02d/%02d %02d:%02d:%02d", 1900+tm.tm_year, 1+tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

   fbuf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   fbuf3= (float *) malloc(nl_out*ns_S1*nc_out*sizeof(float));
   fbuf4= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));
   ibuf= (int *) malloc(nl_out*ns_S1*sizeof(int));
   sbuf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   bbuf= (char *) malloc(nl_out*ns_S1*sizeof(char));
   dbuf= (double *) malloc(nl_out*sizeof(double));

   if ((retval = nc_create(outfile, NC_NETCDF4, &ncid))) ERR(retval);
   printf("Opened %s\n", outfile);

   sprintf(ctmp,"%s", argv[1]);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "1B.GPM.GMI", strlen(ctmp), ctmp))) ERR(retval);
   if ( have_gprof > 0 ) {
     sprintf(ctmp,"%s", gprof_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "2A.GMI.GPROF", strlen(ctmp), ctmp))) ERR(retval);
   }
   if ( have_model > 0 ) {
     sprintf(ctmp,"%s", model_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "MERRA2_file", strlen(ctmp), ctmp))) ERR(retval);
   }
   sprintf(ctmp,"%f", clat);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_lat", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", clon);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_lon", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_date", strlen(cdate), cdate))) ERR(retval);
   sprintf(ctmp,"%06d", irev);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "orbit_rev", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "satname", strlen(satname), satname))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "sensor", strlen(sensor), sensor))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "processing_start", strlen(cdate_start), cdate_start))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "processing_end", strlen(cdate_end), cdate_end))) ERR(retval);

   /*if ( kuka == 1 ) sprintf(ctmp,"%s", "MS");
   else sprintf(ctmp,"%s", "NS");
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "NS_or_MS", strlen(ctmp), ctmp))) ERR(retval);*/

   sprintf(ctmp,"%d", DB_MAXREC);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "max_database_records", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%d", N2MAX);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "max_sorted_database_records", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", WT2MIN);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "min_weight_threshold", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", DB_RAINFRAC);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "min_database_rain_fraction", strlen(ctmp), ctmp))) ERR(retval);
   if ( have_model > 0 ) {
     sprintf(ctmp,"%d", MAX_T2M_DIFF);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "max_t2m_difference", strlen(ctmp), ctmp))) ERR(retval);
     sprintf(ctmp,"%d", MAX_TQV_DIFF);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "max_tqv_difference", strlen(ctmp), ctmp))) ERR(retval);
   }

   if ((retval = nc_def_dim(ncid, "nscan", nl_out, &scan_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "npix", ns_S1, &pix_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "ntb", nc_out, &tb_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "nem", NEM, &em_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "nlev", NLEV_NS, &lev_dimid))) ERR(retval);

   if ((retval = nc_def_grp(ncid, "NS", &grpid_NS))) ERR (retval);
   if ((retval = nc_def_grp(ncid, "MS", &grpid_MS))) ERR (retval);
   if ( have_gprof > 0 ) if ((retval = nc_def_grp(ncid, "GPROF", &grpid1))) ERR (retval);
   if ( have_model > 0 ) if ((retval = nc_def_grp(ncid, "MERRA2", &grpid2))) ERR (retval);

   if ((retval = nc_def_var(ncid, "time", NC_DOUBLE, 1, &scan_dimid, &secs_varid))) ERR(retval);
   sprintf(units,"seconds since 1 Jan 1970");
   if ((retval = nc_put_att_text(ncid, secs_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(ncid, "epc_compare", NC_SHORT, 1, &em_dimid, &nem_compare_varid))) ERR(retval);
   sprintf(units,"Set to 1 for EPC elements used in distance calculation");
   if ((retval = nc_put_att_text(ncid, nem_compare_varid, "units", strlen(units), units))) ERR(retval);

   dimids2[0] = scan_dimid;
   dimids2[1] = pix_dimid;

   /* latitude longitude */
   if ((retval = nc_def_var(ncid, "latitude", NC_FLOAT, 2, dimids2, &lat_varid))) ERR(retval);
   sprintf(units,"degrees");
   if ((retval = nc_put_att_text(ncid, lat_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(ncid, "longitude", NC_FLOAT, 2, dimids2, &lon_varid))) ERR(retval);
   sprintf(units,"degrees");
   if ((retval = nc_put_att_text(ncid, lon_varid, "units", strlen(units), units))) ERR(retval);

   /*---- estimates with the NS database ----*/

   /* mean precip based on combined */
   if ((retval = nc_def_var(grpid_NS, "precip", NC_FLOAT, 2, dimids2, &precip_NS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR+GMI combined estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_NS, precip_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean precip based on radar-only */
   if ((retval = nc_def_var(grpid_NS, "precip2", NC_FLOAT, 2, dimids2, &precip2_NS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_NS, precip2_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip */
   if ((retval = nc_def_var(grpid_NS, "precip_prob", NC_FLOAT, 2, dimids2, &precip_prob_NS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR+GMI combined estimate, in percent");
   if ((retval = nc_put_att_text(grpid_NS, precip_prob_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip2 */
   if ((retval = nc_def_var(grpid_NS, "precip2_prob", NC_FLOAT, 2, dimids2, &precip2_prob_NS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR estimate, in percent");
   if ((retval = nc_put_att_text(grpid_NS, precip2_prob_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* quality flag */
   if ((retval = nc_def_var(grpid_NS, "quality", NC_INT, 2, dimids2, &quality_NS_varid))) ERR(retval);
   sprintf(units,"0=good 1=bad cooords 2=no DBfile located");
   if ((retval = nc_put_att_text(grpid_NS, quality_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean Ts */
   if ((retval = nc_def_var(grpid_NS, "ts", NC_FLOAT, 2, dimids2, &ts_NS_varid))) ERR(retval);
   sprintf(units,"estimated surface temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_NS, ts_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean T2m */
   if ((retval = nc_def_var(grpid_NS, "t2m", NC_FLOAT, 2, dimids2, &t2m_NS_varid))) ERR(retval);
   sprintf(units,"estimated 2-m air temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_NS, t2m_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean column vapor  */
   if ((retval = nc_def_var(grpid_NS, "tqv", NC_FLOAT, 2, dimids2, &tqv_NS_varid))) ERR(retval);
   sprintf(units,"estimated total column water vapor, in mm");
   if ((retval = nc_put_att_text(grpid_NS, tqv_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height 273K */
   if ((retval = nc_def_var(grpid_NS, "h273", NC_SHORT, 2, dimids2, &h273_NS_varid))) ERR(retval);
   sprintf(units,"estimated freezing level height, in meters");
   if ((retval = nc_put_att_text(grpid_NS, h273_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ku-band 20 db level */
   if ((retval = nc_def_var(grpid_NS, "h20dB_ku", NC_SHORT, 2, dimids2, &ku_ztop_NS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ku-band, in meters");
   if ((retval = nc_put_att_text(grpid_NS, ku_ztop_NS_varid, "units", strlen(units), units))) ERR(retval);

   /*---- estimates with the MS database ----*/

   /* mean precip based on combined */
   if ((retval = nc_def_var(grpid_MS, "precip", NC_FLOAT, 2, dimids2, &precip_MS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR+GMI combined estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_MS, precip_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean precip based on radar-only */
   if ((retval = nc_def_var(grpid_MS, "precip2", NC_FLOAT, 2, dimids2, &precip2_MS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_MS, precip2_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip */
   if ((retval = nc_def_var(grpid_MS, "precip_prob", NC_FLOAT, 2, dimids2, &precip_prob_MS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR+GMI combined estimate, in percent");
   if ((retval = nc_put_att_text(grpid_MS, precip_prob_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip2 */
   if ((retval = nc_def_var(grpid_MS, "precip2_prob", NC_FLOAT, 2, dimids2, &precip2_prob_MS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR estimate, in percent");
   if ((retval = nc_put_att_text(grpid_MS, precip2_prob_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* quality flag */
   if ((retval = nc_def_var(grpid_MS, "quality", NC_INT, 2, dimids2, &quality_MS_varid))) ERR(retval);
   sprintf(units,"0=good 1=bad cooords 2=no DBfile located");
   if ((retval = nc_put_att_text(grpid_MS, quality_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean Ts */
   if ((retval = nc_def_var(grpid_MS, "ts", NC_FLOAT, 2, dimids2, &ts_MS_varid))) ERR(retval);
   sprintf(units,"estimated surface temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_MS, ts_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean T2m */
   if ((retval = nc_def_var(grpid_MS, "t2m", NC_FLOAT, 2, dimids2, &t2m_MS_varid))) ERR(retval);
   sprintf(units,"estimated 2-m air temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_MS, t2m_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean column vapor  */
   if ((retval = nc_def_var(grpid_MS, "tqv", NC_FLOAT, 2, dimids2, &tqv_MS_varid))) ERR(retval);
   sprintf(units,"estimated total column water vapor, in mm");
   if ((retval = nc_put_att_text(grpid_MS, tqv_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height 273K */
   if ((retval = nc_def_var(grpid_MS, "h273", NC_SHORT, 2, dimids2, &h273_MS_varid))) ERR(retval);
   sprintf(units,"estimated freezing level height, in meters");
   if ((retval = nc_put_att_text(grpid_MS, h273_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ku-band 20 db level */
   if ((retval = nc_def_var(grpid_MS, "h20dB_ku", NC_SHORT, 2, dimids2, &ku_ztop_MS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ku-band, in meters");
   if ((retval = nc_put_att_text(grpid_MS, ku_ztop_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ka-band 20 db level */
   if ((retval = nc_def_var(grpid_MS, "h20dB_ka", NC_SHORT, 2, dimids2, &ka_ztop_MS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ka-band, in meters");
   if ((retval = nc_put_att_text(grpid_MS, ka_ztop_MS_varid, "units", strlen(units), units))) ERR(retval);


   if ( have_gprof > 0 ) {
     /* GPROF precip */
     if ((retval = nc_def_var(grpid1, "precip", NC_FLOAT, 2, dimids2, &precip_GPROF_varid))) ERR(retval);
     sprintf(units,"estimate from GPROF, in mm/hr");
     if ((retval = nc_put_att_text(grpid1, precip_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF prob of precip */
     if ((retval = nc_def_var(grpid1, "precip_prob", NC_FLOAT, 2, dimids2, &precip_prob_GPROF_varid))) ERR(retval);
     sprintf(units,"probability of precip from GPROF, in percent");
     if ((retval = nc_put_att_text(grpid1, precip_prob_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF T2m index */
     if ((retval = nc_def_var(grpid1, "t2m_index", NC_SHORT, 2, dimids2, &t2m_GPROF_varid))) ERR(retval);
     sprintf(units,"2-m air temperature index from GPROF, in Kelvin");
     if ((retval = nc_put_att_text(grpid1, t2m_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF column vapor index */
     if ((retval = nc_def_var(grpid1, "tqv_index", NC_BYTE, 2, dimids2, &tqv_GPROF_varid))) ERR(retval);
     sprintf(units,"total column water vapor index from GPROF, in mm");
     if ((retval = nc_put_att_text(grpid1, tqv_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF frozen precip */
     if ((retval = nc_def_var(grpid1, "frozen_precip", NC_FLOAT, 2, dimids2, &frozen_precip_GPROF_varid))) ERR(retval);
     sprintf(units,"frozen precip estimate from GPROF, in mm/hr");
     if ((retval = nc_put_att_text(grpid1, frozen_precip_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF surface class */
     if ((retval = nc_def_var(grpid1, "sfc_class", NC_BYTE, 2, dimids2, &sfc_GPROF_varid))) ERR(retval);
     sprintf(units,"TELSEM surface class index");
     if ((retval = nc_put_att_text(grpid1, sfc_GPROF_varid, "units", strlen(units), units))) ERR(retval);
   }

   if ( have_model > 0 ) {
     /* actual MERRA2 Ts */
     if ((retval = nc_def_var(grpid2, "ts", NC_FLOAT, 2, dimids2, &ts_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 surface temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, ts_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 T2m */
     if ((retval = nc_def_var(grpid2, "t2m", NC_FLOAT, 2, dimids2, &t2m_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 2-m air temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, t2m_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 column vapor  */
     if ((retval = nc_def_var(grpid2, "tqv", NC_FLOAT, 2, dimids2, &tqv_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 total column water vapor, in mm");
     if ((retval = nc_put_att_text(grpid2, tqv_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 273K level */
     if ((retval = nc_def_var(grpid2, "h273", NC_SHORT, 2, dimids2, &h273_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 freezing level, in meters");
     if ((retval = nc_put_att_text(grpid2, h273_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

   }



   dimids3[0] = scan_dimid;
   dimids3[1] = pix_dimid;
   dimids3[2] = tb_dimid;

   dimids4[0] = scan_dimid;
   dimids4[1] = pix_dimid;
   dimids4[2] = em_dimid;

   dimids5[0] = scan_dimid;
   dimids5[1] = pix_dimid;
   dimids5[2] = lev_dimid;

   /* TB */
   if ((retval = nc_def_var(ncid, "Tb", NC_FLOAT, 3, dimids3, &tb_varid))) ERR(retval);
   sprintf(units,"Kelvin");
   if ((retval = nc_put_att_text(ncid, tb_varid, "units", strlen(units), units))) ERR(retval);

   /* reconstructed EMIS vector everywhere */
   if ((retval = nc_def_var(ncid, "emis", NC_FLOAT, 3, dimids4, &emis_varid))) ERR(retval);
   sprintf(units,"- - - - - - - - - Kelvin mm");
   if ((retval = nc_put_att_text(ncid, emis_varid, "units", strlen(units), units))) ERR(retval);

   /* EPC vector everywhere */
   if ((retval = nc_def_var(ncid, "pc_emis", NC_FLOAT, 3, dimids4, &pc_emis_varid))) ERR(retval);
   sprintf(units,"principal components of the emis vector");
   if ((retval = nc_put_att_text(ncid, pc_emis_varid, "units", strlen(units), units))) ERR(retval);


   /* top candidate Ku profiles */
   if ((retval = nc_def_var(grpid_NS, "z_ku", NC_SHORT, 3, dimids5, &z_ku_NS_varid))) ERR(retval);
   sprintf(units,"ku-band DPR profile from top-ranked candidate, dB scaled by 100");
   if ((retval = nc_put_att_text(grpid_NS, z_ku_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* top candidate TB */
   if ((retval = nc_def_var(grpid_NS, "Tb1", NC_FLOAT, 3, dimids3, &tb1_NS_varid))) ERR(retval);
   sprintf(units,"TB from top ranked profile, Kelvin");
   if ((retval = nc_put_att_text(grpid_NS, tb1_NS_varid, "units", strlen(units), units))) ERR(retval);



   /* top candidate Ku Ka profiles */
   if ((retval = nc_def_var(grpid_MS, "z_ku", NC_SHORT, 3, dimids5, &z_ku_MS_varid))) ERR(retval);
   sprintf(units,"ku-band DPR profile from top-ranked candidate, dB scaled by 100");
   if ((retval = nc_put_att_text(grpid_MS, z_ku_MS_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(grpid_MS, "z_ka", NC_SHORT, 3, dimids5, &z_ka_MS_varid))) ERR(retval);
   sprintf(units,"ka-band DPR profile from top-ranked candidate, dB scaled by 100");
   if ((retval = nc_put_att_text(grpid_MS, z_ka_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* top candidate TB */
   if ((retval = nc_def_var(grpid_MS, "Tb1", NC_FLOAT, 3, dimids3, &tb1_MS_varid))) ERR(retval);
   sprintf(units,"TB from top ranked profile, Kelvin");
   if ((retval = nc_put_att_text(grpid_MS, tb1_MS_varid, "units", strlen(units), units))) ERR(retval);


   /* on the fly compression for the 3-D vars */
   if ((status= nc_def_var_deflate(ncid, tb_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(ncid, emis_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(ncid, pc_emis_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(grpid_NS, z_ku_NS_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(grpid_NS, tb1_NS_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(grpid_MS, z_ku_MS_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(grpid_MS, z_ka_MS_varid, 0, 1, 4))) ERR(status);
   if ((status= nc_def_var_deflate(grpid_MS, tb1_MS_varid, 0, 1, 4))) ERR(status);


   if ((retval = nc_enddef(ncid))) ERR(retval);


   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     dbuf[n]= secs_S1_buf[isat];
     n++;
   }
   if ((retval = nc_put_var_double(ncid, secs_varid, dbuf))) ERR(retval);

   if ((retval = nc_put_var_short(ncid, nem_compare_varid, em_compare))) ERR(retval);

   count2[0] = nl_out;
   count2[1] = ns_S1;
   start2[0] = 0;
   start2[1] = 0;

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       fbuf[idx1]= lat_S1_buf[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, lat_varid, start2, count2, fbuf))) ERR(retval);

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       fbuf[idx1]= lon_S1_buf[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, lon_varid, start2, count2, fbuf))) ERR(retval);

   /* estimated precip and precip2 for the NS and MS databases */

   if ((retval = nc_put_vara_float(grpid_NS, precip_NS_varid, start2, count2, precip_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip2_NS_varid, start2, count2, precip2_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip_prob_NS_varid, start2, count2, precip_prob_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip2_prob_NS_varid, start2, count2, precip2_prob_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_int(grpid_NS, quality_NS_varid, start2, count2, qual_NS_buf))) ERR(retval);

   if ((retval = nc_put_vara_float(grpid_MS, precip_MS_varid, start2, count2, precip_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip2_MS_varid, start2, count2, precip2_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip_prob_MS_varid, start2, count2, precip_prob_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip2_prob_MS_varid, start2, count2, precip2_prob_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_int(grpid_MS, quality_MS_varid, start2, count2, qual_MS_buf))) ERR(retval);

   /* estimated ts, t2m, tqv, h273, Ztop */

   if ((retval = nc_put_vara_float(grpid_NS, ts_NS_varid, start2, count2, ts_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, t2m_NS_varid, start2, count2, t2m_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, tqv_NS_varid, start2, count2, tqv_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, h273_NS_varid, start2, count2, h273_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, ku_ztop_NS_varid, start2, count2, ku_ztop_epc_NS_buf))) ERR(retval);

   if ((retval = nc_put_vara_float(grpid_MS, ts_MS_varid, start2, count2, ts_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, t2m_MS_varid, start2, count2, t2m_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, tqv_MS_varid, start2, count2, tqv_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, h273_MS_varid, start2, count2, h273_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, ku_ztop_MS_varid, start2, count2, ku_ztop_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, ka_ztop_MS_varid, start2, count2, ka_ztop_epc_MS_buf))) ERR(retval);

   /* GPROF variables if provided */

   if ( have_gprof > 0 ) {
     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= surfacePrecipitation_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, precip_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= probabilityOfPrecip_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, precip_prob_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         sbuf[idx1]= temp2mIndex_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_short(grpid1, t2m_GPROF_varid, start2, count2, sbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         bbuf[idx1]= totalColumnWaterVaporIndex_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_schar(grpid1, tqv_GPROF_varid, start2, count2, bbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= frozenPrecipitation_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, frozen_precip_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         bbuf[idx1]= sfc_GPROF_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_schar(grpid1, sfc_GPROF_varid, start2, count2, bbuf))) ERR(retval);

   }

   /* Actual MERRA2 variables, if provided */

   if ( have_model > 0 ) {
     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= ts_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, ts_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= t2m_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, t2m_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= tqv_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, tqv_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         sbuf[idx1]= h273_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_short(grpid2, h273_MERRA2_varid, start2, count2, sbuf))) ERR(retval);

   }


   /*--- GMI TB and the top-ranked TB database entries for the NS and MS searches */

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = nc_out;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       for (k= 0; k < nc_out; k++) {
         idx= isat*ns_S1*nc_out + jsat*nc_out + k;
         idx1= n*ns_S1*nc_out + jsat*nc_out + k;
         fbuf3[idx1]= tb_S1S2_buf[idx];
       }
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, tb_varid, start3, count3, fbuf3))) ERR(retval);

   if ((retval = nc_put_vara_float(grpid_NS, tb1_NS_varid, start3, count3, tb_1_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, tb1_MS_varid, start3, count3, tb_1_MS_buf))) ERR(retval);

   /*--- emissivity vector and its EPC ---*/

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = NEM;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   if ((retval = nc_put_vara_float(ncid, emis_varid, start3, count3, emis_S1S2_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(ncid, pc_emis_varid, start3, count3, pc_emis_S1S2_buf))) ERR(retval);

   /*--- top-ranked Ku and Ka database entries for the NS and MS searches */

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = NLEV_NS;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   if ((retval = nc_put_vara_short(grpid_NS, z_ku_NS_varid, start3, count3, z_ku_1_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, z_ku_MS_varid, start3, count3, z_ku_1_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, z_ka_MS_varid, start3, count3, z_ka_1_MS_buf))) ERR(retval);



   if ((retval = nc_close(ncid))) ERR(retval);
   printf("Closed %s\n", outfile);

   /*-------------------------------------------------------------------------*/

   free(yyyy_S1_buf);
   free(mm_S1_buf);
   free(dd_S1_buf);
   free(hh_S1_buf);
   free(mn_S1_buf);
   free(ss_S1_buf);

   free(tb_S1_buf);
   free(tb_S2_buf);

   free(lat_S1_buf);
   free(lon_S1_buf);
   free(lat_S2_buf);
   free(lon_S2_buf);

   free(sclat_S1_buf);
   free(sclon_S1_buf);

   free(qual_S1_buf);
   free(qual_S2_buf);

   free(tb_S1S2_buf);
   free(pc_emis_S1S2_buf);
   free(emis_S1S2_buf);

   free(z_ku_1_NS_buf);
   free(z_ku_1_MS_buf);
   free(z_ka_1_MS_buf);
   free(tb_1_NS_buf);
   free(tb_1_MS_buf);

   free(fbuf3);
   free(fbuf4);
   free(fbuf);
   free(sbuf);
   free(bbuf);
   free(dbuf);
   free(ibuf);


   exit(0);
}


