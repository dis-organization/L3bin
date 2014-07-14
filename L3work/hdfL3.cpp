#include <math.h>
#include "hdf.h"

//#include <Rcpp.h>
//using namespace Rcpp;

#define NUMROWS     2160
#define PI  3.1415926535897932384626433832795029L
#define BLIST_FIELDS "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set"
#define BLIST_SIZE  19

static int    basebin[NUMROWS];
//static int16    numbin[NUMROWS];
static short    numbin[NUMROWS];
//static float64  latbin[NUMROWS];
static double  latbin[NUMROWS];
//static int32    totbins;
static int   totbins;

#define PREC_SIZE   8
static uint8    paramrec[PREC_SIZE];
static char *param_fields;
static float32  sum,sum_sq;
static VOIDP    paramptrs[] = {&sum,&sum_sq};


int binlist(int x);

void    initbin(void);

int main() {
  return binlist(1);
}
//// [[Rcpp::export]]
int binlist(int x) {
 int     numbins = 0, *binnums = NULL;
  int     i;
  int     file_id,vdata_ref, vdata_id,numrecs,pvd_id;
  char PARAM[] = "chlor_a";
  initbin();
  
  /* Open the HDF file. */
  file_id = Hopen("S1998001.L3b_DAY_CHL.main", DFACC_READ, 0);
  if(file_id == FAIL){
    fprintf(stderr,"-E- %s line %d: Hopen(%s,DFACC_READ,0) failed.\n",
      __FILE__,__LINE__,"S1998001.L3b_DAY_CHL.main");
    return(EXIT_FAILURE);
  }
  /* Initialize the Vdata interface. */
  if(Vstart(file_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: Vstart(%d) failed.\n",
	    __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }
  
  /* Open the "BinList" Vdata. */
  vdata_ref = VSfind(file_id,"BinList");
  if(vdata_ref == FAIL){
    fprintf(stderr,"-E- %s line %d: VSfind(%d,\"BinList\") failed.\n",
	    __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }
  
  vdata_id = VSattach(file_id, vdata_ref, "r");
  if(vdata_id == FAIL){
    fprintf(stderr,"-E- %s line %d: VSattach(%d,%d,\"r\") failed.\n",
      __FILE__,__LINE__,file_id,vdata_ref);
    return(EXIT_FAILURE);
  }
  /* Find out how many bins are stored in this file. */
  numrecs = VSelts(vdata_id);
  if(numrecs == FAIL){
    fprintf(stderr,"-E- %s line %d: VSelts(%d) failed.\n",
	    __FILE__,__LINE__,vdata_id);
    return(EXIT_FAILURE);
  }
  
  // this doesn't seem to be necessary?
  /* Set up to read the fields in the BinList Vdata records. */
//  if(VSsetfields(vdata_id,BLIST_FIELDS) == FAIL){
//    fprintf(stderr,"-E- %s line %d: VSsetfields(%d,%s) failed.\n",
//	    __FILE__,__LINE__,vdata_id,BLIST_FIELDS);
//    return(EXIT_FAILURE);
//  }

  /* Open the parameter-specific Vdata. */
  vdata_ref = VSfind(file_id,PARAM);
  if(vdata_ref == 0){
    fprintf(stderr,"-E- %s line %d: VSfind(%d,\"%s\") failed.\n",
	    __FILE__,__LINE__,file_id,PARAM);
    return(EXIT_FAILURE);
  }
  pvd_id = VSattach(file_id, vdata_ref, "r");
  if(pvd_id == FAIL){
    fprintf(stderr,"-E- %s line %d: VSattach(%d,%d,\"r\") failed.\n",
	    __FILE__,__LINE__,file_id,vdata_ref);
    return(EXIT_FAILURE);
  }
    /* Set up to read the fields in the parameter-specific Vdata records. */
  {
    int len;
    len = 2*strlen(PARAM) + strlen("_sum,") + strlen("_sum_sq") + 1;
    param_fields = (char *)malloc(len);
    if(param_fields == NULL){
      fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
        __FILE__,__LINE__);
      return(EXIT_FAILURE);
    }
    strcpy(param_fields,PARAM);
    strcat(param_fields,"_sum,");
    strcat(param_fields,PARAM);
    strcat(param_fields,"_sum_sq");
  }
  if(VSsetfields(pvd_id,param_fields) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSsetfields(%d,%s) failed.\n",
	    __FILE__,__LINE__,pvd_id,param_fields);
    return(EXIT_FAILURE);
  }

  printf(param_fields);
  printf("\n");





  /* Output a header record to identify the fields written out below. */
  printf("%80s%15.15s %15.15s\n"," ",PARAM,PARAM);
  printf("    bin centerlat  centerlon");
  printf("     north     south       west       east");
  printf("    n   N         sum_obs sum_squared_obs          weight");
  printf("  time_trend_bits                     l2_flag_bits sel\n");
  printf("------- --------- ----------");
  printf(" --------- --------- ---------- ----------");
  printf(" ---- --- --------------- --------------- ---------------");
  printf(" ---------------- -------------------------------- ---\n");
 
 // for(i = 0; i < numbins; i++){
  for(i = 0; i < 10; i++){
    int   recno;

   // recno = binsearch(binnums[i],vdata_id,numrecs);
    recno = i;
    if(recno >= 0){
      double    n,s,w,e,clat,clon;

      /*
      Read the sum and sum-of-squares for the
      the specified parameter for this bin.
      */
      if(VSseek(pvd_id,recno) == FAIL){
        fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n",
  	__FILE__,__LINE__,pvd_id,recno);
	return(EXIT_FAILURE);
      }
      if(VSread(pvd_id,paramrec,1,FULL_INTERLACE) != 1){
        fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
	fprintf(stderr,"VSread(%d,paramrec,1,FULL_INTERLACE) failed.\n",
		pvd_id);
        return(EXIT_FAILURE);
      }
      /*
      VSfpack() sets the global sum and sum_sq variables
      via the paramptrs pointer array.
      */
      if(
	 VSfpack(
		 pvd_id,_HDF_VSUNPACK,param_fields,paramrec,PREC_SIZE,1,NULL,paramptrs
		 )
	 == FAIL){
        fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
	fprintf(stderr,"VSfpack(%d, ...) failed.\n", pvd_id);
	return(EXIT_FAILURE);
      }

      /* Get the geographical coordinates associated with this bin. */
    //  bin2latlon(binnums[i],&clat,&clon);
    //  bin2bounds(binnums[i],&n,&s,&w,&e);

      /* Output the results. */
    //  printf("%7d %9.5f %10.5f %9.5f %9.5f %10.5f %10.5f ",
	  //   binnums[i],clat,clon,n,s,w,e);
    //  printf("%4d %3d ",nobs,nscenes);
    //  printf("% .8e % .8e % .8e ",sum,sum_sq,weights);
      printf("% .8e % .8e % .8e ",sum,sum_sq,0.0);
   //   printf("%.16s %.32s ",bitstr16(time_rec),bitstr32(flags_set));
   //   printf("%3d",sel_cat);
      printf("\n");
    }

  }

  if(VSdetach(pvd_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSdetach(%d) failed.\n",
	    __FILE__,__LINE__,pvd_id);
    return(EXIT_FAILURE);
  }
  if(VSdetach(vdata_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSdetach(%d) failed.\n",
	    __FILE__,__LINE__,vdata_id);
    return(EXIT_FAILURE);
  }
  if(Vend(file_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: Vend(%d) failed.\n",
	    __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }
  if(Hclose(file_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: Hclose(%d) failed.\n",
	    __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }

  free(param_fields);
  free(binnums);






   return 2;
}



/*
The following functions are based on the pseudocode found in Appendix A of:

Campbell, J.W., J.M. Blaisdell, and M. Darzi, 1995:
Level-3 SeaWiFS Data Products: Spatial and Temporal Binning Algorithms.
NASA Tech. Memo. 104566, Vol. 32,
S.B. Hooker, E.R. Firestone, and J.G. Acker, Eds.,
NASA Goddard Space Flight Center, Greenbelt, Maryland
*/

void initbin(void){
  int   row;

  basebin[0] = 1;
  for(row=0; row<NUMROWS; row++){
    latbin[row] = ((row + 0.5)*180.0/NUMROWS) - 90.0;
    //numbin[row] = (int16)(2*NUMROWS*cos(latbin[row]*PI/180.0) + 0.5);
    numbin[row] = (short)(2*NUMROWS*cos(latbin[row]*PI/180.0) + 0.5);
    if(row > 0){
      basebin[row] = basebin[row - 1] + numbin[row - 1];
    }
  }
  totbins = basebin[NUMROWS - 1] + numbin[NUMROWS - 1] - 1;
}

//int16 lat2row(double lat){
short lat2row(double lat){
//  int16 row;
  short row;

//  row = (int16)((90 + lat)*NUMROWS/180.0);
  row = (short)((90 + lat)*NUMROWS/180.0);
  if(row >= NUMROWS) row = NUMROWS - 1;
  return(row);
}
//
//int32 rowlon2bin(int16 row, double lon){
//  int16 col;
//  int32 bin;
//
//  lon = constrain_lon(lon);
//  col = (int16)((lon + 180.0)*numbin[row]/360.0);
//  if(col >= numbin[row]) col = numbin[row] - 1;
//  bin = basebin[row] + col;
//  return(bin);
//}
//
//int32 latlon2bin(double lat, double lon){
//  int16 row, col;
//  int32 bin;
//
//  /* Constrain latitudes to [-90,90] and longitudes to [-180,180]. */
//  lat = constrain_lat(lat);
//  lon = constrain_lon(lon);
//
//  row = lat2row(lat);
//  col = (int16)((lon + 180.0)*numbin[row]/360.0);
//  if(col >= numbin[row]) col = numbin[row] - 1;
//  bin = basebin[row] + col;
//  return(bin);
//}
//
//void bin2latlon(int32 bin, double *clat, double *clon){
//  int16 row;
//
//  row = NUMROWS - 1;
//  if(bin < 1) bin = 1;
//  while(bin < basebin[row]) row--;
//  *clat = latbin[row];
//  *clon = 360.0*(bin - basebin[row] + 0.5)/numbin[row] - 180.0;
//}
//
//void bin2bounds(
//		int32   bin,
//		double  *north,
//		double  *south,
//		double  *west,
//double  *east
//		){
//  int16     row;
//  double    lon;
//
//  row = NUMROWS - 1;
//  if(bin < 1) bin = 1;
//  while(bin < basebin[row]) row--;
//  *north = latbin[row] + 90.0/NUMROWS;
//  *south = latbin[row] - 90.0/NUMROWS;
//  lon = 360.0*(bin - basebin[row] + 0.5)/numbin[row] - 180.0;
//  *west = lon - 180.0/numbin[row];
//  *east = lon + 180.0/numbin[row];
//}
//
//double constrain_lat(double lat){
//  if(lat >  90) lat =  90;
//  if(lat < -90) lat = -90;
//  return(lat);
//}
//
//double constrain_lon(double lon){
//  while(lon < -180) lon += 360;
//  while(lon >  180) lon -= 360;
//  return(lon);
//}

//l3bindump.c functions
//
//
//int32 binsearch(int32 bin, int32 vdata_id, int32 numrecs){
//  int32 lo, hi, mid;
//
//  lo = 0;
//  hi = numrecs - 1;
//  while(lo <= hi){
//    mid = (lo + hi)/2;
//    if(VSseek(vdata_id,mid) == FAIL){
//      fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n",
//        __FILE__,__LINE__,vdata_id,mid);
//      exit(EXIT_FAILURE);
//    }
//    if(VSread(vdata_id,blistrec,1,FULL_INTERLACE) != 1){
//      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
//      fprintf(stderr,"VSread(%d,blistrec,1,FULL_INTERLACE) failed.\n",
//        vdata_id);
//      exit(EXIT_FAILURE);
//    }
//    /*
//    VSfpack() sets the global bin_num variable (and others)
//    via the bufptrs pointer array.
//    */
//    if(
//       VSfpack(
//	       vdata_id,_HDF_VSUNPACK,BLIST_FIELDS,blistrec,BLIST_SIZE,1,NULL,bufptrs
//	       )
//       == FAIL){
//      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
//      fprintf(stderr,"VSfpack(%d, ...) failed.\n", vdata_id);
//      exit(EXIT_FAILURE);
//    }
//    if     (bin < bin_num) hi = mid - 1;
//    else if(bin > bin_num) lo = mid + 1;
//    else                   return(mid);
//  }
//  return(-1);
//}
//
//char *bitstr16(int16 n){
//  static char   str[17];
//  int       i;
//
//  str[16]=0;
//  for(i = 0; i < 16; i++){
//    if(n & (1 << i)) str[i] = '1';
//    else             str[i] = '0';
//  }
//  return(str);
//}
//
//char *bitstr32(int32 n){
//  static char   str[33];
//  int       i;
//
//  str[32] = 0;
//  for(i = 0; i < 32; i++){
//    if(n & (1 << i)) str[i] = '1';
//    else             str[i] = '0';
//  } 
//  return(str);
//}