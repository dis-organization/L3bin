#include <math.h>
#include "hdf.h"

#include <Rcpp.h>
using namespace Rcpp;

#define NUMROWS     2160
#define BLIST_FIELDS "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set"
#define BLIST_SIZE  19
static uint8    blistrec[BLIST_SIZE];
static int    bin_num;
static int flags_set;
static short    nobs,nscenes,time_rec;
static float32  weights;
static uint8    sel_cat;
static VOIDP    bufptrs[] = {
  &bin_num,&nobs,&nscenes,&time_rec,&weights,&sel_cat,&flags_set
};


static int    basebin[NUMROWS];
static short    numbin[NUMROWS];
static double  latbin[NUMROWS];
static int   totbins;

#define PREC_SIZE   8
static uint8    paramrec[PREC_SIZE];
static char *param_fields;
static float32  summ,sum_sq;
static VOIDP    paramptrs[] = {&summ,&sum_sq};

void    initbin(void);
int   binsearch(int bin, int vdata_id, int numrecs);

//' Say my name. 
//'
//' Short message from sponsor. 
//' @export
// [[Rcpp::export]]
List rcpp_hello_world() {
  CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
  NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  List z            = List::create( x, y, PI) ;
  return z ;
}


//' Basic L3 bin files 
//'
//' Read from L3 bin. 
//' @param filename path to L3 bin OC file (HDF4)
//' @export
// [[Rcpp::export]]
List binlist(CharacterVector filename) {
  //"/home/mdsumner/Git/L3bin/L3work/S1998001.L3b_DAY_CHL.main";
  std::string fname = Rcpp::as<std::string>(filename);
  
  initbin();
  
  const char * c_filename = fname.c_str();
  
  /* Open the HDF file. */
  int     numbins = 0, *binnums = NULL;
  int     i;
  int file_id;
  int     vdata_ref, vdata_id,numrecs,pvd_id;
  char PARAM[] = "chlor_a";
  
  
  /* The entire globe has been selected. */
  int b; 
  binnums = (int *)realloc(binnums,totbins*sizeof(int));
  if(binnums == NULL){
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    return(EXIT_FAILURE);
  }
  for(b = 1; b <= totbins; b++){ binnums[numbins++] = b;
  }
  
  
  
  file_id = Hopen(c_filename, DFACC_READ, 0);
  if(file_id == FAIL){
    fprintf(stderr,"-E- %s line %d: Hopen(%s,DFACC_READ,0) failed.\n", __FILE__, __LINE__, c_filename);
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
  
  NumericVector parsum(numrecs);
  NumericVector parssq(numrecs);
  
  // this doesn't seem to be necessary?
  /* Set up to read the fields in the BinList Vdata records. */
  if(VSsetfields(vdata_id,BLIST_FIELDS) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSsetfields(%d,%s) failed.\n",
	    __FILE__,__LINE__,vdata_id,BLIST_FIELDS);
    return(EXIT_FAILURE);
  }

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


  int   recno;
  int parcount = -1;
for(i = 0; i < numbins; i++){
//for(i = 0; i < 10000; i++){

  
  recno = binsearch(binnums[i],vdata_id,numrecs);
  //recno = i;
  if(recno >= 0){
    parcount++;
    //      double    n,s,w,e,clat,clon;
    
    /*
    Read the sum and sum-of-squares for the
    the specified parameter for this bin.
    */
    
    if(VSseek(pvd_id,recno) == FAIL){
      fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n", __FILE__,__LINE__,pvd_id,recno);
      
      return(EXIT_FAILURE);
    }
    if(VSread(pvd_id,paramrec,1,FULL_INTERLACE) != 1){
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"VSread(%d,paramrec,1,FULL_INTERLACE) failed.\n", pvd_id);
      return(EXIT_FAILURE);
    }
    /*
    VSfpack() sets the global sum and sum_sq variables
    via the paramptrs pointer array.
    */
    if(VSfpack(pvd_id,_HDF_VSUNPACK,param_fields,paramrec,PREC_SIZE,1,NULL,paramptrs)	 == FAIL){
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"VSfpack(%d, ...) failed.\n", pvd_id);
      return(EXIT_FAILURE);
    }
    
     //NumericVector yyy   = NumericVector::create(vdata_id, numrecs, summ, sum_sq, recno) ;
     //return yyy;
//
//      /* Get the geographical coordinates associated with this bin. */
//    //  bin2latlon(binnums[i],&clat,&clon);
//    //  bin2bounds(binnums[i],&n,&s,&w,&e);
//
//      /* Output the results. */
//    //  printf("%7d %9.5f %10.5f %9.5f %9.5f %10.5f %10.5f ",
//	  //   binnums[i],clat,clon,n,s,w,e);

     parsum[parcount] = summ;
     parssq[parcount] = sum_sq;
//    //  printf("%4d %3d ",nobs,nscenes);
//    //  printf("% .8e % .8e % .8e ",sum,sum_sq,weights);
//      printf("% .8e % .8e % .8e ",sum,sum_sq,0.0);
//   //   printf("%.16s %.32s ",bitstr16(time_rec),bitstr32(flags_set));
//   //   printf("%3d",sel_cat);
//      printf("\n");
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



  List z            = List::create(parsum, parssq) ;
  
  //NumericVector y   = NumericVector::create( 0.0, 1.0) ;
  return z ;
}




 int binsearch(int bin, int vdata_id, int numrecs){
  int lo, hi, mid;

  lo = 0;
  hi = numrecs - 1;
  while(lo <= hi){
    mid = (lo + hi)/2;
    if(VSseek(vdata_id,mid) == FAIL){
      fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n",
        __FILE__,__LINE__,vdata_id,mid);
      exit(EXIT_FAILURE);
    }
    if(VSread(vdata_id,blistrec,1,FULL_INTERLACE) != 1){
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"VSread(%d,blistrec,1,FULL_INTERLACE) failed.\n",
	      vdata_id);
      exit(EXIT_FAILURE);
    }
    /*
    VSfpack() sets the global bin_num variable (and others)
    via the bufptrs pointer array.
    */
    if(
       VSfpack(
	       vdata_id,_HDF_VSUNPACK,BLIST_FIELDS,blistrec,BLIST_SIZE,1,NULL,bufptrs
	       )
       == FAIL){
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"VSfpack(%d, ...) failed.\n", vdata_id);
      exit(EXIT_FAILURE);
    }
    if     (bin < bin_num) hi = mid - 1;
    else if(bin > bin_num) lo = mid + 1;
    else                   return(mid);
  }
  return(-1);
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


