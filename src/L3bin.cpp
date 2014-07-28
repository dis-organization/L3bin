// source taken with permission from Norman Kuring (NASA)
// http://oceancolor.gsfc.nasa.gov/DOCS/OCSSW/l3bindump_8c_source.html
// and from smigen.cpp (NASA, Seadas)
// http://oceancolor.gsfc.nasa.gov/DOCS/OCSSW/smigen_8cpp_source.html

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

//' Parameters from initbin. 
//'
//' Initbin. 
//' @export
// [[Rcpp::export]]
List initlist() {
  initbin();
  List alist = List::create(); 
  alist["totbins"] = totbins; 
  alist["NUMROWS"] = NUMROWS;
//  alist["basebin"] = basebin;
//  alist["latbin"] = latbin;
  return alist;
}

//' Longitude and latitude from bin number.  
//'
//' Generate longitude and latitude coordinates from bin number. 
//' @export
// [[Rcpp::export]]
List bin2lonlat(IntegerVector bins){
  short row;
  initbin();
  
  int bin, ibin;
  int n = bins.size();
  NumericVector clon(n);
  NumericVector clat(n);
  
  for (ibin = 0; ibin < n; ibin++) {
    row = NUMROWS - 1;
    bin = bins[ibin];
    while(bin < basebin[row]) row--;
    clat[ibin] = latbin[row];
    clon[ibin] = 360.0*(bin - basebin[row] + 0.5)/numbin[row] - 180.0;
  }
  List lonlat = List::create();
  lonlat["x"] = clon;
  lonlat["y"] = clat;
  return lonlat;
}


//' Basic L3 bin files 
//'
//' Read from L3 bin. 
//' @param filename path to L3 bin OC file (HDF4)
//' @export
// [[Rcpp::export]]
List binlist(CharacterVector filename, CharacterVector vname) {
  // file name
  std::string fname = Rcpp::as<std::string>(filename);
  const char * c_filename = fname.c_str();
  
  // intialize the bin structure
  initbin();
  
  /* counters */
  int b, i, jvar, nvar, numrecs, numbins = 0, *binnums = NULL;
  /* file refs */
  int file_id, vdata_ref, vdata_id,pvd_id;
  
  /* The entire globe */
  binnums = (int *)realloc(binnums,totbins*sizeof(int));
  if(binnums == NULL){
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n", __FILE__,__LINE__);
    return(EXIT_FAILURE);
  }
  for(b = 1; b <= totbins; b++) { 
    binnums[numbins++] = b;
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
  
  /* Rcpp structures for shared binlist */
  IntegerVector bin(numrecs);
  IntegerVector nobservations(numrecs);
  IntegerVector nsc(numrecs);
  NumericVector wghts(numrecs);
  // NumericVector time(numrecs);
  
  /* Set up to read the fields in the BinList Vdata records. */
  if(VSsetfields(vdata_id,BLIST_FIELDS) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSsetfields(%d,%s) failed.\n", __FILE__,__LINE__,vdata_id,BLIST_FIELDS);
    return(EXIT_FAILURE);
  }
  
  // separate loops for vdatat_id and pvd_id, because they seek to different places
  for (i = 0; i < numrecs; i++ ) {
    
    if(VSseek(vdata_id,i) == FAIL){
      fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n", __FILE__,__LINE__,vdata_id,i);
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
    if(VSfpack(vdata_id,_HDF_VSUNPACK,BLIST_FIELDS,blistrec,BLIST_SIZE,1,NULL,bufptrs)== FAIL) {
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"VSfpack(%d, ...) failed.\n", vdata_id);
      exit(EXIT_FAILURE);
    }
    
    // populate the Rcpp vectors
    bin[i] = bin_num;
    nobservations[i] = nobs;
    nsc[i] = nscenes;
    wghts[i] = weights;
    // time[i] = time_rec;   
  }   
  
  // collect Rcpp vectors in a list
  List z  = List::create() ;
  z["bin_num"] = bin;
  z["nobs"] = nobservations;
  z["nscenes"] = nsc;
  z["weights"] = wghts;
  
  // loop over all parameters, get the sum/ssq for each
  nvar = vname.size();
  for (jvar = 0; jvar < nvar; jvar++) {
    // reset the sum/ssq Rcpp vectors each time
    NumericVector parsum(numrecs);
    NumericVector parssq(numrecs);
    
    std::string vstrname  = Rcpp::as<std::string>(vname[jvar]);
    const char * PARAM = vstrname.c_str();
    
    /* Open the parameter-specific Vdata. */
    vdata_ref = VSfind(file_id,PARAM);
    if(vdata_ref == 0){
      fprintf(stderr,"-E- %s line %d: VSfind(%d,\"%s\") failed.\n", __FILE__,__LINE__,file_id,PARAM);
      return(EXIT_FAILURE);
    }
    pvd_id = VSattach(file_id, vdata_ref, "r");
    if(pvd_id == FAIL){
      fprintf(stderr,"-E- %s line %d: VSattach(%d,%d,\"r\") failed.\n", __FILE__,__LINE__,file_id,vdata_ref);
      return(EXIT_FAILURE);
    }
    /* Set up to read the fields in the parameter-specific Vdata records. */
    {
      int len;
      len = 2*strlen(PARAM) + strlen("_sum,") + strlen("_sum_sq") + 1;
      param_fields = (char *)malloc(len);
      if(param_fields == NULL){
        fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n", __FILE__,__LINE__);
        return(EXIT_FAILURE);
      }
      strcpy(param_fields,PARAM);
      strcat(param_fields,"_sum,");
      strcat(param_fields,PARAM);
      strcat(param_fields,"_sum_sq");
    }
    
    if(VSsetfields(pvd_id,param_fields) == FAIL){
      fprintf(stderr,"-E- %s line %d: VSsetfields(%d,%s) failed.\n", __FILE__,__LINE__,pvd_id,param_fields);
      return(EXIT_FAILURE);
    }
    
    for (i = 0; i < numrecs; i++ ) {     
      /*
      Read the sum and sum-of-squares for the
      the specified parameter for this bin.
      */
      
      if(VSseek(pvd_id,i) == FAIL){
        fprintf(stderr,"-E- %s line %d: VSseek(%d,%d) failed.\n", __FILE__,__LINE__,pvd_id,i);        
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
      if(VSfpack(pvd_id,_HDF_VSUNPACK,param_fields,paramrec,PREC_SIZE,1,NULL,paramptrs)   == FAIL){
        fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
        fprintf(stderr,"VSfpack(%d, ...) failed.\n", pvd_id);
        return(EXIT_FAILURE);
      }
      
      parsum[i] = summ;
      parssq[i] = sum_sq;
      
    }
    // sub-list for each parameter (could be done better)
    
    List z2  = List::create() ;
    z2["sum"] = parsum;
    z2["ssq"] = parssq; 
    z[vstrname] = z2;
    
    if(VSdetach(pvd_id) == FAIL){
      fprintf(stderr,"-E- %s line %d: VSdetach(%d) failed.\n", __FILE__,__LINE__,pvd_id);
      return(EXIT_FAILURE);
    }
    
  }
  
  // finish up, close the file
  if(VSdetach(vdata_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: VSdetach(%d) failed.\n", __FILE__,__LINE__,vdata_id);
    return(EXIT_FAILURE);
  }
  if(Vend(file_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: Vend(%d) failed.\n", __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }
  if(Hclose(file_id) == FAIL){
    fprintf(stderr,"-E- %s line %d: Hclose(%d) failed.\n", __FILE__,__LINE__,file_id);
    return(EXIT_FAILURE);
  }
  
  free(param_fields);
  free(binnums);  
  return z;
}


char *bitstr16(int16 n){
  static char   str[17];
  int       i;
  
  str[16]=0;
  for(i = 0; i < 16; i++){
    if(n & (1 << i)) str[i] = '1';
    else             str[i] = '0';
  }
  return(str);
}

char *bitstr32(int32 n){
  static char   str[33];
  int       i;
  
  str[32] = 0;
  for(i = 0; i < 32; i++){
    if(n & (1 << i)) str[i] = '1';
    else             str[i] = '0';
  } 
  return(str);
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

