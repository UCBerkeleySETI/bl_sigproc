#include <stdio.h>
#include <string.h>
#include "filterbank.h"

/*
Lspec filterbank processor - read in a lspec binary file and output filterbank format 

A. Siemion
siemion@berkeley.edu
*/


//prototypes for 16bit->8bit and 32bit->8bit quantization functions

unsigned char quantize(double d, double mean, double sigma)
{
	return (unsigned char) (    (    (d - (mean - (4 * sigma))) / (8 * sigma)   ) * 255.0   );
}

//static inline quant16()
//static inline quant32()
void filterbankmlm_help() { 
	printf("usage is...\n");
}


static inline unsigned short bswap_16(unsigned short x);
static inline unsigned int bswap_32(unsigned int x);
static inline unsigned long long bswap_64(unsigned long long x);
unsigned int slice32(unsigned int value,unsigned int width,unsigned int offset);
unsigned int slice64(unsigned long int value,unsigned int width,unsigned int offset);



main (int argc, char *argv[])
{
  FILE *input, *xpowerfile;
  int   i,j,k,l;
  char file[120];

  char  fil_value[1];	
  char payload[100000];
  unsigned int bump, packet_len, packet_saved_len, ip, s;
  long int          ts_sec, ts_usec; /* packet timestamp seconds, microcseconds */

  /* number of clock cycles between accumulations */
  unsigned long int accum = 131072;

  int bytesread = 1;
  unsigned char xpolpower[1024];
  unsigned char ypolpower[1024];
  unsigned short xpolsquared[1024];
  unsigned short ypolsquared[1024];

  unsigned int value=0;
  signed char crossreal=0;
  signed char crossimag=0;  
  int dropped_packet_cnt = 0;
  int speccnt = 0;
  int err=0;
  unsigned long int fields = 0;
  unsigned long int pkt_index=0;
  unsigned long int pkt_index_prev=0;
  unsigned int kurtout;
  unsigned int kurtsigma;
  double xpowersqmean=0;
  double ypowersqmean=0;
  double xpowersqstd=0;
  double ypowersqstd=0;
  float floatval;
  unsigned short int first_packet = 0;
  unsigned int reversed=0;
  unsigned int nfiles;
  unsigned int equalize_all;
  float f, fcentral;
  int        sumall=0;           /* sum pols to single output */
  int 		quant=0;			/*quantize to eight bits */
  char       file_prepend[120];
  char		 message[120];	
  char quantval;
  /* check number of command-line arguments */
  if (argc<3) {
	    fprintf(stderr, "usage is: %s <input.pcap> <output_prefix>\n", argv[0]);
	    return(0);
  } else {

  }
    
    
    
  /* check number of command-line arguments */
  if (argc<2) {
    filterbankmlm_help();
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }
 
  /* print help if necessary */
  if (help_required(argv[1])) {
    filterbankmlm_help();
    exit(0);
  }
 
  /* set up default global variables */
  hanning=hamming=zerolagdump=swapout=sumifs=headerless=headerfile=0;
  invert_band=clip_threshold=headeronly=0;
  time_offset=start_time=final_time=tstart=0.0;
  obits=-1;
  fcentral=0;
  do_vanvleck=compute_spectra=1;
  strcpy(ifstream,"YYYY");

  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i])) {
        nfiles++;
        i++;
  }
  if (!nfiles) error_message("no input files supplied on command line!");

  /* now parse any remaining command line parameters */
  if (argc>nfiles) {
    i=nfiles+1;
    while (i<argc) {
	  if (strings_equal(argv[i],"-mjd")) {
          /* get the fractional start mjd */
          tstart=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-eq")) {
        /* equalize channels and ibobs */
          equalize_all=1;
      } else if (strings_equal(argv[i],"-sumpols")) {
        /* sum polarizations */
          sumall=1;
      } else if (strings_equal(argv[i],"-kurtout")) {
        /* output kurtosis */
          kurtout=1;
      } else if (strings_equal(argv[i],"-kurtreject")) {
        /* apply kurtosis rfi rejection */
          kurtsigma=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-prepend")) {
        strcpy(file_prepend,argv[++i]);
	  } else if (strings_equal(argv[i],"-freq")) {
          /* get the center frequency */
          fcentral=atof(argv[++i]);
	  } else if (strings_equal(argv[i],"-tsamp")) {
          /* get the center frequency */
          tsamp=atof(argv[++i]);
	  } else if (strings_equal(argv[i],"-quant")) {
          /* quantize to eight bits */
          quant=1;
      } else if (strings_equal(argv[i],"-ra")) {
        src_raj=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-dec")) {
        src_dej=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-reverse")) {
        reversed=1;
      } else {
        /* unknown argument passed down - stop! */
          filterbankmlm_help();

        sprintf(message,"unknown argument (%s) passed to filterbank.",argv[i]);
        error_message(message);
      }
      i++;
    }
  }




  /* hard code some values for now */
  machine_id=13;    //machine id 13 leuschner obs lspec
  telescope_id=11;  //telescope id 11 leuschner observatory
  data_type=1;
  nchans=512;
  foff = 400.0/nchans * -1.0;        // XXX hack
  printf("frequency offset is: %f\n", foff);

  if (fcentral==0) {
  	fcentral=205.00;
  	printf("warning... center freq not set, defaulting to 1420 MHz\n");
  }
  fch1=fcentral  - ((nchans/2.0)-0.5)*foff;    // A bit of a guess
 printf("frequency of channel 1 is: %f\n", fch1);
  if (tsamp==0) {
	  tsamp=0.00999936;
      printf("warning... sample time not set, defaulting to 1.25 milliseconds\n");  		
  }
  if(tstart == 0.0) {
  	 tstart = 53000.0;
  	  printf("warning... tstart not set, defaulting to mjd = 53000\n");  		

  } 
  nifs=1;
  obits = 8;
  nbits = 8; 
 printf("tstart: %f\n", tstart);
    
    
      /* open up input file */
    strcpy(file,argv[1]);
    input=fopen(file,"rb");

  
  	/* open up output files */
    
    sprintf(file,"%s_xpower.fil", file_prepend);
    xpowerfile=fopen(file,"wb");


    
    


//	 if(tstart == 0) tstart=40587.0+((double)ts_sec+((double)ts_usec)/1000000.0)/86400.0;

	 first_packet = 1;
	 filterbank_header(xpowerfile);



	while(fread(&fil_value, sizeof(char), 1, input)) fwrite(&fil_value, sizeof(char), 1, xpowerfile);


       

  update_log("finished");
  close_log();


  fprintf(stderr,"Completed Successfully\n");
  fprintf(stderr, "Spectra count: %d\n", speccnt);
  fprintf(stderr, "Dropped packet count: %d\n", dropped_packet_cnt);

  
  fclose(xpowerfile);

  return(0);

}


static inline unsigned short bswap_16(unsigned short x) {
  return (x>>8) | (x<<8);
}

static inline unsigned int bswap_32(unsigned int x) {
  return (bswap_16(x&0xffff)<<16) | (bswap_16(x>>16));
}

static inline unsigned long long bswap_64(unsigned long long x) {
  return (((unsigned long long)bswap_32(x&0xffffffffull))<<32) | (bswap_32(x>>32));
}

/* extract variable bit length values from a 32bit word */
unsigned int slice32(unsigned int value,unsigned int width,unsigned int offset)
{
    value = value << (32 - (width + offset));
    value = value >> (32 - width);
    return value;
}

/* extract variable bit length values (32 bit max!) from a 64bit word */
unsigned int slice64(unsigned long int value,unsigned int width,unsigned int offset)
{
    value = value << (64 - (width + offset));
    value = value >> (64 - width);
    return value;
}

