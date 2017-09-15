#include <string.h>

/*
  FILTERBANK  - convert raw data from a pulsar machine into "filterbank data"
  a stream of samples: [s1,c1] [s1,c2]...[s1,cn] [s2,c1]...[s2,cn]..... etc
  where the indices s1, s2.... are individual sample readouts and c1, c2.... 
  are the individual channels. For each readout, there are "n" channels.
  Before writing out a stream of samples, filterbank sends a header string
  which is intercepted by other programs (dedisperse, fold etc) downstream.
  (N.B. The header can be omitted using the -headerfile command-line option)
*/

#include "filterbank.h"

main (int argc, char *argv[])
{
  int   i,j,k,nfiles,fileidx,fileidx_start,inputdata,opened=0;
  char  message[80], chan_name[4], fil_value[1];
  FILE  *mult_out[11][4];
  int   mask[11][4];
  int   numsamps,summed=0;
  int   bpp_headersize = 32768;
  float f, fcentral;
  int   sample_end,fileidx_final,byte_end,bytefinal;
  int   data_size,sample_skip,bytestart,ns,blocksize;
  double  sample_final,scantime;
  double  sample_start;
  /* now actually convert the raw data into filterbank format */
  unsigned int ibob=48, channel=2; /* A=3 B=2 C=1 D=0 */
  int          reversed=1;         /* reverse the band */
  /* maybe later change to figure out the 1601 by itself (if 0, then 1600, then remove 2) */
  int          sample_rate=1600;   
  int          skip_packet_2=0;    /* skip packet number 2 for 1601 Hz problem */
  int          start_at_zero=0;    /* wait for counter to reset before start */
  int          check_order=0;      /* check if all iBOBs are arriving in order */
  /* command-line changeable: */
  int          equalize_all=0;     /* equalize channels and telescopes */
  int          sumall=0;           /* sum all ibobs and channels to single output */
  char         file_prepend[120];


  /* ata transient pcap variables */
  char cblock[128], payload[100000];
  unsigned int block[128], equalize_sum[128][11][4];
  float        equalize_coeff[128][11][4], average_power=0; 
  float        temp_power=0.0;
  unsigned int bump, packet_len, packet_saved_len, ip, s;
  unsigned int oldaccumulation[11], accumulation_number[11];
  int          equalize_init=0, first_packet=1;
  int          ts_sec, ts_usec; /* packet timestamp seconds, microcseconds */
  long long    counter[11], prev_counter, last_counter[11];
  
  /* set and clear arrays */
  strcpy(chan_name,"DCBA"); /* A=3 B=2 C=1 D=0 */
  strcpy(file_prepend, "");
  for (i=0; i<11; i++) { accumulation_number[i]=counter[i]=0; oldaccumulation[i]=-1; }
  memset(block,        0, 128*     sizeof(unsigned int));
  memset(equalize_sum, 0, 128*11*4*sizeof(unsigned int));
  memset(mask,         0, 11*4*sizeof(int));
  

  /* check number of command-line arguments */
  if (argc<2) {
    filterbank_help();
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }
 
  /* print help if necessary */
  if (help_required(argv[1])) {
    filterbank_help();
    exit(0);
  }
 
  /* set up default global variables */
  hanning=hamming=zerolagdump=swapout=sumifs=headerless=headerfile=0;
  invert_band=clip_threshold=headeronly=0;
  time_offset=start_time=final_time=tstart=0.0;
  obits=-1;
  do_vanvleck=compute_spectra=1;
  strcpy(ifstream,"XXXX");

  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i])) {
        nfiles++;
        i++;
  }
  if (!nfiles) error_message("no input files supplied on command line!");
  fileidx=1;

  /* now parse any remaining command line parameters */
  if (argc>nfiles) {
    i=nfiles+1;
    while (i<argc) {
      if (strings_equal(argv[i],"-o")) {
        /* get and open file for output */
        strcpy(outfile,argv[++i]);
        output=fopen(outfile,"wb");
        opened=1;
      /* ata filterbank options */
      } else if (strings_equal(argv[i],"-mjd")) {
          /* get the fractional start mjd */
          tstart=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-eq")) {
        /* equalize channels and ibobs */
          equalize_all=1;
      } else if (strings_equal(argv[i],"-sumall")) {
        /* equalize channels and ibobs */
          sumall=1;
      } else if (strings_equal(argv[i],"-skip2")) {
        /* skip spectra 2 for funny 1pps stuff */
        skip_packet_2=1; 
      } else if (strings_equal(argv[i],"-prepend")) {
        strcpy(file_prepend,argv[++i]);
      } else {
        /* unknown argument passed down - stop! */
          filterbank_help();

        sprintf(message,"unknown argument (%s) passed to filterbank.",argv[i]);
        error_message(message);
      }
      i++;
    }
  }




  /* hard code some values for now */
  machine_id=9;    //-1
  telescope_id=9;  //-1
  data_type=1;
  nchans=128;
  foff=838.860800/4/nchans * -1.0;        // XXX hack
  fcentral=1430.0;
  fch1=fcentral  - ((nchans/2)+0.5)*foff;    // A bit of a guess
  nbits=8;
  tsamp=1.0/(double)sample_rate;
  nifs=1;
  strcpy(ifstream,"YYYY");
  obits = 32;
  nbits = 32;

    /* open up input file */
    strcpy(inpfile,argv[1]);
    input=open_file(inpfile,"rb");


    /* broadcast the header */
    /* to the single output file */
    filterbank_header(output);

while(fread(&fil_value, sizeof(char), 1, input)) fwrite(&fil_value, sizeof(char), 1, output);


       

  update_log("finished");
  close_log();
  fclose(output);
  exit(0);
  fprintf(stderr,"Completed Successfully\n");

}

