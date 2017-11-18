#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
  SPLICE2 - splice several filterbank files starting on same time stamp together
  for use to analyse data from cloned machines sampling several parts of the 
  band. e.g. multi-WAPPs or multi Breakthrough Listen banks. 
  Initial version (Nov 2002) required contiguous bands. 
  New version     (Jun 2003) requires only high-lo frequency ordering
  Modified        (Jul 2005) to write out to file if -o file option is given
  Modified 		  (Apr 2016) splice2 - write out merged contiguous filterbank using just frch1 + froff
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sigproc.h"
#include "header.h"
FILE *output;
main (int argc, char **argv)
{
  int i=1, j, k, nfiles=0, *numbt, schans=0, *nchan, ngaps = 0;
  long int nbytes, headerbytes;
  FILE *input[32];
  char *block;
  char *zeros;
  double *stamp, *frch1, *froff, *frmhz;
  double hackvalue = 0;
  output=stdout;
  /* print help if necessary */
  if (argc<=1 || help_required(argv[1])) {
    /*splice_help();*/
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }

  /* open up files */
  while (i<argc) {
    if (file_exists(argv[i])) {
      input[nfiles]=open_file(argv[i],"rb");
      nfiles++;
    } else if (strings_equal(argv[i],"-o")) {
      output=open_file(argv[++i],"wb");
    } else if (strings_equal(argv[i],"-h")) {
      hackvalue = atof(argv[++i]);
    }
    i++;
  }


  
  /* read in headers and check time stamps */
  stamp = (double *) malloc(nfiles*sizeof(double));
  frch1 = (double *) malloc(nfiles*sizeof(double));
  froff = (double *) malloc(nfiles*sizeof(double));
  numbt = (int *) malloc(nfiles*sizeof(int));
  nchan = (int *) malloc(nfiles*sizeof(int));
  for (i=0; i<nfiles; i++) {
    headerbytes = read_header(input[i]);
    //fprintf(stderr, "read %d bytes from header of file %s\n",headerbytes, argv[i+1]); 
    if (headerbytes > 0) {
      stamp[i]=tstart;
      frch1[i]=fch1;
      froff[i]=foff;
      numbt[i]=nbits;
      nchan[i]=nchans;
      schans+=nchans;
    } else {
      error_message("problem reading header parameters");
    }
    if (data_type != 1) 
      error_message("input data are not in filterbank format!");
    if (stamp[i] != stamp[0]) 
      error_message("start times in input files are not identical!");
    if (numbt[i] != numbt[0])
      error_message("number of bits per sample in input files not identical!");
    if (i>0) {
      //if (frch1[i] >= frch1[i-1]) error_message("input files not ordered in descending frequency!");

      if (frch1[i] != frch1[i-1] + nchan[i-1] * froff[i]) {
      	printf("mismatch! current start: %f last ending freq: %f\n", frch1[i], frch1[i-1] + nchan[i-1] * froff[i-1]);
      	printf("Will pad with zeros for %i banks...\n", ngaps = ngaps + (int)( ((frch1[i-1] + nchan[i-1] * froff[i-1]) - frch1[i]) / fabsf(froff[i-1] * nchan[i-1])) );
		
		}
    }

  }


  send_string("HEADER_START");
  send_int("machine_id",machine_id);
  send_int("telescope_id",telescope_id);
  send_int("data_type",1);

/*
  send_string("FREQUENCY_START");
  send_int("nchans",schans);
  frmhz = (double *) malloc(sizeof(double)*schans);
  k=0;
  for (i=0; i<nfiles; i++) {
    for (j=0; j<nchans; j++) {
      frmhz[k]=frch1[i]+j*froff[i];
      send_double("fchannel",frmhz[k++]);
    }
  }
  send_string("FREQUENCY_END");
*/ 
 

  /* files are in descending frequency order, so just send ch1 and offset from file 0 */
  send_double("fch1",frch1[0] + hackvalue);
  send_double("foff",froff[0]);
  send_int("nchans",nchans * (nfiles + ngaps));

  if (!strings_equal(source_name,"")) {
    send_string("source_name");
    send_string(source_name);
  }
  send_coords(src_raj,src_dej,az_start,za_start);
  send_int("nbits",nbits);
  send_double("tstart",tstart);
  send_double("tsamp",tsamp);
  send_int("nifs",nifs);
  send_string("HEADER_END");

  nbytes = (long int) nchans* (long int) nbits/8;
  fprintf(stderr, "nbytes per spectra: %ld\n", nbytes);
  block = (char *) malloc(nbytes);
  zeros = (char *) malloc(nbytes);
  memset(zeros, 0x0, nbytes);
  
  while (1) {
    for (i=0; i<nfiles; i++) {
      headerbytes = fread(block,1,nbytes,input[i]);
      if (feof(input[i])) exit(0);
        if (frch1[i] != frch1[i-1] + nchan[i-1] * froff[i]) {
			k = (int)( ((frch1[i-1] + nchan[i-1] * froff[i-1]) - frch1[i]) / fabsf(froff[i-1] * nchan[i-1]));
			for(j = 0;j<k;j++) fwrite(zeros,nbytes,1,output);
		}
      fwrite(block,nbytes,1,output);
    }
  }
}
