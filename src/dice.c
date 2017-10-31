#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
  DICE - dice up a filterbank file by removing unwanted channels
  Initial version (Feb 2006) adapted from splice.c
  Modified (Feb 23, 2006) to write out fch1 and foff wherever possible
  Modified (Nov 28, 2006) to forcefully put zeros in bad channels and
	write them when force=1 (currently hardwired). The reason for this
	is so that the data can be sub-banded later on. zapped channels are
	written as zeros.
  Modified (Oct 31, 2017)  Adding "fast" chanstart/chanend mode for Breakthrough Listen
    data.
    
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
FILE *output;

main (int argc, char **argv)
{
  int channum, *numbt, kchans=0, ns, nbytes, *keep, simple; 
  int chanstart, chanend;
  unsigned int bit;
  int force=0;
  double *frmhz;
  long int i=1, j, k;
  chanstart=0;
  chanend=0;

  FILE *input,*keepfile;
  

  float *block;
  unsigned char *charblock;
  unsigned char *onebitblock;
  unsigned short *shortblock;

  char outfile[255];
  char keepfilename[255];
  
   output=stdout;
  
  sprintf(keepfilename, "DUMMUFILE");

    if (argc == 1 || help_required(argv[1])) {
		dice_help();			
		exit(0);
	}   

  if (strings_equal(argv[1],"stdin")) {
	input = stdin;
  } else {
        /* open up file */
 
  	if (!file_exists(argv[1])) {
    		fprintf(stderr, "input file %s does not exist...\n", argv[1]);
    		 exit(0);
    }		
  	input=open_file(argv[1],"rb");
  }

  i++;
	 while (i<argc) {
	 
	 
	       if (help_required(argv[i])) {
				dice_help();			
			   exit(0);
		  
		   } else if (strings_equal(argv[i],"-o")) {

			  /* get and open file for output */
			  strcpy(outfile,argv[++i]);
			  if(file_exists(outfile)) {
				  fprintf(stderr,"output file (%s) exists!\n",argv[i]);
				  exit(1);
			  }
			  output=fopen(outfile,"wb");

		  } else if (strings_equal(argv[i],"-keepfile")) {
			 /* start channel */
			 strcpy(keepfilename,argv[++i]);

		  } else if (strings_equal(argv[i],"-chanstart")) {
			 /* start channel */
			 chanstart=atoi(argv[++i]);
			 //if(chanend == 0) chanend = -1;

		  } else if (strings_equal(argv[i],"-chanend")) {
			 /* end channel */
			 chanend=atoi(argv[++i]);
			 //if(chanstart == 0) chanstart = -1;

		  } else if (strings_equal(argv[i],"--force")) {
			 force = 1;

		  } else {

			dice_help();			
			/* unknown argument passed down - stop! */
			fprintf(stderr,"unknown argument (%s) passed to filterbank.",argv[i]);
			exit(1);
		}
		i++;
	  }

if(argc < 3) {

			dice_help();			
			/* unknown argument passed down - stop! */
			fprintf(stderr,"Must specify input and output.\n");
			exit(1);

}

if(output == NULL) {
			dice_help();			
			/* unknown argument passed down - stop! */
			fprintf(stderr,"Output file doesn't exist.\n");
			exit(1);
}



  /* read in header info */
  if (!read_header(input)) 
    error_message("problem reading header parameters\n");

  if (data_type != 1) 
    error_message("input data are not in filterbank format!\n");

if (!file_exists(keepfilename)) {
  if (chanstart < 0) 
    error_message("chanstart can't be negative!\n");

  if (chanstart > nchans || chanend > nchans) {
    fprintf(stderr, "chanstart: %d  chanend: %d   nchans: %d\n",chanstart, chanend, nchans);   
    error_message("try again!\n");
 }

  if (chanend == 0) { 
    printf("defaulting to chanend = nchan = %d\n", nchans);
  	chanend = nchans;
  }  
}
	  

  /* open up file to give channel numbers to keep */
  if (!file_exists(keepfilename)) {
	k = 0;

    frmhz = (double *) malloc(sizeof(double)*1);
 	
	frmhz[0] = fch1+(chanstart)*foff;
	fch1 = frmhz[0];
	kchans = chanend - chanstart;
	simple=1;

  } else {
	keepfile=open_file(keepfilename,"r");

    keep=(int *) malloc(nchans*sizeof(int));
	frmhz = (double *) malloc(sizeof(double)*nchans);
	for (i=0;i<nchans;i++) keep[i]=0;


	k=0;
 
	while (fscanf(keepfile,"%d\n",&channum)==1) {
	  keep[channum-1]=1;
	  frmhz[k++]=fch1+(channum-1)*foff;
	  kchans++;
    }
  

	/* do a check to see whether fch1 and foff are sufficient to describe
	   the diced channels */
	fch1=frmhz[0];
	simple=1;
	for (i=0;i<kchans;i++) if (frmhz[i] != fch1+i*foff) simple=0;

	/* force the output to include the zapped channels */
	if (force) {
	  kchans=nchans;
	  simple=1;
	}

  }
  



  /* broadcast header for output file */
  send_string("HEADER_START");
  send_int("machine_id",machine_id);
  send_int("telescope_id",telescope_id);
  send_int("data_type",1);
  if (simple) {
    send_double("fch1",fch1);
    send_double("foff",foff);
    send_int("nchans",kchans);
  } else {
    send_string("FREQUENCY_START");
    send_int("nchans",kchans);
    for (i=0; i<kchans; i++) send_double("fchannel",frmhz[i]);
    send_string("FREQUENCY_END");
  }


	  
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






if (!file_exists(keepfilename)) {
  /* Process with chanstart/chanend options */

	if(nbits != 32 && nbits != 8) error_message("dice (chanstart/chanend) - currently only works with 1, 8, 16 bit data");
	
	printf("%d\n",(int) (chanend - chanstart)*(nbits/8));
	
    charblock = (unsigned char *) malloc((chanend - chanstart)*(nbits/8));
    i=0;
	while (1) {
		 fseek(input, (nbits/8) * chanstart, SEEK_CUR);
		 if(fread(charblock, 1, (nbits/8) * (chanend - chanstart), input) != 0) {
		 	i++;
		 } else {
		 	break;
		 }
		 fwrite(charblock, 1, (nbits/8) * (chanend - chanstart), output);
		 fseek(input, (nbits/8) * (nchans - chanend), SEEK_CUR);
		 printf("%ld\n", ftell(input));
	}

  
} else {

  /* Process with keepfile */
  
  block = (float *) malloc(nchans*sizeof(float));


  switch (nbits) {
	 case 1:
	   onebitblock = (unsigned char *) malloc(nchans/8); //* sizeof(unsigned char));
	   break;

	 case 8:
	   charblock = (unsigned char *) malloc(nchans*sizeof(unsigned char));
	   break;

	 case 16:
	   shortblock = (unsigned short *) malloc(nchans*sizeof(unsigned short));
	   break;

  }  


  while (read_block(input,nbits,block,nchans)==nchans) {
    k=0;
    if (force) {
      for (i=0; i<nchans; i++) {
	  if (keep[i]) 
		block[k++]=block[i];
	  else
	  block[k++]=0;
      }	
    } else {
      for (i=0; i<nchans; i++) if (keep[i]) block[k++]=block[i];
    }
    
    switch (nbits) {
    case 1:
      for (i=0;i<k/8;i++){
        onebitblock[i]=0;
		for(j=0;j<8;j++){
          if (block[8*i+j]>0){
            bit = 1<<j;
            onebitblock[i] = onebitblock[i]|bit;
          }
        }
      }
      fwrite(onebitblock,1,k/8,output);
      break;
    case 8:
      for (i=0;i<k;i++) charblock[i]=block[i];
      fwrite(charblock,1,k,output);
      break;
    case 16:
      for (i=0;i<k;i++) shortblock[i]=block[i];
      fwrite(shortblock,2,k,output);
      break;

    default:
    error_message("dice (keepfile) - currently only works with 1, 8, 16 bit data");
    }
  }
  
  
  
}





}
