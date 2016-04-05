#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
   kurt_fil: compute kurtosis filterbank for power and power squared spectral time series.   
   usage:
   kurt_fil <power.fil> <powersq.fil> -o kurtosis.fil
 */
 
#include <stdio.h>
#include <stdlib.h>
#include "header.h"
#include "filterbank.h"

#include <string.h>
#include <gsl/gsl_histogram.h>

int wapp_header_size, wapp_incfile_length;
int nbins;
	double period;

unsigned char quantize(float d, float min, float max)
{
    if(d > max) d = max;
    if(d < min) d = min;
    return (unsigned char)(((d - min) / (max-min)) * 255.0);
}


unsigned char quantizev1(double d, double mean, double sigma)
{
	return (unsigned char) (    (    (d - (mean - (3 * sigma))) / (6 * sigma)   ) * 255.0   );
}

main(int argc, char *argv[]) 
{
	FILE *fileptr[50];  //max input files 50
	FILE *output;
    char buffer[8192];
    double spectra_sum[8192];  //max channels 8192
    
    /* for floats */
    float fbuffer[8192];  //max channels 8192
    
    float fspectra_power[8192];
    float fspectra_powersq[8192];
    float fspectra_kurt[8192];
    double dspectra_power_temp[8192];
    double dspectra_powersq_temp[8192];
    
    float scale_factor;
	float max, min;  /* minimum and maximum values for the quantization region */
	
	float quantmin;
	float quantmax;
	
	double lower, upper;
	
	long int decimate = 128;

	int chanstart = 0; /* channel to start at for quantization level computation */
	int chanend = 0;   /* channel to end at for quantization level computation */


	double nelements = 0;
	double powersum = 0.0;
    double discard = 0.000003;

	char filename[80],*telescope,*backend,*datatype,message[80],unit[16], outfile[120];
	int i,j, k, year,month,day,check,rah,ram,ded,dem, nfiles, opened=0;
	long int m,n;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;
    long long mindatasize=0;
	long long qlen=0;
	int writeobsdbline;

	/* get these from first file on command line, apply to output file */
	int output_nchans=0, output_nbits=0, output_nifs=0, input_nbits;
	long long output_datasize=0;
    long long output_headersize=0;

	long long input_datasize=0;
    long long input_headersize=0;
    
    int output_ptr;
	unsigned char quantval;
    
    float nacc;
    double mean=0.0, sigma=0.0;
    
    
    /* Create histograms and set ranges   */
	/* _full for all pulse detections     */
	/* _clean for rfi rejected detections */
	gsl_histogram *spectra_quant = gsl_histogram_alloc (11221);
	gsl_histogram_set_ranges_uniform (spectra_quant, 0, 11220);

    
  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i]) && (i < argc)) {
        //printf("opening %s\n", argv[i]);
  	    nfiles++;
        i++;
  }
  if (!nfiles) {
	  error_message("no input files supplied on command line!");
	  exit(1);
  }

  if (nfiles != 2) {
      error_message("must have only two input files - power and power squared!");
	  exit(1);
  }


  for(i = 0; i<nfiles; i++) {
			  fileptr[i]=open_file(argv[i+1],"rb");
  }


  if (argc>nfiles) {
	 i=nfiles+1;
	 while (i<argc) {
		if (strings_equal(argv[i],"-o")) {
		
		  /* get and open file for output */
		  strcpy(outfile,argv[++i]);
		  if(file_exists(outfile)) {
			  sprintf(message,"output file (%s) exists!",argv[i]);
			  error_message(message);
			  exit(1);
		  }
		  output=fopen(outfile,"wb");
		  opened=1;
	    } else if (strings_equal(argv[i],"-obits")) {
		   /* number of output bits */
		   output_nbits=atoi(argv[++i]);
	    } else if (strings_equal(argv[i],"-qlen")) {
		   /* number of spectra to use for quantization */
		   qlen=atoi(argv[++i]);
	    } else if (strings_equal(argv[i],"-chanstart")) {
		   /* number of spectra to use for quantization */
		   chanstart=atoi(argv[++i]);
	    } else if (strings_equal(argv[i],"-chanend")) {
		   /* number of spectra to use for quantization */
		   chanend=atoi(argv[++i]);
	    } else if (strings_equal(argv[i],"-nacc")) {
		   /* number of spectra to use for quantization */
		   nacc = atoi(argv[++i]);
		} else {

			/* unknown argument passed down - stop! */
			sprintf(message,"unknown argument (%s) passed to %s.",argv[i], argv[0]);
			error_message(message);
			exit(1);
		}
		i++;
	  }
  }

  if (!opened) {
	  error_message("must have an output file (-o <output>)!");
	  exit(1);
  }

 // if (nacc < 1) {
 //	  error_message("nacc must be specifed!");
 //	  exit(1);
 // }

	pulsarcentric=barycentric=0;
	writeobsdbline=0;
	


    
  for(i = 0; i<nfiles; i++) {
		headersize=read_header(fileptr[i]);		
	    rewind(fileptr[i]);
		datasize=sizeof_file(argv[i+1])-headersize;
	    numsamps=nsamples(argv[i+1],headersize,nbits,nifs,nchans);
		if(output_datasize == 0 || datasize < output_datasize) {
			output_datasize = datasize;
		    output_headersize = headersize;
			output_ptr = i;
		}
		if(output_nchans == 0) {
		  output_nchans = nchans; 
		  input_nbits = nbits; 
		  output_nifs = nifs;		
		}
        if(nchans != output_nchans || nbits != input_nbits || nifs != output_nifs) {
			sprintf(message,"channel/if/bit mismatch, exiting...");
			error_message(message);
			exit(1);
        }
        
        printf("outbits: %d inputbits: %d file: %s headersize: %d nbits: %d nifs: %d nchans: %d\n", output_nbits, input_nbits, argv[i+1],headersize, nbits, nifs, nchans); 
  }


  if (output_nbits == 8 && input_nbits == 8) {
		  
		printf("8bit to 8bit not yet implemented.. exiting...\n");
		exit(1);

	} else if (output_nbits == 32 && input_nbits == 32) {

			
		  
		   printf("minimum data size is: %ld\n", output_datasize);
		   
		   printf("will dump: %d\n", (output_datasize / (long long) (output_nifs * output_nchans)));
		   
		   headersize=read_header(fileptr[output_ptr]);
		   rewind(fileptr[output_ptr]);
		 
		   
		   nacc = ((float) round((tsamp * fabs(foff) * 1000000))) * decimate;
		   printf("tsamp: %f foff: %f nacc: %f\n", tsamp, foff, nacc); 
		 
		 
		   printf("header size lead file: %d\n", headersize); 
		   fread(buffer, sizeof(char), headersize, fileptr[output_ptr]);
		   fwrite(buffer, sizeof(char), headersize, output);
		   rewind(fileptr[output_ptr]);
		 
		 
		 
		  /* bump past header for all input files */
			for(i = 0; i<nfiles; i++) {
				 headersize=read_header(fileptr[i]);	
			 }
		   
		 
		 
		  /* rewind and bump past header for all input files */
			for(i = 0; i<nfiles; i++) {
				 rewind( fileptr[i] );
				 headersize=read_header(fileptr[i]);	
			 }
		 
		 
		 /* outer loop, read through all spectra, quantize to 8 bits, write to file */
		 
	  	for (j = 0; j < ((output_datasize/4) / (long long) (output_nifs * output_nchans)); j = j + decimate){


		
			 for(i=0;i<nchans;i++) fspectra_power[i] = 0.0;			
			 for(i=0;i<nchans;i++) fspectra_powersq[i] = 0.0;			
			 for(i=0;i<nchans;i++) dspectra_power_temp[i] = 0.0;			
			 for(i=0;i<nchans;i++) dspectra_powersq_temp[i] = 0.0;			

			 for(m = 0; m < decimate; m++) {	 			
					 /* read one power spectra */			
				   fread(&fspectra_power, sizeof(float), (output_nifs * output_nchans), fileptr[0]);
		 		   for(i=0;i<nchans;i++) dspectra_power_temp[i] = dspectra_power_temp[i] + (double) fspectra_power[i];
		 		   		 		  
					 /* read one powersq spectra */
				   fread(&fspectra_powersq, sizeof(float), (output_nifs * output_nchans), fileptr[1]);			  			 		
			 	   for(i=0;i<nchans;i++) dspectra_powersq_temp[i] = dspectra_powersq_temp[i] + (double) fspectra_powersq[i];
			 
			 }
			
			 for(i=0;i<nchans;i++) fspectra_power[i] = (float) dspectra_power_temp[i];			
			 for(i=0;i<nchans;i++) fspectra_powersq[i] = (float) dspectra_powersq_temp[i];	
				
			
				  for(k=0; k < nchans; k++) {
					 fspectra_kurt[k] = ((nacc + 1)/(nacc-1)) * ( (nacc * fspectra_powersq[k] /powf(fspectra_power[k],2)) - 1);			 
					 //printf("%f \n",fspectra_kurt[k]);
					 //usleep(500000);
				 }
				  fwrite(fspectra_kurt, sizeof(float), nchans, output);	 				 
			 
		
			//for(k=0; k < nchans; k++) {		
				  //quantval = quantize(fspectra_sum[k], quantmin, quantmax);
				  //fwrite(&quantval, sizeof(char), 1, output);	 				 
			 //}
			 			   		
		}

		} else if (output_nbits == 32 && input_nbits == 8) {

		printf("32bit to 8bit not yet implemented.. exiting...\n");
		exit(1);



		} else if (output_nbits == 8 && input_nbits == 32) {
		
		printf("8bit to 32bit not yet implemented.. exiting...\n");
		
		exit(1);

  	
  	} else {

	  error_message("This routine only works on 8 bit or 32 bit filterbank data!");
	  exit(1);
  	
  	}



}
