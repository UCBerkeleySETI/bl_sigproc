#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
   sumfil: sum filterbank format power spectra time series.   
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "filterbank.h"
//#include "sigproc.h"
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
    double spectra_sum[4096];  //max channels 4096
    
    /* for floats */
    float *fbuffer;  
    float *fspectra_sum; 



    float scale_factor;
	float max, min;  /* minimum and maximum values for the quantization region */
	
	float quantmin;
	float quantmax;
	
	double lower, upper;

	int tcollapse = 1;
	int fcollapse = 1;

	int chanstart = 0; /* channel to start at for quantization level computation */
	int chanend = 0;   /* channel to end at for quantization level computation */


	double nelements = 0;
	double powersum = 0.0;
    double discard = 0.000003;

	char filename[80],*telescope,*backend,*datatype,message[80],unit[16], outfile[120];
	long int i, j, k, m, n;
	int year,month,day,check,rah,ram,ded,dem, nfiles, opened=0;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;
    long long mindatasize=0;
	long long qlen=0;
	int writeobsdbline;

	/* get these from first file on command line, apply to output file */
	int output_nchans=0, output_nbits=0, output_nifs=0, input_nbits;
	int input_nchans=0, input_nifs=0;
	
	long long output_datasize=0;
    long long output_headersize=0;

	long long input_datasize=0;
    long long input_headersize=0;
    
    int output_ptr;
	unsigned char quantval;
    
    double mean=0.0, sigma=0.0;
    
    
    /* Create histograms and set ranges   */
	/* _full for all pulse detections     */
	/* _clean for rfi rejected detections */
	gsl_histogram *spectra_quant; //= gsl_histogram_alloc (11221);
	//gsl_histogram_set_ranges_uniform (spectra_quant, 0, 11220);

    
  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i]) && (i < argc)) {
        //printf("opening %s\n", argv[i]);
  	    nfiles++;
        i++;
  }



  for(i = 0; i<nfiles; i++) {
			  fileptr[i]=open_file(argv[i+1],"rb");
  }


  if (argc>nfiles) {
	 i=nfiles+1;
	 while (i<argc) {
	 
	 
	       if (help_required(argv[i])) {
			   sum_fil_help();
			   exit(0);
		   } else if (strings_equal(argv[i],"-o")) {
		
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
			 /* start channel */
			 chanstart=atoi(argv[++i]);

		  } else if (strings_equal(argv[i],"-chanend")) {
			 /* end channel */
			 chanend=atoi(argv[++i]);

		  } else if (strings_equal(argv[i],"-tcollapse")) {
			 /* end channel */
			 tcollapse=atoi(argv[++i]);

		  } else if (strings_equal(argv[i],"-fcollapse")) {
			 /* end channel */
			 fcollapse=atoi(argv[++i]);


		  } else {

			sum_fil_help();
			/* unknown argument passed down - stop! */
			sprintf(message,"unknown argument (%s) passed to filterbank.",argv[i]);
			error_message(message);
			exit(1);
		}
		i++;
	  }
  }

  if (!nfiles) {
	  error_message("no input files supplied on command line!");
	  exit(1);
  }

  if (nfiles > 50) {
      error_message("too many input files supplied on command line (max 50)!");
	  exit(1);
  }
  
  if (!opened) {
	  error_message("must have an output file (-o <output>)!");
	  exit(1);
  }




	pulsarcentric=barycentric=0;
	writeobsdbline=0;
	


    
  for(i = 0; i<nfiles; i++) {
		headersize=read_header(fileptr[i]);		
	    rewind(fileptr[i]);
		datasize=sizeof_file(argv[i+1])-headersize;
	    numsamps=nsamples(argv[i+1],headersize,nbits,nifs,nchans);

		if(output_datasize == 0 || datasize < output_datasize) {
			output_datasize = datasize;
			input_datasize = datasize;
		    output_headersize = headersize;

			output_ptr = i;

		}
		if(output_nchans == 0) {
		  output_nchans = nchans/fcollapse; 
		  input_nchans = nchans;
		  input_nbits = nbits; 
		  output_nifs = nifs;
		  input_nifs = nifs;		
		}
        if(ceil(nchans/fcollapse) != nchans/fcollapse) {
			sprintf(message,"nchan not evenly divisible by fcollapse factor!");
			error_message(message);
			exit(1);
        }
        
        printf("outbits: %d inputbits: %d file: %s headersize: %Ld nbits: %d nifs: %d nchans: %d\n", output_nbits, input_nbits, argv[i+1],headersize, nbits, nifs, nchans); 
  }


  if (output_nbits == 8 && input_nbits == 8) {
		  
		 
		  printf("minimum data size is: %Ld\n", output_datasize);
		  
		  printf("will dump: %Ld\n", (output_datasize / (long long) (output_nifs * output_nchans)));
		  
		  headersize=read_header(fileptr[output_ptr]);
		  rewind(fileptr[output_ptr]);
		
		  printf("header size lead file: %Ld\n", headersize); 
		  fread(buffer, sizeof(char), headersize, fileptr[output_ptr]);
		  fwrite(buffer, sizeof(char), headersize, output);
		  rewind(fileptr[output_ptr]);
		
		 /* bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				headersize=read_header(fileptr[i]);	
			}
		  
		
		/* outer loop, read through 8000 spectra (1024000 spectral points), histogram, calculate sigma for quantizing to 8 bits */
		
		  for (j = 0; j < 8000; j++){
		  
			 /* read n spectra (1 spectra x n files), sum */
			 for(i=0;i<128;i++) spectra_sum[i] = 0.0;
			 for(i = 0; i<nfiles; i++) {
				  fread(buffer, sizeof(char), (output_nifs * output_nchans), fileptr[i]);
				  for(k=0; k<128; k++) spectra_sum[k] = (spectra_sum[k] + ((double) buffer[k]));
			  }
			  //for(i=0; i<128; i++) printf("%f\n", spectra_sum[i]);
			  for(i=0; i<128; i++) gsl_histogram_increment(spectra_quant, spectra_sum[i]);
			  //usleep(5000000);	
		   }	
				
		sigma = gsl_histogram_sigma(spectra_quant);
		mean = gsl_histogram_mean(spectra_quant);
		
		gsl_histogram_free (spectra_quant);
		
		 /* rewind and bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				rewind( fileptr[i] );
				headersize=read_header(fileptr[i]);	
			}
		
		
		/* outer loop, read through all spectra, quantize to 8 bits, write to file */
		
		  for (j = 0; j < (output_datasize / (long long) (output_nifs * output_nchans)); j++){
		  
			 /* read n spectra (1 spectra x n files), sum */
			 for(i=0;i<128;i++) spectra_sum[i] = 0.0;
			 for(i = 0; i<nfiles; i++) {
				  fread(buffer, sizeof(char), (output_nifs * output_nchans), fileptr[i]);
				  for(k=0; k<128; k++) spectra_sum[k] = (spectra_sum[k] + ((double) buffer[k]));
			  }
			  for(i=0; i<128; i++) {	  
				  quantval = quantize(spectra_sum[i], mean, sigma);
				  fwrite(&quantval, sizeof(char), 1, output);	  
				  //printf("%f : %d\n", spectra_sum[i], quantval);
			  }
				  //usleep(5000000);
		
		   }

	} else if (output_nbits == 32 && input_nbits == 32) {

	  printf("input 32, output 32\n");
		  printf("minimum data size is: %Ld\n", output_datasize);
		  
		  printf("will dump: %Ld\n", (output_datasize / (long long) (output_nifs * output_nchans)));
		  

		  rewind(fileptr[output_ptr]);	
		  headersize=read_header(fileptr[output_ptr]);
		
		  machine_id = 10;
		  printf("header size lead file: %Ld\n", headersize); 
		  nbits=32;
		  obits=32;
		  tsamp = tsamp * (double) tcollapse;
		  foff = foff * (double) fcollapse;
		  
		  nchans = output_nchans;
		  strcpy(ifstream,"YYYY");
		  
		  //hanning=hamming=zerolagdump=swapout=sumifs=headerless=headerfile=0;
		  //invert_band=clip_threshold=headeronly=0;
		
		  filterbank_header(output);
		
		
		 /* rewind and bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				rewind( fileptr[i] );
				headersize=read_header(fileptr[i]);	
			}
		
		
		/* outer loop, read through all spectra, quantize to 8 bits, write to file */
		
		  for (j = 0; j < ((input_datasize/4) / (long long) (input_nifs * input_nchans)); j++){


		     memset(fspectra_sum, 0x0, output_nchans * sizeof(float));
			
			 for(i = 0; i < tcollapse; i++) {
				   fread(fbuffer, sizeof(float), (input_nifs * input_nchans), fileptr[0]);
				  
				   for(k = 0; k < nchans; k = k + fcollapse) {   		
				   		for (m = 0; m < fcollapse; m++) fspectra_sum[k/fcollapse] = (fspectra_sum[k/fcollapse] + fbuffer[k + m]);
				   }
				   
			 }					
			

			fwrite(fspectra_sum, sizeof(float), output_nchans, output);	 				 
			 			   		
		   }

		} else if (output_nbits == 32 && input_nbits == 8) {


		  //output_datasize = output_datasize; // cut data output down by a factor of 4;
		  
		  printf("minimum data size is: %Ld\n", output_datasize);
		  
		  printf("will dump: %Ld\n", (output_datasize / (long long) (output_nifs * output_nchans)));
		  
		  headersize=read_header(fileptr[output_ptr]);
		  rewind(fileptr[output_ptr]);
		
		  printf("header size lead file: %Ld\n", headersize); 
		  nbits=32;
		  obits=32;
		  strcpy(ifstream,"YYYY");
		  filterbank_header(output);

		  fread(buffer, sizeof(char), headersize, fileptr[output_ptr]);
		  //fwrite(buffer, sizeof(char), headersize, output);
		  rewind(fileptr[output_ptr]);
		
		 /* bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				headersize=read_header(fileptr[i]);	
			}
		  
		
		/* outer loop, read through 8000 spectra (1024000 spectral points), histogram, calculate sigma for quantizing to 8 bits */
		
		  for (j = 0; j < 8000; j++){
		  
			 /* read n spectra (1 spectra x n files), sum */
			 for(i=0;i<128;i++) spectra_sum[i] = 0.0;
			 for(i = 0; i<nfiles; i++) {
				  fread(buffer, sizeof(char), (output_nifs * output_nchans), fileptr[i]);
				  for(k=0; k<128; k++) spectra_sum[k] = (spectra_sum[k] + ((double) buffer[k]));
			  }
			  //for(i=0; i<128; i++) printf("%f\n", spectra_sum[i]);
			  for(i=0; i<128; i++) gsl_histogram_increment(spectra_quant, spectra_sum[i]);
			  //usleep(5000000);	
		   }	
				
		sigma = gsl_histogram_sigma(spectra_quant);
		mean = gsl_histogram_mean(spectra_quant);
		
		gsl_histogram_free (spectra_quant);
		
		 /* rewind and bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				rewind( fileptr[i] );
				headersize=read_header(fileptr[i]);	
			}
		
		
		/* outer loop, read through all spectra, quantize to 8 bits, write to file */
		
		  for (j = 0; j < (output_datasize / (long long) (output_nifs * output_nchans)); j++){
		  
		  	
			   
			 /* read n spectra (1 spectra x n files), sum */
			 for(i=0;i<128;i++) fspectra_sum[i] = 0.0;
			 for(i = 0; i<nfiles; i++) {
				  fread(buffer, sizeof(char), (output_nifs * output_nchans), fileptr[i]);
				  for(k=0; k<128; k++) fspectra_sum[k] = (fspectra_sum[k] + ((double) buffer[k]));
			  }
			
			  fwrite(&fspectra_sum, sizeof(float), (output_nifs * output_nchans), output); 
				  //usleep(5000000);
		
		   }


		} else if (output_nbits == 8 && input_nbits == 32) {
		/* for now, this is the only mode guaranteed to work!  AS 7/12 */



		  printf("input 32, output 8\n");
		  printf("minimum data size is: %Ld\n", output_datasize);
		  
		  printf("will dump: %Ld\n", (output_datasize / (long long) (output_nifs * output_nchans)));
		  

		  rewind(fileptr[output_ptr]);	
		  headersize=read_header(fileptr[output_ptr]);
		
		  machine_id = 10;
		  printf("header size lead file: %Ld\n", headersize); 
		  nbits=8;
		  obits=8;
		  tsamp = tsamp * (double) tcollapse;
		  foff = foff * (double) fcollapse;
		  
		  nchans = output_nchans;
		  strcpy(ifstream,"YYYY");
		  
		  //hanning=hamming=zerolagdump=swapout=sumifs=headerless=headerfile=0;
		  //invert_band=clip_threshold=headeronly=0;
		
		  filterbank_header(output);
		
		 /* bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
		   		rewind(fileptr[i]);
				headersize=read_header(fileptr[i]);	
			}
		  
		/*
		I was thinking of reading in all the values used to compute the quantization, excluding all the zeros, 
		calculate the min and the max, histogram the values, and then calculate the threshold values that throw away 
		0.001% at the top and 0.001% at the bottom the problem with that
		*/
		
		
		/* step 1: read through quantization region and identify min and max values */
		  /* loop over all the spectra in the quantization region */
		  /* NEED TO CHECK IF qlen > totlen */
		  min = 0;
		  max = 0;
		  
		  if (chanend == 0) chanend = input_nchans; // if chanend not specified, set equal to nchans
		  
		  fbuffer = (float *) malloc (input_nchans * sizeof(float));
		  fspectra_sum = (float *) malloc (output_nchans * sizeof(float));

		  
		  printf("inp_nifs: %d  output_nifs: %d output_nchans %d nfiles %d qlen %Ld input_nchans %d, output_nchans %d\n", input_nifs,output_nifs, output_nchans, nfiles, qlen, input_nchans, output_nchans);
		  for (j = 0; j < qlen; j++){		  
			 /* read n spectra (1 spectra x n files), sum */
			  /* reset spectra_sum to zero */			
		     memset(fspectra_sum, 0x0, output_nchans * sizeof(float));
		
			 for(i = 0; i < tcollapse; i++) {
				   fread(fbuffer, sizeof(float), (input_nifs * input_nchans), fileptr[0]);
			   
				   for(k = 0; k < input_nchans; k = k + fcollapse) {
				   		for (m = 0; m < fcollapse; m++) fspectra_sum[k/fcollapse] = (fspectra_sum[k/fcollapse] + fbuffer[k + m]);
				   }

			 }			

			 for(k=0; k < output_nchans; k++) {		
				 //printf("%f\n", );
				 if(fspectra_sum[k] > max) max = fspectra_sum[k];
				 if((fspectra_sum[k] < min && fspectra_sum[k] > 0.0) || min == 0.0) min = fspectra_sum[k];				  	
			 }
			  
		   }	

		printf("min: %f max: %f\n", min, max);		

	   
		for(i = 0; i<nfiles; i++) {
			 rewind(fileptr[i]);
			 headersize=read_header(fileptr[i]);	
		}
		


		spectra_quant = gsl_histogram_alloc ((int) (ceil(max) - floor(min)));
		gsl_histogram_set_ranges_uniform (spectra_quant, (int) floor(min), (int) ceil(max));

		
		printf("histogram configured with size: %d and bounds %d %d\n",(int) (ceil(max) - floor(min)),(int) floor(min), (int) ceil(max));		

		  for (j = 0; j < qlen; j++){		  
			 /* read n spectra (1 spectra x n files), sum */
			  /* reset spectra_sum to zero */			
		     memset(fspectra_sum, 0x0, output_nchans * sizeof(float));
			
			 for(i = 0; i < tcollapse; i++) {
				   fread(fbuffer, sizeof(float), (input_nifs * input_nchans), fileptr[0]);
				  
				   for(k = 0; k < input_nchans; k = k + fcollapse) {   		
				   		for (m = 0; m < fcollapse; m++) fspectra_sum[k/fcollapse] = (fspectra_sum[k/fcollapse] + fbuffer[k + m]);
				   }
				   
			 }			

			 for(k=0; k < output_nchans; k++) {		
				 if(fspectra_sum[k] > 0.0) {
				 	gsl_histogram_increment(spectra_quant, fspectra_sum[k]);			 
			 		nelements++;
			 	}
			 }
			  
		   }	
		
		  		  
		  
		 printf("got: %f %f in %f elements\n", gsl_histogram_sigma(spectra_quant), gsl_histogram_mean(spectra_quant), nelements);
		
		
		
		/* here we identify the bins corresponding to the inner fraction 1 - 2 x discard */
		
		 i=0;
		 powersum = 0.0;
		 while(powersum/nelements < discard) {
		 	powersum = powersum + gsl_histogram_get(spectra_quant, i);		 
		 	i++;
		 }


		 gsl_histogram_get_range(spectra_quant, i - 1, &lower, &upper);		 
		 quantmin = (float) lower;
		 
		 		 
		  printf("bottom thresh bin %ld val %f\n", i, quantmin);

		 i = (long int) (ceil(max) - floor(min) - 1);
		 powersum = 0.0;
		 while(powersum/nelements < discard) {
		 	powersum = powersum +  gsl_histogram_get(spectra_quant, i);		 
		 	i--;
		 }

		 gsl_histogram_get_range(spectra_quant, i + 1, &lower, &upper);		 
		 quantmax = (float) upper;

		  printf("top thresh bin %ld val %f\n",i, quantmax);

		  gsl_histogram_free (spectra_quant);
		 
		 
		 
		
		 /* rewind and bump past header for all input files */
		   for(i = 0; i<nfiles; i++) {
				rewind( fileptr[i] );
				headersize=read_header(fileptr[i]);	
			}
		
		
		/* outer loop, read through all spectra, quantize to 8 bits, write to file */
		
		  for (j = 0; j < ((input_datasize/4) / (long long) (input_nifs * input_nchans)); j++){


		     memset(fspectra_sum, 0x0, output_nchans * sizeof(float));
			
			 for(i = 0; i < tcollapse; i++) {
				   fread(fbuffer, sizeof(float), (input_nifs * input_nchans), fileptr[0]);
				  
				   for(k = 0; k < nchans; k = k + fcollapse) {   		
				   		for (m = 0; m < fcollapse; m++) fspectra_sum[k/fcollapse] = (fspectra_sum[k/fcollapse] + fbuffer[k + m]);
				   }
				   
			 }					
			
			 for(k=0; k < output_nchans; k++) {		
				  quantval = quantize(fspectra_sum[k], quantmin, quantmax);
				  fwrite(&quantval, sizeof(char), 1, output);	 				 
			 }
			 			   		
		   }

  	
  	} else {

	  error_message("This routine only works on 8 bit or 32 bit filterbank data!");
	  exit(1);
  	
  	}



/*
	if (argc>1) {
		print_version(argv[0],argv[1]);
		if (help_required(argv[1])) {
			header_help();
			exit(0);
		} else if (file_exists(argv[1])) {
			strcpy(filename,argv[1]);
			fileptr=open_file(filename,"rb");
		} else if (!file_exists(argv[1]) && (strncmp(argv[1],"-",1) !=0)) {
			sprintf(message,"Data file: %s not found...\n",argv[1]);
			error_message(message);
			exit(1);
		}
	}
*/

}
