#include <stdio.h>
#include <string.h>
#include "filterbank.h"

/*
Mars Lightning pcap processor - read in a mars lightning pcap capture and output filterbank format 

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
  FILE *input, *xpowerfile, *ypowerfile, *xpowersqfile, *ypowersqfile, *powersumfile, *powersqsumfile;
  int   i,j,k,l;
  char file[120];


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
  machine_id=12;    //-1
  telescope_id=10;  //-1
  data_type=1;
  nchans=1024;
  foff=104.8576/nchans * -1.0;        // XXX hack
  if (fcentral==0) {
  	fcentral=1420.00;
  	printf("warning... center freq not set, defaulting to 1420 MHz\n");
  }
  fch1=fcentral  - ((nchans/2)+0.5)*foff;    // A bit of a guess
  if (tsamp==0) {
	  tsamp=0.00125;
      printf("warning... sample time not set, defaulting to 1.25 milliseconds\n");  		
  }
  nifs=1;
  obits = 8;
  nbits = 8;
    
    
      /* open up input file */
    strcpy(file,argv[1]);
    input=fopen(file,"rb");

  
  	/* open up output files */
    
    sprintf(file,"%s_xpower.fil", file_prepend);
    xpowerfile=fopen(file,"wb");

    sprintf(file,"%s_ypower.fil", file_prepend);
    ypowerfile=fopen(file,"wb");

    sprintf(file,"%s_xpowersq.fil", file_prepend);
    xpowersqfile=fopen(file,"wb");

    sprintf(file,"%s_ypowersq.fil", file_prepend);
    ypowersqfile=fopen(file,"wb");
    
    



 	if(sumall) {

		sprintf(file,"%s_powersum.fil", file_prepend);
		powersumfile=fopen(file,"wb");
	
		sprintf(file,"%s_powersqsum.fil", file_prepend);
		powersqsumfile=fopen(file,"wb");
 	}


/* bump past pcap global header */
bump = fread(&payload, sizeof(char), 24, input);

do {
            bump = fread(&payload, sizeof(char), 8, input); /* skip 8 bytes */

			/* these should be equal, but sometimes tcpdump/gulp/etc doesn't work quite right? */
            bump = fread(&packet_saved_len, sizeof(int), 1, input);           
            bump = fread(&packet_len, sizeof(int), 1, input);            


                 if((packet_len == packet_saved_len) && (packet_len == 8242) && (bump != 0)) {
 					
 					speccnt++;
 					bump = fread(&payload, sizeof(char), packet_len, input);
									 
					/* skip past ip and udp headers, read event records */ 
					memcpy(&pkt_index, payload+42, 8);

					/* first 64 bits is the clock cycle counter value */
					pkt_index = bswap_64(pkt_index);
					

					pkt_index_prev = pkt_index;					
					for (j = (8242 - 8192);j<8242;j=j+8){
						  //copy one event recorc from the payload						 
						  memcpy(&fields, payload + j, 8);						  
						  
						  xpolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,0);
						  ypolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,24);

						  xpolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,32));
 						  ypolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,48));						  
						  
						  xpowersqmean += (double) xpolsquared[(j - (8242 - 8192))/8];
						  ypowersqmean += (double) ypolsquared[(j - (8242 - 8192))/8];
						  						  			
					 }
										 
                 } else if ( bump !=0 ){
                    //skip this packet 
					//fprintf(stderr, "pcap file contains a non-event-record packet - trying to skip\n");
					bump = fread(&payload, sizeof(char), packet_saved_len, input);
					//fprintf(stderr, "got: %d - %d\n", packet_len, packet_saved_len);

                }   
} while (speccnt != 1024);

xpowersqmean = xpowersqmean / (1024*1024); //divide by 1024 1024 channel spectra
ypowersqmean = ypowersqmean / (1024*1024); //divide by 1024 1024 channel spectra

speccnt=0;
rewind(input);
pkt_index_prev=0;
/* bump past pcap global header */
bump = fread(&payload, sizeof(char), 24, input);

do {
            bump = fread(&payload, sizeof(char), 8, input); /* skip 8 bytes */

			/* these should be equal, but sometimes tcpdump/gulp/etc doesn't work quite right? */
            bump = fread(&packet_saved_len, sizeof(int), 1, input);           
            bump = fread(&packet_len, sizeof(int), 1, input);            


                 if((packet_len == packet_saved_len) && (packet_len == 8242) && (bump != 0)) {
 					
 					speccnt++;
 					bump = fread(&payload, sizeof(char), packet_len, input);
									 
					/* skip past ip and udp headers, read event records */ 
					memcpy(&pkt_index, payload+42, 8);

					/* first 64 bits is the clock cycle counter value */
					pkt_index = bswap_64(pkt_index);
					

					pkt_index_prev = pkt_index;					
					for (j = (8242 - 8192);j<8242;j=j+8){
						  //copy one event recorc from the payload						 
						  memcpy(&fields, payload + j, 8);						  
						  
						  xpolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,0);
						  ypolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,24);
						  xpolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,32));
 						  ypolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,48));						  
						  
						  //printf("%f %f\n", (double) xpolsquared[(j - (8242 - 8192))/8], (double) ypolsquared[(j - (8242 - 8192))/8]);
					 	  xpowersqstd += pow( ( ((double) xpolsquared[(j - (8242 - 8192))/8]) - xpowersqmean) , 2);
					 	  ypowersqstd += pow( ( ((double) ypolsquared[(j - (8242 - 8192))/8]) - ypowersqmean) , 2);
	  						  			
					 }
										 
                 } else if ( bump !=0 ){
                    //skip this packet 
					//fprintf(stderr, "pcap file contains a non-event-record packet - trying to skip\n");
					bump = fread(&payload, sizeof(char), packet_saved_len, input);
					//fprintf(stderr, "got: %d - %d\n", packet_len, packet_saved_len);
                }   
} while (speccnt != 1024);

xpowersqstd = pow((xpowersqstd / ((1024*1024)-1)), 0.5); //divide by 1024 1024 channel spectra
ypowersqstd = pow((ypowersqstd / ((1024*1024)-1)), 0.5); 


printf("xpsqmean %f ypsqmean %f xpsqstd %f ypsqstd %f\n", xpowersqmean, ypowersqmean, xpowersqstd, ypowersqstd);
speccnt=0;
rewind(input);
pkt_index_prev=0;

    /* bump past pcap global header */
    bump = fread(&payload, sizeof(char), 24, input);

	//printf("size of unsigned long int is %ld and should be 4\n\n", sizeof(unsigned long long int));


do {
            bump = fread(&payload, sizeof(char), 8, input); /* skip 8 bytes */

			if(!first_packet) {
				 ts_sec = (((unsigned int) payload[0])%256);
				 ts_sec += (((unsigned int) payload[1])%256) * 256;
				 ts_sec += (((unsigned int) payload[2])%256) * 65536;
				 ts_sec += (((unsigned int) payload[3])%256) * 16777216;             
		
				 ts_usec = (((unsigned int) payload[4])%256);
				 ts_usec += (((unsigned int) payload[5])%256) * 256;
				 ts_usec += (((unsigned int) payload[6])%256) * 65536;
				 ts_usec += (((unsigned int) payload[7])%256) * 16777216;
		         if(tstart == 0) tstart=40587.0+((double)ts_sec+((double)ts_usec)/1000000.0)/86400.0;
				 printf("tstart: %f\n", tstart);
				 
				 first_packet = 1;
				 filterbank_header(xpowerfile);
				 filterbank_header(ypowerfile);

				 /* go to 32 bits for the power sq values */
				 obits = 32;
				 nbits = 32;
				 filterbank_header(xpowersqfile);
				 filterbank_header(ypowersqfile);
				 if(sumall) {
					  filterbank_header(powersumfile);
					  filterbank_header(powersqsumfile); 	
				 }

			}

			/* these should be equal, but sometimes tcpdump/gulp/etc doesn't work quite right? */
            bump = fread(&packet_saved_len, sizeof(int), 1, input);           
            bump = fread(&packet_len, sizeof(int), 1, input);            

                 //printf("packet len is %d\n", packet_len);
                 //if we received a 4138 byte packet, its probably a bee2 UDP event dump

                 if((packet_len == packet_saved_len) && (packet_len == 8242) && (bump != 0)) {
 					
 					speccnt++;
 					bump = fread(&payload, sizeof(char), packet_len, input);
				
					 
					/* skip past ip and udp headers, read event records */ 
					memcpy(&pkt_index, payload+42, 8);

					/* first 64 bits is the clock cycle counter value */
					pkt_index = bswap_64(pkt_index);
					
					//printf("%ld\n", pkt_index);
					//printf("packet index delta: %ld\n", (pkt_index - pkt_index_prev)/accum);
					//(pkt_index_prev != 0) && ((pkt_index - pkt_index_prev) != accum)
					if((pkt_index_prev != 0) && ((pkt_index - pkt_index_prev) != accum)){
               				
						  dropped_packet_cnt += ((pkt_index - pkt_index_prev)/accum - 1);
						  err = (((int) pkt_index - (int) pkt_index_prev));
						  
						  
						  /* pad file for dropped spectra... */
						  
						  for(l = 0; l < ((pkt_index - pkt_index_prev)/accum - 1);l++){
							  
							  if(!reversed) { 
							  
								 for(k=512;k<1024;k++){
										  fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
										  fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
										  quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
										  fwrite(&quantval, sizeof(char), 1, xpowersqfile);
										  quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
										  fwrite(&quantval, sizeof(char), 1, ypowersqfile);
								 }
			   
								 for(k=0;k<512;k++){
										 fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
										 fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
										 quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
										 fwrite(&quantval, sizeof(char), 1, xpowersqfile);
										 quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
										 fwrite(&quantval, sizeof(char), 1, ypowersqfile);
								 }
							  
							  } else {
							  
									 /* --------- output reversed band ---------------*/
									 for(k=511;k>=0;k--){
										  fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
										  fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
										  quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
										  fwrite(&quantval, sizeof(char), 1, xpowersqfile);
										  quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
										  fwrite(&quantval, sizeof(char), 1, ypowersqfile);
									 }					
									 
									 for(k=1023;k>=512;k--){
										  fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
										  fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
										  quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
										  fwrite(&quantval, sizeof(char), 1, xpowersqfile);
										  quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
										  fwrite(&quantval, sizeof(char), 1, ypowersqfile);
									 }					
							  
							  }
	  
							  printf("dropped pkts: %d at spec count (strts at 1): %d (corrected)\n", err,speccnt);
						  
						   }
					}
					pkt_index_prev = pkt_index;					
					for (j = (8242 - 8192);j<8242;j=j+8){
						  //copy one event recorc from the payload						 
						  memcpy(&fields, payload + j, 8);
						  
						  
						  xpolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,0);
						  ypolpower[(j - (8242 - 8192))/8] = (unsigned char) slice64(fields,8,24);
						  xpolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,32));
 						  ypolsquared[(j - (8242 - 8192))/8] = (unsigned short) bswap_16(slice64(fields,16,48));

						  //extract record components
						  //event = slice(fields,1,31);
							
						  //pfb_fft_power = ntohl(pfb_fft_power);
						  
						  /*
						  hex dump
						  */
						  
						  
						  //fwrite(&xpolpower, sizeof(char), 1, xpowerfile);
						  //fwrite(&ypolpower, sizeof(char), 1, ypowerfile);
						  //fwrite(&xpolsquared, sizeof(short), 1, xpowersqfile);
						  //fwrite(&ypolsquared, sizeof(short), 1, ypowersqfile);

 			             //printf("X: %u Y: %u X2: %hu Y2: %hu \t\t", xpolpower[(j - (8242 - 8192))/8], ypolpower[(j - (8242 - 8192))/8], xpolsquared[(j - (8242 - 8192))/8], ypolsquared[(j - (8242 - 8192))/8]);

						  /*						  
						  printf("X: 0x%x Y: 0x%x X2: %hu Y2: %hu \t\t", xpolpower, ypolpower, xpolsquared, ypolsquared); 

						  printf("0x%02X", (unsigned char) payload[j]);
						  printf("%02X", (unsigned char) payload[j+1]);
						  printf("%02X", (unsigned char) payload[j+2]);
						  printf("%02X", (unsigned char) payload[j+3]);						  
						  printf("%02X", (unsigned char) payload[j+4]);
						  printf("%02X", (unsigned char) payload[j+5]);
						  printf("%02X", (unsigned char) payload[j+6]);
						  printf("%02X\n", (unsigned char) payload[j+7]);
						  */
						  
						  //printf("0x%16X\n", (unsigned long int) fields);
						  
						  
						  
					 }
					
					 //fft_bin = (fft_bin + 512) % 1024;
					 
					 //Reorder frequency bins in regular order (low frequencies first)
					 if(!reversed) {

							for(k=512;k<1024;k++){
								 fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
								 fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);								 
								 //quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
								 //fwrite(&quantval, sizeof(char), 1, xpowersqfile);
								 //quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
								 //fwrite(&quantval, sizeof(char), 1, ypowersqfile);
								 floatval = (float) xpolsquared[k]; 
								 fwrite(&floatval, sizeof(float), 1, xpowersqfile);
								 floatval = (float) ypolsquared[k]; 								  								 
								 fwrite(&floatval, sizeof(float), 1, ypowersqfile);
							}					
							
							for(k=0;k<512;k++){
								 fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
								 fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
								 //quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
								 //fwrite(&quantval, sizeof(char), 1, xpowersqfile);
								 //quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
								 //fwrite(&quantval, sizeof(char), 1, ypowersqfile);
								 floatval = (float) xpolsquared[k]; 
								 fwrite(&floatval, sizeof(float), 1, xpowersqfile);
								 floatval = (float) ypolsquared[k]; 								  								 
								 fwrite(&floatval, sizeof(float), 1, ypowersqfile);
							}					
					} else {

							/* --------- output reversed band ---------------*/

							for(k=511;k>=0;k--){
								 fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
								 fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
								 //quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
								 //fwrite(&quantval, sizeof(char), 1, xpowersqfile);
								 //quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
								 //fwrite(&quantval, sizeof(char), 1, ypowersqfile);
								 floatval = (float) xpolsquared[k]; 
								 fwrite(&floatval, sizeof(float), 1, xpowersqfile);
								 floatval = (float) ypolsquared[k]; 								  								 
								 fwrite(&floatval, sizeof(float), 1, ypowersqfile);
							}					
							
							for(k=1023;k>=512;k--){
								 fwrite(xpolpower + k, sizeof(char), 1, xpowerfile);
								 fwrite(ypolpower + k, sizeof(char), 1, ypowerfile);
	
								 //quantval = quantize((double) xpolsquared[k], xpowersqmean, xpowersqstd); 
								 //fwrite(&quantval, sizeof(char), 1, xpowersqfile);
								 //quantval = quantize((double) ypolsquared[k], ypowersqmean, ypowersqstd); 								 
								 //fwrite(&quantval, sizeof(char), 1, ypowersqfile);

								 floatval = (float) xpolsquared[k]; 
								 fwrite(&floatval, sizeof(float), 1, xpowersqfile);
								 floatval = (float) ypolsquared[k]; 								  								 
								 fwrite(&floatval, sizeof(float), 1, ypowersqfile);

							}					
							
															
					
					}

                 } else if ( bump !=0 ){
                    //skip this packet 
					//fprintf(stderr, "pcap file contains a non-event-record packet - trying to skip\n");
					bump = fread(&payload, sizeof(char), packet_saved_len, input);
					//fprintf(stderr, "got: %d - %d\n", packet_len, packet_saved_len);

                }   
} while (bump != 0);

/*
If you want to grab the packet time stamp, you might do something like this...

    bump = fread(&payload, sizeof(char), 4, input);
      ts_sec = (((unsigned int) payload[0])%256);
      ts_sec += (((unsigned int) payload[1])%256) * 256;
      ts_sec += (((unsigned int) payload[2])%256) * 65536;
      ts_sec += (((unsigned int) payload[3])%256) * 16777216;    

    bump = fread(&payload, sizeof(char), 4, input);
      ts_usec = (((unsigned int) payload[0])%256);
      ts_usec += (((unsigned int) payload[1])%256) * 256;
      ts_usec += (((unsigned int) payload[2])%256) * 65536;
      ts_usec += (((unsigned int) payload[3])%256) * 16777216;
*/




  /*fclose(output);*/
  fprintf(stderr,"Completed Successfully\n");
  fprintf(stderr, "Spectra count: %d\n", speccnt);
  fprintf(stderr, "Dropped packet count: %d\n", dropped_packet_cnt);

  
  fclose(xpowerfile);
  fclose(ypowerfile);
  fclose(xpowersqfile);
  fclose(ypowersqfile);

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

