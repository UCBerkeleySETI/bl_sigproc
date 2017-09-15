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

int wapp_isalfa;

main (int argc, char *argv[])
{
  int   i,j,k,nfiles,fileidx,fileidx_start,inputdata,opened=0;
  char  message[80], chan_name[4];
  FILE  *mult_out[11][4];
  FILE  *sum_output;
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
  unsigned int block[128];
  float        equalize_coeff[128][11][4], average_power=0; 
  double tacc[11], tdelta[11]; 
  float        temp_power=0.0;
  unsigned int bump, packet_len, packet_saved_len, ip, s;
  unsigned int oldaccumulation[11], accumulation_number[11];
  int          equalize_init=0, first_packet=1, correcterror=0, statonly=0;
  int          ts_sec, ts_usec; /* packet timestamp seconds, microcseconds */
  long long    counter[11], prev_counter, last_counter[11];
  double spectra_sum;
  int delta;
  int delta_counter[11];
  double average=0.0;  
  double start_times[11];
  double end_times[11];
  double start_cnt[11];
  double end_cnt[11];
  
  double skip_acc[11];
  double skip_count[11];
  double skip_old_acc[11];
  
  int antenna_select = -1;
  int ibob_select = -1;

  for (i=0; i<11; i++) { start_times[i]=end_times[i]=start_cnt[i]=end_cnt[i]=-1; }
  
  
  /* set and clear arrays */
  strcpy(chan_name,"DCBA"); /* A=3 B=2 C=1 D=0 */
  strcpy(file_prepend, "");
  
  for (i=0; i<11; i++) { accumulation_number[i]=counter[i]=0; oldaccumulation[i]=-1; }

  for (i=0; i<11; i++) delta_counter[i] = 0;
  
  memset(block,        0, 128*     sizeof(unsigned int));
  memset(mask,         0, 11*4*sizeof(int));
  memset(tacc,         0, 11 * sizeof(double));
  memset(tdelta,         0, 11 * sizeof(double));
  memset(skip_acc,         0, 11 * sizeof(double));
  memset(skip_old_acc,         0, 11 * sizeof(double));
  memset(skip_count,         0, 11 * sizeof(double));
  /* -- mask? -- */
  /* crab observation: */
  /*   2B         2D         4D         5A         5C         7B         7D         8A         8D         10C         10D  */  
//  mask[2][2]=mask[2][0]=mask[4][0]=mask[5][3]=mask[5][1]=mask[7][2]=mask[7][0]=mask[8][3]=mask[8][0]=mask[10][1]=mask[10][0]=1;

  /* mask iBOB 38 */
  mask[0][0]=mask[0][1]=mask[0][2]=mask[0][3]=1;

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
  wapp_isalfa=invert_band=clip_threshold=headeronly=0;
  time_offset=start_time=final_time=tstart=0.0;
  obits=-1;
  tsamp=0.0;
  do_vanvleck=compute_spectra=1;
  strcpy(ifstream,"XXXX");
  machine_id=11;    //-1
  telescope_id=10;  //-1
  data_type=1;
  

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
      } else if (strings_equal(argv[i],"-c")) {
        /* get clip threshold (sigma) */
        clip_threshold=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-s")) {
        /* get starting time (s) */
        start_time=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-r")) {
        /* get time to read (s) this is adjusted below if skipping */
        final_time=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-n")) {
        /* output number of bits per sample to write */
        obits=atoi(argv[++i]);
      } else if (strings_equal(argv[i],"-dt")) {
        /* add a time offset in seconds to tstart */
        time_offset=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-i")) {
        /* flag IBOB to write */
        i++;
        if (atoi(argv[i])<0 || atoi(argv[i])>10) {
          printf("IFstream must lie between 0 and 10\n");
          return;
        }
        ibob_select = atoi(argv[i]);
      } else if (strings_equal(argv[i],"-ant")) {
        /* flag IF stream to write */
        i++;
        if (atoi( argv[i])<0 || atoi(argv[i])>3) {
          printf("antenna stream must lie between 0 and 3\n");
          return;
        }
        antenna_select = atoi(argv[i]);
      } else if (strings_equal(argv[i],"-swapout")) {
        /* perform byte swapping on all output data */
        swapout=1;
      } else if (strings_equal(argv[i],"-floats")) {
        /* write data as floating point numbers */
        obits=32;
      } else if (strings_equal(argv[i],"-sumifs")) {
        /* sum IFs if necessary */
        sumifs=1;
      } else if (strings_equal(argv[i],"-zerolag")) {
        /* zerolagdump used for correlators e.g. WAPP */
        zerolagdump=1;
        obits=32;
      } else if (strings_equal(argv[i],"-rawcfs")) {
        /* write correlation functions only */
        compute_spectra=do_vanvleck=0;
      } else if (strings_equal(argv[i],"-corcfs")) {
        /* write corrected correlation functions */
        compute_spectra=0;
        do_vanvleck=1;
      } else if (strings_equal(argv[i],"-novanvleck")) {
        /* don't apply van vleck correction */
        do_vanvleck=0;
      } else if (strings_equal(argv[i],"-invert")) {
        /* invert the band after FFT */
        invert_band=1;
      } else if (strings_equal(argv[i],"-hamming")) {
        /* Hamming smoothing */
        hamming=1;
        hanning=0;
      } else if (strings_equal(argv[i],"-hanning")) {
        /* Hanning smoothing */
        hanning=1;
        hamming=0;
      } else if (strings_equal(argv[i],"-headerfile")) {
        /* no binary headers but write data to "head" file */
        headerless=headerfile=1;
      } else if (strings_equal(argv[i],"-headeronly")) {
        /* only binary header written */
        headeronly=1;    
      /* ata filterbank options */
      } else if (strings_equal(argv[i],"-mjd")) {
          /* get the fractional start mjd */
          tstart=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-skip2")) {
        /* skip spectra 2 for funny 1pps stuff */
        skip_packet_2=1; 
      } else if (strings_equal(argv[i],"-prepend")) {
        strcpy(file_prepend,argv[++i]);
      } else if (strings_equal(argv[i],"-ra")) {
        src_raj=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-dec")) {
        src_dej=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-az")) {
        az_start=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-el")) {
        za_start=atof(argv[++i]);       
      } else if (strings_equal(argv[i],"-machine_id")) {
        machine_id=atoi(argv[++i]);       
      } else if (strings_equal(argv[i],"-telescope_id")) {
        telescope_id=atoi(argv[++i]);       
      } else if (strings_equal(argv[i],"-correct")) {
        /* attempt to correct dropped packet errors for periodicity searches */
        correcterror=1;
      } else if (strings_equal(argv[i],"-tsamp")) {
        /* set sampling time */
        tsamp=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-statonly")) {
        /* only produce statistics */
        statonly = 1;
      } else {
        /* unknown argument passed down - stop! */
          filterbank_help();

        sprintf(message,"unknown argument (%s) passed to filterbank.",argv[i]);
        error_message(message);
      }
      i++;
    }
  }

  /* if no IF streams selected, flag all */
  if (strings_equal(ifstream,"XXXX")) strcpy(ifstream,"YYYY");

  /* adjust finish time if necessary */
//  if (final_time > 0.0) final_time+=start_time;

  if (!opened) {
    /* no output file selected, use standard output */
    output=stdout;
    strcpy(outfile,"stdout");
  }


  fileidx_start = fileidx;
//  strcpy(inpfile,argv[fileidx]);
//  input=open_file(inpfile,"rb");


  /* hard code some values for now */

  nchans=128;
  foff=838.860800/4/nchans * -1.0;        // XXX hack
  fcentral=1430.0;
  fch1=fcentral  - ((nchans/2)+0.5)*foff;    // A bit of a guess
  nbits=8;
//  tsamp=1.0/(double)sample_rate;
  if (tsamp == 0.0) tsamp = 0.000625;
  nifs=1;

if (obits == 32) nbits = 32;

  /* ------------------------------ */

//  inputdata=typeof_inputdata(input,inpfile);
//  fclose(input);

  /* main loop around input files */
  while (fileidx <= nfiles) {

    /* open up input file */
    strcpy(inpfile,argv[fileidx]);
    input=open_file(inpfile,"rb");

    /* bump past pcap global header */
    bump = fread(&payload, sizeof(char), 24, input);

    /* skip first part if so requested */
    if (start_time) {       
        sample_skip = (int) 11 * sample_rate * start_time; /* Number of iBOBs is 11*/
        fprintf(stderr, "Skipping %6.2f seconds (%d samples)\n",
                start_time, sample_skip);
        for (i=0;i<sample_skip;i++) {
            bump = fread(&payload, sizeof(char), 8, input); /* skip 12 bits */
            bump = fread(&packet_saved_len, sizeof(int), 1, input);           
            bump = fread(&packet_len, sizeof(int), 1, input);            
                 
                 if((packet_len == packet_saved_len) && (packet_len <= 1500)) {
                    //skip this packet using the usual maneuver
                     bump = fread(&payload, sizeof(char), packet_len, input);
                 } else {
                    //skip this packet using a more careful method
                     do {
                        bump = fread(&payload, sizeof(char), 16, input);

                        packet_len = (((unsigned int) payload[12])%256);
                        packet_len += (((unsigned int) payload[13])%256) * 256;
                        
                        packet_saved_len = (((unsigned int) payload[8])%256);
                        packet_saved_len += (((unsigned int) payload[9])%256) * 256;                        
                        fseek(input, -15, SEEK_CUR);
                    } while (packet_len != 575 || packet_saved_len !=575);
                        fseek(input, -1, SEEK_CUR);
                    

                }
   
        }
    }
    /* read the time stamp for the first packet */
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




/* if not overruled on command line */
    if (tstart==0.0 && first_packet) {
        /* 1-1-1970 to MJD */
        //printf("ts_sec is %d  sizeof int is %d\n", ts_sec, sizeof(int));
        tstart=40587.0+((double)ts_sec+((double)ts_usec)/1000000.0)/86400.0;
    }

    if (fileidx == fileidx_start) {
        /* add on a time offset in seconds to the start time */
        tstart+=time_offset/86400.0;
        /* broadcast the header */
          
       /* to the multiple output files */
		if(!statonly){
			for (i=0; i<11; i++) {
				for (j=0; j<4; j++) {
					/* open them first */
					sprintf(outfile,"%sibob%02d%c.fil",file_prepend, i, chan_name[j]);
	//                  fprintf(stderr, "going to open %d %d: %s \n", i, j, outfile); fflush(stderr);
					 mult_out[i][j]=fopen(outfile,"wb");
					 /* write header */
					 filterbank_header(mult_out[i][j]);
				 }
			 }
		}

        if (headeronly) exit(0);
    }
    

    /* equalizing? */
    if (equalize_all) {
        equalize_init = 1;
    } else {
        for (i=0; i<128; i++) {
            for (j=0; j<4; j++) {
                for (k=0; k<11; k++) {
                    equalize_coeff[i][k][j] = 1.0;
                }
            }
        }
    }

       
//    fprintf(stderr, "size=%d \n", sizeof(char));
    
    /* run over the spectra, headers, etc */
    while (bump != 0) {

        /* read pcap per-packet header */
//      bump = fread(&payload, sizeof(char), 16, input);
//      packet_len = (((unsigned int) payload[12])%256);
//      packet_len += (((unsigned int) payload[13])%256) * 256;
//      packet_len += (((unsigned int) payload[14])%256) * 65536;
//      packet_len += (((unsigned int) payload[15])%256) * 16777216;

        if (!first_packet) {
            bump = fread(&payload, sizeof(char), 16, input); /* skip 16 bytes */
         
         
            ts_sec = (((unsigned int) payload[0])%256);
            ts_sec += (((unsigned int) payload[1])%256) * 256;
            ts_sec += (((unsigned int) payload[2])%256) * 65536;
            ts_sec += (((unsigned int) payload[3])%256) * 16777216;             
   
            ts_usec = (((unsigned int) payload[4])%256);
            ts_usec += (((unsigned int) payload[5])%256) * 256;
            ts_usec += (((unsigned int) payload[6])%256) * 65536;
            ts_usec += (((unsigned int) payload[7])%256) * 16777216;
   
   
   
            packet_len = (((unsigned int) payload[12])%256);
            packet_len += (((unsigned int) payload[13])%256) * 256;

            packet_saved_len = (((unsigned int) payload[8])%256);
            packet_saved_len += (((unsigned int) payload[9])%256) * 256;


            //bump = fread(&packet_len, sizeof(int), 1, input);
            //printf("not on first packet: %d\n", packet_len);
        } else {
            /* for first packet 8 bits of time stamp have been read farther up */
            
            bump = fread(&payload, sizeof(char), 8, input); /* skip 4 bits */


            packet_len = (((unsigned int) payload[4])%256);
            packet_len += (((unsigned int) payload[5])%256) * 256;

            packet_saved_len = (((unsigned int) payload[0])%256);
            packet_saved_len += (((unsigned int) payload[1])%256) * 256;
           
            //bump = fread(&packet_len, sizeof(int), 1, input);
            first_packet = 0;
            //printf("on first packet: %d\n", packet_len);
        }

        /* process if correct size */
        if(packet_len == 575 && packet_saved_len == 575) {
            bump = fread(&payload, sizeof(char), packet_len, input);

//          fprintf(stderr, "bump = %d \n", bump); fflush(stderr);
            ip = (((unsigned int) payload[29])%256);
            


            
            /* accumulation number for this IP */
            if( ((int) (ip-38) >= 0) &&  ((int) (ip-38) < 12)) {
            
            accumulation_number[ip-38] = (((unsigned int) payload[46])%256) +
                (((unsigned int) payload[47])%256) * 256   + 
                (((unsigned int) payload[48])%256) * 65536 +
                (((unsigned int) payload[49])%256) * 16777216;

          //fprintf(stderr, "ip-38=%d acc_no=%d counter=%d packet_len=%d, packet_saved_len=%d\n", ip-38, 
          //        (int) accumulation_number[ip-38], (int) counter[ip-38], packet_len, packet_saved_len); fflush(stderr);


  if(start_times[ip-38] == -1) start_times[ip-38]=(double) ts_sec + ((double) ts_usec/1000000.);
  if(start_cnt[ip-38] == -1) start_cnt[ip-38] = (double) accumulation_number[ip-38];
  
  end_times[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
  end_cnt[ip-38] = (double) accumulation_number[ip-38];

 
            /* only start after the first counter reset */
            if (start_at_zero) {
                if ((oldaccumulation[ip-38]==-1) && (accumulation_number[ip-38]==0)) {
                    start_at_zero=0;
                    fprintf(stderr, "Started output at first zero counter\n");  fflush(stderr);
                } else {
                    continue;
                }
            }

               

            /* quit if this packet exceeds the max read time */
            if (final_time) {       
                if (tsamp*counter[ip-38] >= final_time) {
                    fprintf(stderr, "Quitting after %fs\n",
                            (float) tsamp*counter[ip-38]);
                    break;
                }
            }
                    
            /* missing packets? */
            
            delta = accumulation_number[ip-38] - oldaccumulation[ip-38];
            if ( (accumulation_number[ip-38] - oldaccumulation[ip-38] != 1) &&
                 (accumulation_number[ip-38] - oldaccumulation[ip-38] != -1600) &&
                 (oldaccumulation[ip-38] != -1) ) {
                
                fprintf(stderr, 
                        "Error in packet order: ibob=%d at counter=%d, time=%6.1f (new=%d  old=%d delta=%d)\n",
                        ip-38, (int) counter[ip-38], (float) tsamp*counter[ip-38],
                        accumulation_number[ip-38], oldaccumulation[ip-38], 
                        (accumulation_number[ip-38] - oldaccumulation[ip-38]));
						if( delta > 2500 ) printf("absurdly high delta: %d\n", delta);
                        if ( delta > 0  && delta < 2500) {
								if( (correcterror == 1) ) 
								{
									   fprintf(stderr, "Correcting...%d\n",  (accumulation_number[ip-38] - oldaccumulation[ip-38]) );
									   
									   //reset time counter
									   tacc[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
									   
									  for(k = 1;k < (accumulation_number[ip-38] - oldaccumulation[ip-38]);k++){ 
											if(!statonly){		
											   if (obits==32) {
												for (j=0; j<4; j++) { /* for all channels in this iBOB */
														if (! reversed) {
															/* normal order */
															for (i=64 ;i < 128;i++) {
																temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
																fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
															}
															for (i=0  ;i < 64 ;i++) {
																temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
																fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]); 
															}
														} else {
															/* --------- output reversed band ---------------*/
															for (i=63 ;i >= 0 ;i--) {
																temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
																fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
															}
															for (i=127;i >= 64;i--) {
																temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
																fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
															}
														}
												}
											
											} else {
								
												for (j=0; j<4; j++) { /* for all channels in this iBOB */
														if (! reversed) {
															/* normal order */
															for (i=64 ;i < 128;i++) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
															for (i=0  ;i < 64 ;i++) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
														} else {
															/* --------- output reversed band ---------------*/
															for (i=63 ;i >= 0 ;i--) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
															for (i=127;i >= 64;i--) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
														}
												}
								
											}                
									  	  }
									  }
								}            
						
    					}
    	    } else {



                    if( tacc[ip-38] == 0.0 ) {
                       tacc[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
				    } else {
				       //printf("%f %f\n", (double) ts_sec + ((double) ts_usec/1000000.), tacc[ip-38]);
				       tdelta[ip-38] += ( ((double) ts_sec + ((double) ts_usec/1000000.)) - tacc[ip-38] );				    
					    delta_counter[ip-38]++;
				       tacc[ip-38] = (double) ts_sec + (((double) ts_usec)/1000000.);
					}


					if(skip_packet_2){
						 if( (accumulation_number[ip-38] == 2) && (accumulation_number[ip-38] - oldaccumulation[ip-38] == 1) ){											 	
							 /* first time through */
							 if(skip_old_acc[ip-38] == 0.0) {
								 skip_old_acc[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
							 } else {
							 	 /* catch dropped 2nd packet */
							 	 if (  (((double) ts_sec + ((double) ts_usec/1000000.)) - skip_old_acc[ip-38] )  < 1.2){
									  skip_acc[ip-38] = skip_acc[ip-38] +  (((double) ts_sec + ((double) ts_usec/1000000.)) - skip_old_acc[ip-38] ) ;
									  skip_old_acc[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
									  skip_count[ip-38] = skip_count[ip-38] + 1;	
								 } else {
								 	  skip_old_acc[ip-38] = (double) ts_sec + ((double) ts_usec/1000000.);
								 }
							 }
						 }	
					}



            }
            oldaccumulation[ip-38] = accumulation_number[ip-38];

            /* what's the order like? */
            if (check_order && (counter[ip-38] == prev_counter - 1)) { /* was: -1 */
                fprintf(stderr, "\nOops, ibob out of order: ibob=%d:   counter=%d, prev_counter=%d\n",
                        ip-38, (int) counter[ip-38], (int) prev_counter);
                fprintf(stderr, "  Last counters:\n");
                for (i=0;i<11;i++) {
                    fprintf(stderr, "    ibob %d: %d\n", i, (int) last_counter[i]);
                }
                exit(-1);
            }
            prev_counter = counter[ip-38];
            last_counter[ip-38] = counter[ip-38];
            counter[ip-38]++;
//          fprintf(stderr, "ip-38=%d counter=%d\n ", ip-38, (int) counter[ip-38]); fflush(stderr);

            /* occasionally print progress */
            if ( (ip-38==0) && (counter[ip-38]+1)%1000==1) { 
                fprintf(stderr, "\r Writing packet=%d time=%6.1f", 
                        (int) counter[ip-38], (float) tsamp*counter[ip-38]);
            }

			/* normal order */
			
			/* investigate off-by-one */
			
			/*
			if((ip-38 == 2) && (accumulation_number[ip-38]==0 || accumulation_number[ip-38]==1600 || accumulation_number[ip-38]==1599 || accumulation_number[ip-38]==1)){
				spectra_sum=0;
				for (i=0  ;i < 128 ;i++) {
				   printf("%d ", (unsigned int) payload[(63+(i*4)+j)] );
				   spectra_sum += (unsigned int) payload[(63+(i*4)+j)];

				}
				printf("\nibob: %i acc: %i sum: %f\n", ip-38, accumulation_number[ip-38], spectra_sum);
			}
			*/
			
            /* skip packet number 2 until 1600/1601 packets per second problem is fixed */
            if (skip_packet_2 && accumulation_number[ip-38]==2) continue;

            /* -- what to do with packet? read for equalization, for summing or for writing -- */
            
  

            if (obits==32) {
                for (j=0; j<4; j++) { /* for all channels in this iBOB */
                        /* ibob sends channels 64-127 first, then 0-63 */

                      //for (i=63 ;i >= 0 ;i--) fprintf(stderr, "%d==%i ", i, (signed int) payload[(63+(i*4)+j)]);
                      //for (i=127;i >= 64;i--) fprintf(stderr, "%d==%u ", i, (unsigned int) payload[(63+(i*4)+j)]);
                      //fprintf(stderr, "\n"); fflush(stderr);
					if(!statonly) {
                        if (! reversed) {
                            /* normal order */
                            for (i=64 ;i < 128;i++) {
                                temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
                                fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
//                                fprintf(stderr, "%d==%f ", i, temp_power);
                            }
                            for (i=0  ;i < 64 ;i++) {
                                temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
                                fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]); 
                            }
                        } else {
                            /* --------- output reversed band ---------------*/
                            for (i=63 ;i >= 0 ;i--) {
                                temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
                                fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
//                                fprintf(stderr, "%d==%f ", i, temp_power);
                            }
                            for (i=127;i >= 64;i--) {
                                temp_power = ((float) payload[(63+(i*4)+j)]) * equalize_coeff[i][ip-38][j];
                                fwrite(&temp_power,sizeof(float),1,mult_out[ip-38][j]);
                            }
                        }
                	}
                }
            
            } else {
				if(!statonly) {
					  for (j=0; j<4; j++) { /* for all channels in this iBOB */
							  /* ibob sends channels 64-127 first, then 0-63 */
							
							//for (i=63 ;i >= 0 ;i--) fprintf(stderr, "%d==%i ", i, (signed int) payload[(63+(i*4)+j)]);
							//for (i=127;i >= 64;i--) fprintf(stderr, "%d==%u ", i, (unsigned int) payload[(63+(i*4)+j)]);
							//fprintf(stderr, "\n"); fflush(stderr);
							  if (! reversed) {
								  /* normal order */
								  for (i=64 ;i < 128;i++) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
								  for (i=0  ;i < 64 ;i++) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
							  } else {
								  /* --------- output reversed band ---------------*/
								  for (i=63 ;i >= 0 ;i--) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
								  for (i=127;i >= 64;i--) fwrite(&payload[(63+(i*4)+j)],sizeof(char),1,mult_out[ip-38][j]);
							  }
					  }
				}
	
            }
                                         
            }
            
        } else {
            /* wrong packet length ??   */
            /* not sure about this code */
            if(bump != 0) {
                if((packet_len == packet_saved_len) && (packet_len <= 1500) && (packet_len != 0)) {
                //skip this packet using the usual maneuver
                     fprintf(stderr,"skipping packet with len=%d  saved_len=%d at ip=%d\n", packet_len, packet_saved_len, ip);
                     fflush(stderr);
                     bump = fread(&payload, sizeof(char), packet_len, input);
                     fprintf(stderr, "done..\n");
                     fflush(stderr);
                } else {
                //skip this packet using a more careful method
                     fprintf(stderr,"Messed up packet! skipping packet with len=%d  saved_len=%d at ip=%d\n", packet_len, packet_saved_len, ip);
                     do {
                        bump = fread(&payload, sizeof(char), 16, input);

                        packet_len = (((unsigned int) payload[12])%256);
                        packet_len += (((unsigned int) payload[13])%256) * 256;
                        
                        packet_saved_len = (((unsigned int) payload[8])%256);
                        packet_saved_len += (((unsigned int) payload[9])%256) * 256;                        
                        fseek(input, -15, SEEK_CUR);
                    } while (packet_len != 575 || packet_saved_len !=575);
                        fseek(input, -1, SEEK_CUR);
                    

                }
		    }



        }
    }

    fileidx++;
    fclose(input);
  }

  /* all done, update log, close all files and exit normally */
  if(!statonly) { for (i=0; i<11; i++) for (j=0; j<4; j++) close(mult_out[i][j]); }
  printf("\n");
  for (i=0; i<11; i++){
    printf("IBOB: %d - %10Lf pkts received in %6.12Lf sec. - avg acc time: %6.12Lf sec. - avg deltaT:  %6.12Lf sec. \n", i, (long double) (end_cnt[i]-start_cnt[i]), (long double) (end_times[i] - start_times[i]),  (long double) ((end_times[i] - start_times[i]) / (end_cnt[i]-start_cnt[i])), (long double) tdelta[i]/delta_counter[i]   );
	average = average +  (tdelta[i]/delta_counter[i]);
  }
    
  printf("\n\ntotal average: %6.12Lf sec on mjd: %f\n", (long double) (average/11.0), tstart);
  
  if(skip_packet_2) {
	  	for (i=0; i<11; i++){
  			printf("\n\naverage between 1pps: %6.12Lf sec on mjd: %f\n", (long double) (skip_acc[i]/skip_count[i]), tstart);
  		}
  }

  
  update_log("finished");
  close_log();
  /*fclose(output);*/
  fprintf(stderr,"Completed Successfully\n");
  exit(0);
}
