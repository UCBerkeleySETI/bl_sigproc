#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
   dedisperse  - dedisperses raw filterbank data or folded pulse profiles
 */
#include "dedisperse.h"
#include <string.h>
FILE *input, *output;
char  inpfile[80], outfile[80], ignfile[80];

/* global variables describing the operating mode */
int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
double refrf,userdm,fcorrect;
float clipvalue,jyf1,jyf2;
int fftshift;
int doSubZero;
/**
 * doInputClip is a series of binary flags!
 * Clip below specified: 00000001 1
 * Clip above specified: 00000010 2
 * Clip below median   : 00000100 4
 * Clip above median   : 00001000 8
 * 
 * Clearly one can have above and below (3 or 15, or some other bizare combination!).
 * M Keith 2007
 **/
int doInputClip;
FILE *inputClipLowFile, *inputClipHighFile;
float *inputClipLow,*inputClipHigh;


#include "wapp_header.h"
#include "key.h"
struct WAPP_HEADER *wapp;
struct WAPP_HEADER head;

main (int argc, char *argv[])
{
	/* local variables */
	char string[80];
	int i,j,useroutput=0,autooutput=0,nfiles,fileidx,sigproc,scan_number,subscan=0;
	doBaselineFile= 0;

	/* check number of command-line arguments and print help if necessary */
	if (argc<2) {
		dedisperse_help();
		exit(0);
	} else {
		print_version(argv[0],argv[1]);
	}

	/* print help if necessary */
	if (help_required(argv[1])) {
		dedisperse_help();
		exit(0);
	}
	/* work out how many files are on the command line */
	i=1;
	nfiles=0;
	while(file_exists(argv[i])) {
		nfiles++;
		i++;
	}

	/* set up default globals */
	userbins=usrdm=asciipol=stream=clipping=swapout=headerless=0;
	sumifs=wapp_inv=wapp_off=barycentric=0;
	nobits=32;
	ascii=1;
	fftshift=1;
	profnum1 = 0;
	profnum2 = 1000;
	nbands=baseline=1;
	clipvalue=refrf=userdm=fcorrect=0.0;
	jyfactor=jyf1=jyf2=1.0;
	refdm=-1.0;
	doInputClip=doSubZero=0;
	output=NULL;
	strcpy(ignfile,"");

	/* now parse any remaining command line parameters */
	if (argc>nfiles) {
		i=nfiles+1;
		while (i<argc) {
			if (strings_equal(argv[i],"-d")) {
				/* get DM from the user */
				userdm=atof(argv[++i]);
				usrdm=1;
			} else if (strings_equal(argv[i],"-b")) {
				if (strings_equal(argv[++i],"nchans")) {
					nbands=0;
				} else {
					/* get number of subbands */
					nbands=atoi(argv[i]);
				}
			} else if (strings_equal(argv[i],"-B")) {
				/* get output number of bits */
				nobits=atoi(argv[++i]);
			} else if (strings_equal(argv[i],"-o")) {
				/* get and open file for output */
				strcpy(outfile,argv[++i]);
				useroutput=1;
			} else if (strings_equal(argv[i],"-autooutput")) {
				/* open output files automatically as mjd.frq.scn */
				autooutput=1;
			} else if (strings_equal(argv[i],"-c")) {
				/* set clipvalue */
				clipvalue=atof(argv[++i]);
				baseline=clipping=1;
			} else if (strings_equal(argv[i],"-j")) {
				/* get user-supplied Jansky calibration factor */
				jyfactor=atof(argv[++i]);
			} else if (strings_equal(argv[i],"-J")) {
				/* get user-supplied Jansky calibration factors */
				jyf1=atof(argv[++i]);
				jyf2=atof(argv[++i]);
			} else if (strings_equal(argv[i],"-old")) {
				/* revert to the old method of profile dedispersion ! */
				fftshift=0;
			} else if (strings_equal(argv[i],"-f")) {
				/* get optional reference frequency for dedispersion */
				if (strings_equal(argv[++i],"mid")) {
					/* will set reference frequency to band centre */
					refrf = -1.0; 
				} else {
					/* will set user-supplied reference frequency */
					refrf = atof(argv[i]);
				}
			} else if (strings_equal(argv[i],"-F")) {
				/* get optional frequency to correct value in header */
				fcorrect = atof(argv[++i]);
			} else if (strings_equal(argv[i],"-n")) {
				/* get number of bins */
				userbins=atoi(argv[++i]);
			} else if (strings_equal(argv[i],"-i")) {
				/* file containing channel numbers to ignore */
				strcpy(ignfile,argv[++i]);
			} else if (strings_equal(argv[i],"-p")) {
				/* get which profile numbers to use - starts at 0 */
				profnum1=atoi(argv[++i]);
				profnum2=atoi(argv[++i]);
			} else if (strings_equal(argv[i],"-wappinvert")) {
				/* Make wapp frequency channel increments negative */
				wapp_inv=1;
			} else if (strings_equal(argv[i],"-wappoffset")) {
				/* Puts fmid between two middle channels; otherwise it's in the
				 * 	middle of one of those channels (for use with 
				 * 	pre-52900 data ONLY) */
				wapp_off=1;
			} else if (strings_equal(argv[i],"-icliplow")) {
				/* set clipvalue */
				inputClipLowFile=fopen(argv[++i],"r");
				if(!inputClipLowFile){
					fprintf(stderr,"ERROR: requested low clip file %s does not exist.\n",argv[i-1]);
					exit(1);	
				}

				doInputClip = doInputClip | 1; /* set bit 1 */
			} else if (strings_equal(argv[i],"-icliphigh")) {
				/* set clipvalue */
				inputClipHighFile=fopen(argv[++i],"r");
				if(!inputClipHighFile){
					fprintf(stderr,"ERROR: requested high clip file %s does not exist.\n",argv[i-1]);
					exit(1);        
				}

				doInputClip = doInputClip | 2; /* set bit 2 */
			} else if (strings_equal(argv[i],"-icliplow-median")) {
                                /* set clipvalue */
                                doInputClip = doInputClip | 4; /* set bit 3 */
                        } else if (strings_equal(argv[i],"-icliphigh-median")) {
                                /* set clipvalue */
                                doInputClip = doInputClip | 8; /* set bit 4 */
                        } else if (strings_equal(argv[i],"-swapout")) {
				/* perform byte swapping on all output data */
				swapout=1;
			} else if (strings_equal(argv[i],"-headerless")) {
				/* no headers */
				headerless=1;
			} else if (strings_equal(argv[i],"-nobaseline")) {
				/* no baseline subtraction */
				baseline=0;
			} else if (strings_equal(argv[i],"-baseline-file")) {
				/* use an external file for scaling channels */
				strcpy(baselineFile,argv[++i]);
				doBaselineFile = 1;
				baselineValidStart=baselineValidEnd=0;
			} else if (strings_equal(argv[i],"-subzero")) {
				/* do the Zero DM subtraction */
				doSubZero = 1;
			} else if (strings_equal(argv[i],"-sumifs")) {
				/* Sum IFs for final profile */
				sumifs = 1;;
			} else if (strings_equal(argv[i],"-epn")) {
				/* write EPN profiles when dedispersing folded data */
				ascii=0;
			} else if (strings_equal(argv[i],"-ascii")) {
				/* write ascii profiles when dedispersing folded data */
				ascii=1;
			} else if (strings_equal(argv[i],"-asciipol")) {
				/* write data as ASCII numbers for Jim's polarization code */
				asciipol=1;
			} else if (strings_equal(argv[i],"-stream")) {
				/* write data as ASCII streams */
				stream=1;
			} else {
				/* unknown argument passed down - stop! */
				dedisperse_help();
				sprintf(string,"unknown argument (%s) passed to %s",argv[i],argv[0]);
				error_message(string);
			}
			i++;
		}
	}

	if (!useroutput && !autooutput) {
		/* no output file selected, use standard output */
		output=stdout;
		strcpy(outfile,"stdout");
	}

	if (!nfiles) {
		strcpy(inpfile,"stdin");
		nfiles=1;
	}
	fileidx=1;


	/* main loop around input files */
	while (fileidx <= nfiles) {

		/* open input datafile and get type of file from the input data */
		if (strings_equal(inpfile,"stdin")) {
			input=stdin;
		} else {
			strcpy(inpfile,argv[fileidx]);
			input=open_file(inpfile,"rb");
		}

		if (autooutput) {
			/* get scan number from the file name's last four digits */
			scan_number=atoi(strtok(&inpfile[strlen(inpfile)-7],"."));
		}

		/* read in the header to establish what the input data are... */
		while (!feof(input)) {
			sigproc=read_header(input);
			if (sigproc) {
				if (data_type == 3) {
					if (autooutput) {
						subscan++;
						sprintf(outfile,"%.0f.%.0f.%03d%03d.prf",floor(tstart),fch1+foff*((float)nchans)/2.0,scan_number,subscan);
						output=open_file(outfile,"w");
					} else {
						if (output != stdout && fileidx == 1) output=open_file(outfile,"w");
					}
					machine2prf(input,output);
				} else {
					if (foff > 0.0) 
						error_message("dedisperse can't handle low->high frequency ordering!");
					if (fileidx == 1) {
						/*
						 * Need to read the list of values to clip input data on, if we are doing that.
						 */
						if(doInputClip & 1){ /* Low Clip */
							inputClipLow = (float*)malloc(nchans * sizeof(float));
							for(j = 0; j < nchans; j++){
								if(feof(inputClipLowFile)){
									fprintf(stderr,"Warning: Low clip file only lists %d channels, when there are %d in the input\n",j,nchans);
									inputClipLow[j] = 0;
									continue;
								}

								fscanf(inputClipLowFile,"%f\n",inputClipLow+j);
							}
							fclose(inputClipLowFile);
						}

						if(doInputClip & 2){ /* High Clip */
                                                        inputClipHigh = (float*)malloc(nchans * sizeof(float));
                                                        for(j = 0; j < nchans; j++){
								if(feof(inputClipHighFile)){
									fprintf(stderr,"Warning: High clip file only lists %d channels, when there are %d in the input\n",j,nchans);
									inputClipHigh[j] = 0;
									continue;
								}
                                                                fscanf(inputClipHighFile,"%f\n",inputClipHigh+j);
                                                        }
							fclose(inputClipHighFile);
                                                }
						if(doInputClip & 4){ /* Low Median Clip */
							inputClipLow = (float*)malloc(nchans * sizeof(float));
						}
						if(doInputClip & 8){ /* High Median Clip */
                                                        inputClipHigh = (float*)malloc(nchans * sizeof(float));
                                                }

						/* this is filterbank data */
						if (output!=stdout) output=open_file(outfile,"wb");
						/* open up logfile */
						open_log("dedisperse.monitor");
						/* broadcast header */
						dedisperse_header();
					}
					/* dedisperse */
					dedisperse_data(input,output);
					/* close log files if on last input file */
					if (fileidx == nfiles) {
						update_log("finished");
						close_log();
						/* we need to free any clip arrays that are in use! */
						if(doInputClip & 5){ /* Any Low Clip */
							free(inputClipLow);
                                                }

                                                if(doInputClip & 10){ /* Any High Clip */
							free(inputClipHigh);
                                                }

					}
				}
			} else {
				/* this is something else - possibly folded data - check it out.... */
				switch(typeof_inputdata(input,argv[1])) {
					case 2:
					case 4:
					case 6:
						/* SIGPROC/PSPM/WAPP/BPP timing-mode data */
						if (output != stdout && fileidx == 1) output=open_file(outfile,"w");
						machine2prf(input,output);
						exit(0);
						break;
					default:
						error_message("input data file is of unknown origin!!!");
				}
			}
			if (autooutput) fclose(output);
		}
		fileidx++;
		fclose(input);
	}
	exit(0);
}
