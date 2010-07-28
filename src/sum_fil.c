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
#include "sigproc.h"
#include <string.h>

int wapp_header_size, wapp_incfile_length;
int nbins;
	double period;

main(int argc, char *argv[]) 
{
	FILE *fileptr[50];
	char filename[80],*telescope,*backend,*datatype,message[80],unit[16];
	int i,j,year,month,day,check,rah,ram,ded,dem;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;

	int writeobsdbline;


  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i])) {
        nfiles++;
        i++;
  }
  if (!nfiles) {
	  error_message("no input files supplied on command line!");
	  exit(1);
  }

  if (nfiles > 50) {
      error_message("too many input files supplied on command line (max 50)!");
	  exit(1);
  }





	fileptr=stdin;
	strcpy(filename,"stdin");
	strcpy(rawdatafile,"stdin");
	pulsarcentric=barycentric=0;

	writeobsdbline=0;

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