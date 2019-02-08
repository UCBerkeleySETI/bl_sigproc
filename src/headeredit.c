#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
   header - show header parameters in a data file
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
	
	
void headeredit_help() /*includefile*/
{
  puts("");
  puts("headeredit  - modify header parameters of filterbank data\n");
  puts("usage: headeredit {filename} -{options}\n");
  puts("filename is the filterbank data file \n");
  puts("options:\n");
  puts("-fch1 <double>       - set frequency of channel 1 in MHz");
  puts("-foff <double>       - set channel bandwidth in MHz");
  puts("-nchans <int>        - set number of channels");
  puts("-tstart <double>     - set time stamp of first sample (MJD)");
  puts("-tsamp <double>      - set sample time (us)");
  puts("-nbits <int>         - set number of bits per sample");
  puts("-src_raj <double>    - set RA (J2000, HHMMSS.ssss)");
  puts("-src_dej <double>    - set DEC (J2000, DDMMSS.ssss)");  
  puts("");
}
main(int argc, char *argv[]) 
{
	FILE *fileptr;
	char filename[80],*telescope,*backend,*datatype,message[80],unit[16];
	int i=0,j=0;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;
	char sname[80];

	if (argc>1) {
		if (help_required(argv[1])) {
			headeredit_help();
			exit(0);
		} else if (file_exists(argv[1])) {
			strcpy(filename,argv[1]);
			fileptr=open_file(filename,"r+b");
		} else if (!file_exists(argv[1]) && (strncmp(argv[1],"-",1) !=0)) {
			sprintf(message,"Data file: %s not found...\n",argv[1]);
			error_message(message);
			exit(1);
		}
	}




i=2;

	if (argc>2) {


		/* check command-line parameters */ 
		while (i<argc) {
		
			if (strings_equal(argv[i],"-fch1")) {
				fch1=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "fch1")) fwrite(&fch1,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-foff")) {
				foff=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "foff")) fwrite(&foff,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-nchans")) {
				nchans = (int) strtol(argv[++i], NULL, 10);
				if (move_to_keyword(fileptr, "nchans")) fwrite(&nchans,sizeof(int),1,fileptr);  

			} else if (strings_equal(argv[i],"-tstart")) {
				tstart=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "tstart")) fwrite(&tstart,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-tsamp")) {
				tsamp=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "tsamp")) fwrite(&tsamp,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-nbits")) {
				nbits = (int) strtol(argv[++i], NULL, 10);
				if (move_to_keyword(fileptr, "nbits")) fwrite(&nbits,sizeof(int),1,fileptr);  

			} else if (strings_equal(argv[i],"-src_raj")) {
				src_raj=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "src_raj")) fwrite(&src_raj,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-src_dej")) {
				src_dej=strtod(argv[++i], NULL);
				if (move_to_keyword(fileptr, "src_dej")) fwrite(&src_dej,sizeof(double),1,fileptr);  

			} else if (strings_equal(argv[i],"-source_name")) {
				strcpy(sname,argv[++i]);
				if (move_to_keyword(fileptr, "source_name")) fwrite(&sname,sizeof(char),14,fileptr); // Fixed for source_name=FindMe			 
			} else {
				sprintf(message,"unknown argument (%s) passed to header",argv[i]);
				error_message(message);
			}
			i++;
		}



		exit(0);
	}


}


void get_string(FILE *inputfile, int *nbytes, char string[])
{
  int nchar;
  strcpy(string,"ERROR");
  fread(&nchar, sizeof(int), 1, inputfile);
  *nbytes=sizeof(int);
  if (feof(inputfile)) {
  	fprintf(stderr,"Error: reached end of file...\n");
  	exit(0);
  }
  if (nchar>80 || nchar<1) return;
  fread(string, nchar, 1, inputfile);
  string[nchar]='\0';
  *nbytes+=nchar;
}



/* attempt to scan to a keyword in a filterbank file */
int move_to_keyword(FILE *inputfile, char *keyword) /* includefile */
{
  int nbytes,totalbytes;
  char string[80];
  rewind(inputfile);

  
  /* try to read in the first line of the header */
  get_string(inputfile,&nbytes,string);
  if (!strings_equal(string,"HEADER_START")) {
  	  fprintf(stderr,"This file is not in filterbank format...\n");
	/* the data file is not in standard format, rewind and return */
	rewind(inputfile);
	return 0;
  }
  /* store total number of bytes read so far */
  totalbytes=nbytes;

  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string(inputfile,&nbytes,string);
    if (strings_equal(string,"HEADER_END")) {
    	fprintf(stderr,"Couldn't find keyword: %s\n", keyword);
		
		rewind(inputfile);
		return 0;
	}
    if (strings_equal(string,keyword)) {
        fprintf(stderr,"Moving to keyword: %s\n", keyword);

      return totalbytes + nbytes;
    } 

    totalbytes+=1;

    //rewind to 1 byte ahead of previous read, try again
    fseek(inputfile, -1 * (nbytes - 1), SEEK_CUR);    
	if (totalbytes > 1000) {
    	fprintf(stderr,"Couldn't find keyword: %s in 1000 bytes\n", keyword);
		rewind(inputfile);
		return 0;
    } 

  } 

}





