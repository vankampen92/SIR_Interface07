/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                           POWER SPECTRUM ANALYSIS                         */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "pwsp.h"

void Saving_to_File_powerSpc(char Prefix[], float Time[], float Number[], 
			     int Points, int File_Number)
{
  int i;
  char file[12]; 
  FILE *fp;
  
  file[0]='\0';    
  fitxer(file, Prefix, File_Number, ".dat");
  fp = fopen(file, "w");
  for(i=1; i<=Points; i++)
    fprintf(fp,"%f\t %f\n",Time[i], Number[i]);
  fclose(fp);
}

void Saving_to_File_autocorr(char Prefix[], float Number[], 
			     int Points, int File_Number)
{
  int i;
  char file[12]; 
  FILE *fp;
  
  file[0]='\0';    
  fitxer(file, Prefix, File_Number, ".dat");
  fp = fopen(file, "w");
  for(i=0; i<Points; i++)
    fprintf(fp,"%d\t %f\n",i, Number[i+1]);
  fclose(fp);
}

void Showing_powerSpc(float Time[], float Number[], int Points)
{
  int i;
  for(i=1; i<=Points; i++)
    printf("%f\t %f\n",Time[i], Number[i]);
}














