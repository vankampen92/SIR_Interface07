/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*	                                                                     */
/*                          (CONSTANT COMMUNITY SIZE)                        */
/*                                                                           */
/*                             David Alonso, 2000 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "SIR.h"

void ShowVillageContent(Community *pVil, int i)
{ 
    /* */
  int k;
  printf(" Village No %d\n", i);
  printf(" Number of Inhabitants\t\t%d\n",pVil->N);
  printf(" Total Number of Infected Individuals:\t\t%d\n",pVil->n);
  printf(" Total number of Susceptival Individuals:\t\t%d\n",pVil->m);
  printf(" Number of Infected Individuals for each Strain:\n");
  printf(" I[%d] = %d\n",1,pVil->I[1]);
  printf(" Patch position: \t\t\t");
  printf("(%f, %f)\n", pVil->position.x, pVil->position.y);
  getchar();
}

void ShowVillage_StrainEvolution(Community *pVil, int i, double time)
{

  printf(" Time: %g\n",time); 
  printf(" Village Label %d: \t",i);
  printf("%d  ",pVil->I[0]); 
  printf("%d  ",pVil->I[1]);
  printf("%d  ",pVil->I[2]);
  printf("\n");

}
 
void SaveTimeSeries(Community *pVil, int i, double te)
{
  /* Data are saved in POPULATION fraction */
  FILE *pf;
  char file[14];
  float fraction;
  
  /*Susceptible */
  file[0]='\0';   
  fitxer(file, "sus", i, ".dat");
  pf = fopen(file,"a");
  fraction = pVil->I[0]/(float)pVil->N;
  fprintf(pf,"%f\n", fraction); 
  fclose(pf); 

  /* Infective */
  file[0]='\0';   
  fitxer(file, "inf", i, ".dat");
  pf = fopen(file,"a");
  fraction = pVil->I[1]/(float)pVil->N;
  fprintf(pf,"%f\n", fraction);
  fclose(pf);
}

void write_vectors(float *time, float *x, float *y, int nP)
{
  int i;
  
  for(i = 1; i<=nP; i++){
    printf("%f\t%f\t%f\n",time[i],x[i],y[i]);
  }
  //getchar();
}

void SaveVillageContent(char Pre[], Community *pVil, 
			int k_Type, int i, double te)
{ 
  /* Population Fraction are saved. 
     Since some files are openen in "a" append mode, make sure
     you erase before starting a new simulation */

  FILE *pf;
  char file[14];
  int k,j, n_0, n_1;
  float fraction;
  /* */
  
  file[0]='\0';   
  fitxer(file, Pre, i, ".dat");
  pf = fopen(file,"a");
  n_0 = pVil->I[0];
  n_1 = pVil->I[1];

  if(k_Type == 0){
    fprintf(pf,"%g\t%d\t", te, pVil->I[0]);
    fprintf(pf,"%d\t", pVil->I[1]);
    fprintf(pf,"\n");
  }
  else{
    fraction = pVil->I[0]/(float)pVil->N;
    fprintf(pf,"%g\t%f\t", te, fraction);
    fraction = pVil->I[1]/(float)pVil->N;
    fprintf(pf,"%f\t", fraction);
    fprintf(pf,"\n");
  }
  fclose(pf); 
  /* getchar(); */
}
 





