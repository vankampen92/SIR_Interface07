/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                                  SIR MODEL                                */
/*                      Computing the Analytic Power Spectrum                */
/*                             David Alonso, 2005 (c)                        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../../../SIR.h"
#include "../../../SIR_Analytic_General.h"

int No_of_Points;
double a,d,g,mu,b,Imm;
double d_0,d_1;
double b_0,b_1;
int POPULATION;
int COHERENCE;
int NUMBER_of_PLOTS;
float r;

void eps_WGW_image_coherence(char *File_Name, float **Data, int width, int length, float Max_z, float Min_z);
void vertical_Scale_coherence(float **Data, int NPTS, float *Max, float *Min);

int main(int argc, char **argv)
{
  double Psi, Fi; /* Populations Fractions */
  double d_Rel, b_Rel, value, coherence, relative_coherence, f_peak, f_semi;
  int point_Damping, point_Stochas, point_Stable;
  int i,j,k;  
  ParamSet Parameters;
  FILE *fp, *fp_1, *fp_2, *F_TEX;
  float **Data_1;
  float **Data_2;
  float Max_z, Min_z;
  double imm[5] = {0.000002,0.00001,0.001,0.01,0.1};
  char name[30];
  char pre_1[5] = "over";
  char pre_2[5] = "rel.";
  char suf[5] = ".eps";

  /* Initial settings and default values * * * * * * * * * * * * * * * * * * *  */
  No_of_Points =  10;  NUMBER_of_PLOTS = 1;
  POPULATION = 100000; COHERENCE = 0; 
  b = 1.175; a = 0.; d = 5.5e-5; g = 1./13.; mu = 0.; /* Measles as a default */
  b_0 = 0.; b_1 = 30.; d_0 = 0.; d_1 = 1.;            /* Range of Dimensionless values 
						         to Scan */
  Imm = 1.e-5;  /* External Transmission (Immigration), where Imm is given in days^(-1) */
  r = 0.1;
  /* END (Initial settings and default values) * * * * * * * * * * * * * ** * * */
  
  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);
  imm[0] = Imm;

  settingParameterStruct(&Parameters);
  modelReport("report.txt");

  Data_1 = matrix(0,No_of_Points, 0,No_of_Points);
  Data_2 = matrix(0,No_of_Points, 0,No_of_Points);
  
  Fixed_Points(&Parameters, &Fi, &Psi);  Stability(&Parameters); 
  Fixed_Points_General(&Parameters, &Fi, &Psi);  Stability_General(&Parameters);
 
  F_TEX = fopen("./Fig_eps/test_Fig.tex", "w"); 
  printf("Make sure the file ""./Fig_eps/test_Fig.tex"" exits!!!\n");
  printf("...otherwise you get a beautiful ""segmentation fault""\n\n");
  preambul_TEX_File(F_TEX);

  for(k=0; k<NUMBER_of_PLOTS; k++){
    if(k == 0) {
      fp = fopen("overall.dat", "w");      /* Overall amplification */
      fp_1 = fopen("coherence.dat", "w");  /* Coherence around the peak */
      fp_2 = fopen("rel_coh.dat", "w");    /* File (x,y,z) where z is the relative coherence */
    }
    Parameters.Imm = imm[k]; 
    Parameters.Beta = b; Parameters.Alpha = a; Parameters.Delta = d; Parameters.Gamma = g; Parameters.Mu = mu;  
    /* Making parameters Dimensionless quatities... */ 
    changingTimeScale(&Parameters, 1./g);
    printf("Dimensionless Parameaters: b/g = %g; g/g = %g; Imm/g = %g\n", 
	   Parameters.Beta, Parameters.Gamma, Parameters.Imm); 
    printf("Calculating the borders in the parameter space. \n");
    printf("NOTE: This calculation assumes dimensionless parameters.\n");
    
    for(i=0; i<No_of_Points; i++){
      for(j=0; j<No_of_Points; j++){
	
	point_Damping = point_Stochas = point_Stable = 0;
	
	d_Rel = d_0 + (i+1)* (d_1 - d_0)/(double)No_of_Points;
	b_Rel = b_0 + (j+1)* (b_1 - b_0)/(double)No_of_Points;
	
	re_settingParamStruct(&Parameters, d_Rel, b_Rel);
	//eigen_Values(&Parameters);
	/* 
	   Overall amplification is calculated only in 
	   the region where there is stability and damping.
	*/
	point_Stable = Condition_Stability(&Parameters);
	value = 0.;
	if(point_Stable == 1){
	  point_Damping = 2*Condition(&Parameters, 0.25);	
	  if(point_Damping == 2){
	    //point_Stochas = 4*Condition(&Parameters, 0.5);
	    point_Stochas = 4*Exact_Condition_Peak(&Parameters, 1);
	    //value = overall_Amplification(&Parameters, 1);
	    value = overall_Analytic(&Parameters, 1);
	  }
	}
	point_Stochas += point_Stable + point_Damping;
	/* 
	   The variable $value$ can actually have one out of four values:
	   0: Point is instable;
	   1: Point is stable;
	   3: Point is stable and present damping oscillations;
	   7: Point is stable and present both daming and stochastic amplification;
	*/
	if(COHERENCE == 1 || COHERENCE == 2){
	  relative_coherence = 0.;
	  if(point_Stochas == 7){
	    f_peak = Resonance_Frequency_Peak(&Parameters, 1);
	    f_semi = 2. * (double)r * f_peak;
	    relative_coherence = 100.*coherence_Analytic(&Parameters, 1, f_semi, f_peak)/value;
	  }
	  Data_2[i][j] = (float)relative_coherence;
	  if(k == 0){
	    fprintf(fp_1,"%g\t%g\t%g\n", d_Rel, b_Rel, coherence);
	    fprintf(fp_2,"%g\t%g\t%g\n", d_Rel, b_Rel, relative_coherence);
	  }
	}
	if(COHERENCE == 0 || COHERENCE == 2){
	  if(k == 0) fprintf(fp,"%g\t%g\t%g\n", d_Rel, b_Rel, value);
	  Data_1[i][j] = (float)value;
	}
      }
    }
    if(COHERENCE == 0 || COHERENCE == 2){ 
      vertical_Scale(Data_1, No_of_Points, &Max_z, &Min_z);
      //name_Ordered("OBar", k, ".eps", name);
      //eps_vertical_Bar(name, No_of_Points, No_of_Points, Max_z, Min_z);
      name_Ordered("OBar", k, ".dat", name);
      dat_vertical_Bar(name, No_of_Points, No_of_Points, Max_z, Min_z);

      name_Ordered(pre_1, k, suf, name);
      eps_WGW_image(name, Data_1, No_of_Points, No_of_Points, Max_z, Min_z);
      append_Fig_TEX(F_TEX, name, imm[k]);      
    }
    if(COHERENCE == 1 || COHERENCE == 2){
      vertical_Scale(Data_2, No_of_Points, &Max_z, &Min_z);
      //name_Ordered("CBar", k, suf, name);
      //eps_vertical_Bar(name, No_of_Points, No_of_Points, Max_z, Min_z);
      name_Ordered("CBar", k, ".dat", name);
      dat_vertical_Bar(name, No_of_Points, No_of_Points, Max_z, Min_z);
            
      name_Ordered(pre_2, k, suf, name);
      eps_WGW_image(name, Data_2, No_of_Points, No_of_Points, Max_z, Min_z);
      append_Fig_TEX(F_TEX, name, imm[k]);      
    }
    if(k == 0){
      fclose(fp);
      fclose(fp_1);
      fclose(fp_2);
    }
  }  
  
  close_TEX_File(F_TEX);
  free_matrix(Data_1, 0,No_of_Points, 0,No_of_Points);
  free_matrix(Data_2, 0,No_of_Points, 0,No_of_Points);						      
  modelReport("report_End.txt");
  printf("\nEnd of progam\n");
  return (0);
}

void eps_WGW_image_coherence(char *File_Name, float **Data, int width, int length, float Max_z, float Min_z)
{
  int k;
  double neww, newl, rat, xorigin, yorigin;
  FILE *FP;
  int i,j;
  float z, zero_offset;
  unsigned char grey;

  zero_offset = 0.05; 
  FP = fopen(File_Name,"w");
 
  /*  k = 1(2), 16(256) grey different colors */
  k = 2;
  
  neww = (double)width; newl = (double)length;
  /* scale size */
  rat = 500.0/neww; neww = 500.0; newl = rat*newl;
  
  if(newl>700.0)
    {
      rat = 700.0/newl; newl = 700.0; neww = rat*neww;
    }

  /* 8.5x11 page with no margins is 612x792 points */
  xorigin = (612 - neww)/2.0; yorigin = (792 - newl)/2.0;
  /* defines compliance */
  fprintf(FP,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  /* defines size of image */
  fprintf(FP,"%%%%BoundingBox: %f %f %f %f\n",
	  xorigin-20, yorigin-20,xorigin+neww+5,yorigin+newl+5);
  /* save present graphics state */
  fprintf(FP,"gsave\n");
 
  //axes_Definition(FP, xorigin, yorigin);
  draw_box_around_image(FP, xorigin, yorigin, neww, newl);
  
  /* define procedure "bufstr" to read strings of "width" characters */
  fprintf(FP,"/bufstr %d string def\n\n", k*width);
  /* move to origin of image */
  fprintf(FP,"%f %f translate\n",xorigin,yorigin);
  /* define the scale (default is 1 pt) */
  fprintf(FP,"%f %f scale\n\n",neww,newl);
  /* width, height, and number of bits per pixel */
  fprintf(FP,"%d %d %d\n",width,length,4*k);
  /* a matrix definition */
  fprintf(FP,"[%d 0 0 %d 0 0]\n",width,length);
  /* guidelines on how to read the image data */
  fprintf(FP,"{currentfile bufstr readhexstring pop} bind image\n");  

  for(i=0; i<width; i++)
    for(j=0; j<length; j++){
      z = (Data[j][i] - Min_z)/(Max_z - Min_z);

      if(z > 0)
	grey = (unsigned char) (255.0 * (z + zero_offset)/(zero_offset + 1.));
      else
	grey = (unsigned char) (255.0 * (0. + zero_offset)/(zero_offset + 1.));

      fprintf(FP,"%.2x",grey);
      if((j+1)%60 == 0) fprintf(FP,"\n");
    }

  fprintf(FP,"\n");
  fprintf(FP,"grestore\n");
  fprintf(FP,"showpage\n");

  fclose(FP);
}

void vertical_Scale_coherence(float **Data, int NPTS, float *Max, float *Min)
{
  int i,j;
  
  for(i=0; i<NPTS; i++)
    for(j=0; j<NPTS; j++)
      if(Data[i][j] > 0.){ 
	*Min = Data[0][0];
	break;
      }
  
  *Max = Data[0][0]; 
  for(i=0; i<NPTS; i++)
    for(j=0; j<NPTS; j++){
      if(Data[i][j] > 0.){
	*Max = MAX(*Max, Data[i][j]);
	*Min = MIN(*Min, Data[i][j]);
      }
    }

  printf("Max Value: %f\tMin Value: %f\n", *Max, *Min);
}














