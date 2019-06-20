#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  
  if(argc!=2)
  {   
    printf("ERROR:  Usage: ./Practice3 <input file>");      
  }
  
  else{
  FILE * in;
  char tempbuff[400];

  in  = fopen(argv[1], "r");

  if(in == NULL){
    printf("file %s not found", argv[1]);
    }
    
  char tmpstr1[16];
  char tmpstr2[16];
  char InputDir[16];
  char InputFile[16];
  char OutputDir[16];
  char OutputFile[16];

  int InputFormat;
  int OutputFormat;
  int Nparticles;
  int nmax;
  int lmax;
  int cov_matrix;
  int BootstrapSample;
  int Nsampling;
  double r_s;
  

  fgets(tempbuff,400,in);
  sscanf(tempbuff, "%15s  %15s", tmpstr1, tmpstr2);

  
  if (strcmp(tmpstr1,"InputDir")==0){
     strcpy(InputDir, tmpstr2);

  printf("%s\n",  InputDir);  
  }

  fgets(tempbuff,400,in);
  sscanf(tempbuff, "%15s  %15s", tmpstr1, tmpstr2);

  if (strcmp(tmpstr1,"InputFile")==0) {
     strcpy(InputFile, tmpstr2);

  printf("%s\n",  InputFile);  
  }


  if (strcmp(tmpstr1,"InputFormat")==0) {
     InputFormat = atoi(tmpstr2);

  printf("%d\n",  InputFormat);  
  }

  
  if (strcmp(tmpstr1,"OutputDir")==0) {
     strcpy(OutputDir, tmpstr2);

  printf("%s\n",  OutputDir);  
  }


  if (strcmp(tmpstr1,"OutputFile")==0) {
     strcpy(OutputFile, tmpstr2);

  printf("%s\n",  OutputFile);  
  }

  
  if (strcmp(tmpstr1,"OutputFormat")==0) {
     OutputFormat =atoi(tmpstr2);

  printf("%d\n",  OutputFormat);  
  }


  if (strcmp(tmpstr1,"Nparticles")==0) {
     Nparticles =  atoi(tmpstr2);

  printf("%d\n",  Nparticles);  
  }


  if (strcmp(tmpstr1,"Rs")==0) {
     r_s =  atof(tmpstr2);

  printf("%f\n",  r_s);  
  }


  if (strcmp(tmpstr1,"nmax")==0) {
     nmax = atoi(tmpstr2);

  printf("%d\n",  nmax);  
  }

  if (strcmp(tmpstr1,"lmax")==0) {
     lmax = atoi(tmpstr2);

  printf("%d\n",  lmax);  
  }


  if (strcmp(tmpstr1,"CovMatrix")==0) {
     cov_matrix = atoi(tmpstr2);

  printf("%d\n",  cov_matrix);  
  }


  if (strcmp(tmpstr1,"BootstrapSample")==0) {
     BootstrapSample = atoi(tmpstr2);

  printf("%d\n",  BootstrapSample);  
  }


  if (strcmp(tmpstr1,"NSampling")==0) {
     Nsampling =atoi(tmpstr2);

  printf("%d\n",  Nsampling);  
  }
  

  fclose(in);
 
  }


          
  //memset(dest, '\0', sizeof(dest));
  //strcpy(src, "This is tutorialspoint.com");
  //strcpy(dest, src);

  //printf("Final copied string : %s\n", src);
  //printf("Final copied string : %s\n", dest);
                       
  return(0);
}
