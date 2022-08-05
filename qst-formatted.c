#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_eigenvaluesv2.h"

int main(){
int k,i,j,N =7;
double thop=1,E=0,t;
FILE * fileout1;
FILE * fileout2;
FILE * fileout3;
FILE * fileout4;
FILE * fileout5;
fileout1 = fopen("eigvals", "w");
fileout2 = fopen("matrix", "w");
fileout3 = fopen("eigvects", "w");
fileout4 = fopen("psi", "w");
fileout5 = fopen("prob", "w");

double mat[N][N];
double newMat[N][N];
double mat_r[N][N];
double mat_i[N][N];
double psi[N];
double psi_r[N];
double psi_i[N];
double prob[N];

for (i = 0; i < N; i++){
   for (j =0; j<N; j++){
      mat[i][j]=0.0;
      newMat[N][N] = 0.0;
      mat_r[N][N] = 0.0;
      mat_i[N][N] =0.0;
   }
}

for (i =0; i<N; i++){
      psi_r[i] = 0.0;
      psi_i[i] = 0.0;
      psi[N] = 0.0;
      prob[N] = 0.0;
}
psi_r[0] = 1.0;
psi[0] = 1.0;
prob[0] = 1.0;

// mat[0][0] = 11;
// mat[0][1] = 8;
// mat[1][0] = 8;
// mat[1][1] = -1;

for (i = 0; i < N; i++){
    mat[i][i] = E;
}


for (i = 1; i < N; i++){
   thop = sqrt(i * (N-i));
   mat[i][i-1] = thop;
   mat[i-1][i] = thop;
}


for(i = 0; i < N; i++)
{
   for(j = 0; j < N; j++)
   {
      newMat[i][j] = mat[i][j];
   }
}

printf("\n Matrix to be diagonalized:\n");
for (i=0; i<N; i++){
   printf("\n ");
   for (j=0; j<N; j++){
      fprintf(fileout2,"%8.4lf ", mat[i][j]);
   }
   fprintf(fileout2,"\n ");
}
printf("\n ");

double val_r[N],val_i[N],vecs[N][N];
get_eigenvalues(N, mat, val_r, val_i, vecs,1);

fprintf(fileout4,"%/lf", psi[0]);
for(t = 0; t < 4; t = t + 0.01) 
{
   for (int i = 0; i < N; i++){
   for (int j =0; j<N; j++){
      mat_r[N][N] = 0.0;
      mat_i[N][N] =0.0;
   } }
   for(k = 0; k < N; k++) 
   {
      for(int i = 0; i < N; i++)
      {
         for(int j = 0; j < N; j++)
          {
          newMat[i][j] += ((vecs[k][i] * vecs[k][j]) * val_r[k]);
            // double norm;
            // double denomSum;
            // denomSum=0;
            // for(int h = 0; h < N; h++) 
            //    {
            //    denomSum += pow(vecs[k][h],2); 
            //    }
            // norm = 1/denomSum;
            mat_r[i][j] += ((vecs[k][i] * vecs[k][j]) * cos(val_r[k] * t));
            mat_i[i][j] -= ((vecs[k][i] * vecs[k][j]) * sin(val_r[k] * t));
            psi_r[i] += ((mat_r[i][j]*psi_r[i]) - (mat_i[i][j]*psi_i[i]));
            psi_i[i] += ((mat_r[i][j]*psi_i[i]) + (mat_i[i][j]*psi_r[i]));
            psi[i] = psi_r[i] + psi_i[i];
            fprintf(fileout4,"%12.8lf", psi[i]);
            prob[i] = (psi_r[i]  * psi_r[i]) + (psi_i[i] * psi_i[i]);
            fprintf(fileout5,"%12.8lf", prob[i]);
          }
      }
   }
}







printf("\nEigenvalues: \n");
for (i = 0; i < N;i++){
fprintf(fileout1,"%12.6lf",val_r[i]);
}

printf("\nEigenVectors: \n");
for (j = 0; j < N;j++){
   for (i = 0; i < N; i++){
      
      fprintf(fileout3, "%8.4lf",vecs[i][j]);
     
   }
      
   fprintf(fileout3,"\n");
   }
fclose(fileout1);
fclose(fileout2);
fclose(fileout3);
fclose(fileout4);
fclose(fileout5);
return 0;
}
