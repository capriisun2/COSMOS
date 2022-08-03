#include <stdlib.h>
#include <stdio.h>
#include "get_eigenvaluesv2.h"            
int main(){
int i,j,N;
double thop1,thop2,E1;
FILE * fileout1;
FILE * fileout2;
fileout1 = fopen("eigHist15", "w");
fileout2 = fopen("matrix", "w");

printf("\n Enter N,thop1,thop2,E1:  ");
scanf("%i %lf %lf %lf",&N,&thop1,&thop2,&E1);

double mat[N][N];

for (i =0; i<N; i++){
   for (j =0; j<N; j++){
      mat[i][j]=0.0;
   }
}


for (i=0; i<N; i++){
if(i % 2 == 0) 
	{ mat[i][i] = E1;}
else
	{ mat[i][i] = E1;}
}


mat[0][1]    = thop1;
for (i=1; i<N-1; i++){
   mat[i][i-1] = thop1;
   mat[i][i+1] = thop1;
}
for (i=0; i<5; i++){
   mat[i+5][i] = thop2;
}
for(i=5; i<N; i++){
   mat[i-5][i] = thop2;
}
mat[(N/2)-1][(N/2)] = 0.0;
mat[(N/2)][(N/2)-1] = 0.0;
mat[N-1][N-2] = thop1;

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

printf("\nEigenvalues: \n");
for (i=0; i<N;i++){
fprintf(fileout1,"%12.6lf %12.6lf\n",val_r[i],val_i[i]);
}

printf("\nEigenVectors: \n");
for (j=0; j<N;j++){
   for (i=0; i<N;i++){
      printf("%8.4lf ",vecs[i][j]);
   }
   printf("\n");
}
fclose(fileout1);
fclose(fileout2);
return 0;
}
