#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix_functions.h"
#include "get_eigenvaluesv2.h"

int main(){
int i,j,N =33;
double thop=0.9,E=2.4;
FILE * fileout1;
FILE * fileout2;
FILE * fileout3;
fileout1 = fopen("eigvals", "w");
fileout2 = fopen("matrix", "w");
fileout3 = fopen("eigvecs", "w");

double mat[N][N];
double newMat[N][N];

for (i =0; i<N; i++){
   for (j =0; j<N; j++){
      mat[i][j]=0.0;
   }
}


for (i=0; i<N; i++){
    mat[i][i] = E;
}


mat[0][1]    = thop;
for (i=1; i<N-1; i++){
   mat[i][i-1] = thop;
   mat[i][i+1] = thop;
}

mat[N-1][N-2] = thop;

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


// for(i = 0; i < N; i++)
// {
//     for(j = 0; j < N; j++)
//     {
//         newMat[i][j] = vecs[i] * vecs [j];
//     }
// }


printf("\nEigenvalues: \n");
for (i=0; i<N;i++){
fprintf(fileout1,"%12.6lf",val_r[i]);
}

printf("\nEigenVectors: \n");
for (j=0; j<N;j++){
   for (i=0; i<N;i++){
      fprintf(fileout3, "%8.4lf ",vecs[i][j]);
   }
   printf("\n");
}
fclose(fileout1);
fclose(fileout2);
fclose(fileout3);
return 0;
}