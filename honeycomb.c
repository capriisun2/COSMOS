//
// Created by eli on 7/25/22.
//
#include <stdlib.h>
#include <stdio.h>
#include "matrix_functions.h"
int main(){
FILE * fileout1;
FILE * fileout2;
fileout1=fopen("honeycomb_eig","w");
fileout2=fopen("honeycomb_mat","w");
    int L = 80;
    int H = 80;
    int N = L*H;
    Mat testmat,vecs;
    Vec vals_real, vals_im;

    makeMatrix(N,N,&testmat);
    makeMatrix(N,N,&vecs);
    makeVector(N,&vals_real);
    makeVector(N,&vals_im);
//    randomMatrix(&testmat);

double thop = 0.7,E = 6;
int i,j,k;
/* printf("\n Enter thop,E:  ");
scanf("%lf %lf",&thop,&E);*/

for (i=0; i<N; i=i+1){
for (j=0; j<N; j=j+1){
testmat.Matrix[i][j]    = 0.0; }}

for (i=0; i<N; i=i+1){
testmat.Matrix[i][i]    = E; }

for(i=0; i < L; i++) { //Rows
    k = i%2;

for(j = i*L; j < (i+1) * L; j++) { //Columns
    if(j % 2 == k && j < N-1) {
    testmat.Matrix[j][j+1] = thop;
    testmat.Matrix[j+1][j] = thop;
    }
    if(j < N-L) {
        testmat.Matrix[j][j+L] = thop;
        testmat.Matrix[j+L][j] = thop;
    }
}
}

   printf("Input matrix :\n");
   for (i=0; i<N; i=i+1){
   for (j=0; j<N; j=j+1){
   fprintf(fileout2, "%6.2lf",testmat.Matrix[i][j]);} 
   fprintf(fileout2,"\n"); 
   }

    get_eigenvalues(&testmat,&vals_real,&vals_im,&vecs,1);

//    printMatrix(&testmat); 

    printf("Real eigenparts :\n");
    for (i=0;i<N;i++){
    fprintf(fileout1,"%12.6lf \n",vals_real.Vector[i]);
    }
//    printVector(&vals_real);
//    printf("Imaginary eigenparts :\n");
//    printVector(&vals_im);
//    printf("Eigenvectors :\n");
//    printMatrix(&vecs);
  
    fclose(fileout1);
    fclose(fileout2);
    return 0;
}
