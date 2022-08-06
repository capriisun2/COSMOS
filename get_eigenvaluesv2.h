#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
Written by ENB
V1.0: 7/19/22 - initial
V1.1 = 7/22/22 - add sorting, use optimium work size
To use, call: get_eigenvalues(int N, double mat[N][N], double vals_real_return[N],double vals_im_return[N],double vecs_return[N][N], sortflag)

N and mat are inputs, vals and vecs are outputs. If sortflag==1, the results are sorted by the magnitude of the real component of the eigenvalue.

This header requires OpenBLAS.

To compile file "test.c", run: gcc test.c -lopenblas
*/


//Import DGEEV
extern void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );
                
                
int cmp(const void *pa, const void *pb){
	const double *a = *(const double **)pa;
	const double *b = *(const double **)pb;
	if (a[0] == b[0]){
		return 0;
		}
	if (a[0] > b[0]){
		return 1;
		}
	if (a[0] < b[0]){
		return -1;
		}
	
}
void get_eigenvalues(int N, double mat[N][N], double vals_real_return[N],double vals_im_return[N],double vecs_return[N][N],int sortflag){
	//COPY FOR MEMORY SAFETY
	int LWORK=-1;
	double work_test[1];
	double mat_copy[N*N];
	double vecs_return_copy[N*N];
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			mat_copy[N*(i)+j]=mat[j][i];
		}
	}
	
	int info;
	dgeev_( "N","V", &N, mat_copy, &N, vals_real_return, vals_im_return, NULL, &N, vecs_return_copy, &N, work_test, &LWORK, &info );
        LWORK=work_test[0];
        double *work = malloc (sizeof (double) * LWORK);
        memset (work, 0, sizeof (double) * LWORK);
        
        dgeev_( "N","V", &N, mat_copy, &N, vals_real_return, vals_im_return, NULL, &N, vecs_return_copy, &N, work, &LWORK, &info );
        
        free(work);
        for (int i = 0; i<N; i++){
        	for(int j = 0; j<N; j++){
         		vecs_return[i][j] = vecs_return_copy[N*(i)+j];
         	}
         }
        if (sortflag == 1){
        	double **eigsortarr;
        	eigsortarr = malloc(N * sizeof(double*));
    		for (int i = 0; i < N; i++){
        		eigsortarr[i] = malloc(2 * sizeof(double));
        		eigsortarr[i][0]=vals_real_return[i];
        		eigsortarr[i][1] = i;
   		 }
   		 
        	qsort(eigsortarr,N,sizeof(eigsortarr[0]),cmp);
        	
        	double vals_im_copy[N];
        	
        	//double vecs_copy[N][N];
        	double **vecs_copy;
        	vecs_copy = malloc(N * sizeof(double*));
    		for (int i = 0; i < N; i++){
    			vals_im_copy[i] = vals_im_return[i];
        		vals_real_return[i]=eigsortarr[i][0];
        		vecs_copy[i] = malloc(N * sizeof(double));
        		for (int j=0; j<N;j++){
        			vecs_copy[i][j] = vecs_return[i][j];
        		}
        		
   		 }
        	
        	
        	
        	for (int i=0; i<N;i++){
        		vals_im_return[i] = vals_im_copy[(int) eigsortarr[i][1]];
        		for (int j=0; j<N;j++){
        			vecs_return[i][j] = vecs_copy[(int) eigsortarr[i][1]][j];
        		}
        	}
        	 for (int i = 0; i < N; i++){
        		free(eigsortarr[i]);
   		 }
   		free(eigsortarr);
   		free(vecs_copy);
        }
        if (info !=0){
         printf("%s","error");
         }
 
}

