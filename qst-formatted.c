#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "get_eigenvaluesv2.h"

int main() {
  int k, i, j, N = 7;
  double thop = 1, E = 0, t;
  FILE* eigvalslog;
  FILE* matrixlog;
  FILE* eigvectslog;
  FILE* psilog;
  FILE* problog;
  eigvalslog = fopen("out/eigvals", "w");    // fileout1
  matrixlog = fopen("out/matrix", "w");      // fileout2
  eigvectslog = fopen("out/eigvects", "w");  // fileout3
  psilog = fopen("out/psi", "w");            // fileout4
  problog = fopen("out/prob", "w");          // fileout5

  double mat[N][N];
  double newMat[N][N];
  double mat_r[N][N];
  double mat_i[N][N];
  double psi[N];
  double psi_r[N];
  double psi_i[N];
  double prob[N];
  double mat[N][N];
  double newMat[N][N];
  double mat_r[N][N];
  double mat_i[N][N];
  double psi[N];
  double psi_r[N];
  double psi_i[N];
  double newpsi_r[N];
  double newpsi_i[N];
  double prob[N];

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      mat[i][j] = 0.0;
      newMat[N - 1][N - 1] = 0.0;
      mat_r[N - 1][N - 1] = 0.0;
      mat_i[N - 1][N - 1] = 0.0;
    }
  }

  for (i = 0; i < N; i++) {
    psi_r[i] = 0.0;
    psi_i[i] = 0.0;
    psi[N - 1] = 0.0;
    prob[N - 1] = 0.0;
  }
  psi_r[0] = 1.0;
  psi[0] = 1.0;
  prob[0] = 1.0;

  for (i = 0; i < N; i++) {
    mat[i][i] = E;
  }

  for (i = 1; i < N; i++) {
    thop = sqrt(i * (N - i));
    mat[i][i - 1] = thop;
    mat[i - 1][i] = thop;
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      newMat[i][j] = mat[i][j];
    }
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      fprintf(matrixlog, "%8.4lf ", mat[i][j]);
    }
    fprintf(matrixlog, "\n ");
  }

  double val_r[N], val_i[N], vecs[N][N];
  get_eigenvalues(N, mat, val_r, val_i, vecs, 1);

  fprintf(psilog, "%lf", psi[0]);
  for (t = 0; t < 4; t = t + 0.01) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        mat_r[N - 1][N - 1] = 0.0;
        mat_i[N - 1][N - 1] = 0.0;
      }
    }
    for (k = 0; k < N; k++) {
      for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
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
          psi_r[i] += ((mat_r[i][j] * psi_r[i]) - (mat_i[i][j] * psi_i[i]));
          psi_i[i] += ((mat_r[i][j] * psi_i[i]) + (mat_i[i][j] * psi_r[i]));
          psi[i] = psi_r[i] + psi_i[i];
          fprintf(psilog, "%12.8lf", psi[i]);
          prob[i] = (psi_r[i] * psi_r[i]) + (psi_i[i] * psi_i[i]);
          fprintf(problog, "%12.8lf", prob[i]);
        }
      }
    }
  }

  for (i = 0; i < N; i++) {
    fprintf(eigvalslog, "%12.6lf", val_r[i]);
  }

  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      fprintf(eigvectslog, "%8.4lf", vecs[i][j]);
    }

    fprintf(eigvectslog, "\n");
  }
  fclose(eigvalslog);
  fclose(matrixlog);
  fclose(eigvectslog);
  fclose(psilog);
  fclose(problog);
  return 0;
  for (t = 0; t < 4; t = t + 0.01) {
    fprintf(problog, "timestep %lf \n", t);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        mat_r[i][j] = 0.0;
        mat_i[i][j] = 0.0;
      }
    }
    for (k = 0; k < N; k++) {
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
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
        }
      }
    }
    double norm = 0;
    for (int i = 0; i < N; i++) {
      newpsi_r[i] = 0.0;
      newpsi_i[i] = 0.0;
      for (int j = 0; j < N; j++) {
        newpsi_r[i] += ((mat_r[i][j] * psi_r[j]) - (mat_i[i][j] * psi_i[j]));
        newpsi_i[i] += ((mat_r[i][j] * psi_i[j]) + (mat_i[i][j] * psi_r[j]));
      }
      psi_r[i] = newpsi_r[i];
      psi_i[i] = newpsi_i[i];
      prob[i] = (psi_r[i] * psi_r[i]) + (psi_i[i] * psi_i[i]);
      fprintf(problog, "%6i %12.8lf \n", i, prob[i]);
    }
  }

  printf("\nEigenvalues: \n");
  for (i = 0; i < N; i++) {
    fprintf(eigvalslog, "%12.6lf", val_r[i]);
  }

  printf("\nEigenVectors: \n");
  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      fprintf(eigvectslog, "%8.4lf", vecs[i][j]);
    }

    fprintf(eigvectslog, "\n");
  }
  fclose(eigvalslog);
  fclose(matrixlog);
  fclose(eigvectslog);
  fclose(psilog);
  fclose(problog);
  return 0;
}
