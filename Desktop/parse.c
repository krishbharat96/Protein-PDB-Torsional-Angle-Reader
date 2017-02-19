#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N_FRAMES  framenum
#define N_ATOMS   Atomsnum // actual number + 2
#define ATOM_C0   C0num
#define ATOM_N1   N1num
#define ATOM_Ca1  CAnum
#define ATOM_C1   C1num
#define ATOM_N2   N2num 	 

#define PI    3.14159

float angle(float x[], float y[], float z[], int n1, int n2, int n3, int n4);
float angle2(float x[], float y[], float z[], int n1, int n2, int n3, int n4);
void crossprod(float res[], float vec1[], float vec2[]);

int main(int argc, char* argv[]) {
  char filename[100];
  FILE *infile;

  strcpy(filename, argv[1]);
  infile = fopen(filename, "rb");
  if (!infile) {
    fprintf(stderr, "ERROR: unable to open %s.\n", filename);
    exit(1);
  }
/* read atoms.txt here and populate variables atom_c0, n1, ca1, atom_c1, atom_ca1, atom_n2 */

  float x[N_ATOMS], y[N_ATOMS], z[N_ATOMS];
  float phi, psi;

  int nC0, nN1, nCa1, nC1, nN2;
  nC0 = (ATOM_C0)-1;
  nN1 = (ATOM_N1)-1;
  nCa1 = (ATOM_Ca1)-1;
  nC1 = (ATOM_C1)-1;
  nN2 = (ATOM_N2)-1;

  fseek(infile, 280, SEEK_SET);

  int i;
  for (i=1; i<=N_FRAMES; i++) {
    fseek(infile, 56, SEEK_CUR);
    fread(x, 1, (N_ATOMS*4), infile);
    fread(y, 1, (N_ATOMS*4), infile);
    fread(z, 1, (N_ATOMS*4), infile);

    phi = angle2(x, y, z, nC0, nN1, nCa1, nC1);
    psi = angle2(x, y, z, nN1, nCa1, nC1, nN2);

    if (isnan(phi) == 1) {
      fprintf(stderr, "ERROR: NaN phi at frame %d.\n", i);
    } else if (isnan(psi) == 1) {
      fprintf(stderr, "ERROR: NaN psi at frame %d.\n", i);
    } else {
      printf("%d\t%.2f\t%.2f\n", i, phi, psi);
    }
    //printf("%f %f %f\n", x[nC0]-30.0, y[nC0]-30.0, z[nC0]-30.0);
  }

  fclose(infile);
  return 0;
}

// diheral angle -- this should not be used!
// ref: http://www.bio.net/mm/xtal-log/1997-February/002899.html
// ref: Glusker et al., "Crystal Structure Analysis for Chemists and Biologists", pp. 465-469
float angle(float x[], float y[], float z[], int n1, int n2, int n3, int n4) {
  float d12, d13, d14, d23, d24, d34;

  d12 = sqrt((x[n1] - x[n2])*(x[n1] - x[n2]) + (y[n1] - y[n2])*(y[n1] - y[n2]) + (z[n1] - z[n2])*(z[n1] - z[n2]));
  d13 = sqrt((x[n1] - x[n3])*(x[n1] - x[n3]) + (y[n1] - y[n3])*(y[n1] - y[n3]) + (z[n1] - z[n3])*(z[n1] - z[n3]));
  d14 = sqrt((x[n1] - x[n4])*(x[n1] - x[n4]) + (y[n1] - y[n4])*(y[n1] - y[n4]) + (z[n1] - z[n4])*(z[n1] - z[n4]));
  d23 = sqrt((x[n2] - x[n3])*(x[n2] - x[n3]) + (y[n2] - y[n3])*(y[n2] - y[n3]) + (z[n2] - z[n3])*(z[n2] - z[n3]));
  d24 = sqrt((x[n2] - x[n4])*(x[n2] - x[n4]) + (y[n2] - y[n4])*(y[n2] - y[n4]) + (z[n2] - z[n4])*(z[n2] - z[n4]));
  d34 = sqrt((x[n3] - x[n4])*(x[n3] - x[n4]) + (y[n3] - y[n4])*(y[n3] - y[n4]) + (z[n3] - z[n4])*(z[n3] - z[n4]));

  float P, Q;

  P = pow(d12,2) * ( pow(d23,2) + pow(d34,2) - pow(d24,2)) + \
      pow(d23,2) * (-pow(d23,2) + pow(d34,2) + pow(d24,2)) + \
      pow(d13,2) * ( pow(d23,2) - pow(d34,2) + pow(d24,2)) - \
      2 * pow(d23,2) * pow(d14,2);

  Q = (d12 + d23 + d13) * ( d12 + d23 - d13) * \
      (d12 - d23 + d13) * (-d12 + d23 + d13) * \
      (d23 + d34 + d24) * ( d23 + d34 - d24) * \
      (d23 - d34 + d24) * (-d23 + d34 + d24);
/*
  float vec1[3], vec2[3];

  vec1[0] = x[n2] - x[n1];
  vec1[1] = y[n2] - y[n1];
  vec1[2] = z[n2] - z[n1];

  vec2[0] = x[n3] - x[n4];
  vec2[1] = y[n3] - y[n4];
  vec2[2] = z[n3] - z[n4];

  float inner, abs1, abs2;

  inner = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
  abs1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
  abs2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
*/
  return 180.0*acos(P/sqrt(Q))/PI;
}

// ref: http://www.math.fsu.edu/~quine/MB_10/6_torsion.pdf
float angle2(float x[], float y[], float z[], int n1, int n2, int n3, int n4) {
  float b1[3], b2[3], b3[3];

  b1[0] = x[n1] - x[n2];
  b1[1] = y[n1] - y[n2];
  b1[2] = z[n1] - z[n2];

  b2[0] = x[n2] - x[n3];
  b2[1] = y[n2] - y[n3];
  b2[2] = z[n2] - z[n3];

  b3[0] = x[n3] - x[n4];
  b3[1] = y[n3] - y[n4];
  b3[2] = z[n3] - z[n4];

  float b12[3], b23[3];
  crossprod(b12, b1, b2);
  crossprod(b23, b2, b3);

  float b212[3];
  crossprod(b212, b2, b12);

  float arg1, arg2;
  //arg1 = b23[0]*b12[0] + b23[1]*b12[1] + b23[2]*b12[2];
  //arg2 = (b23[0]*b212[0] + b23[1]*b212[1] + b23[2]*b212[2]) / sqrt(b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2]);
  arg1 = -(b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2])*(b1[0]*b3[0] + b1[1]*b3[1] + b1[2]*b3[2]) + (b1[0]*b2[0] + b1[1]*b2[1] + b1[2]*b2[2])*(b2[0]*b3[0] + b2[1]*b3[1] + b2[2]*b3[2]);
  arg2 = sqrt(b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2])*(b1[0]*b23[0] + b1[1]*b23[1] + b1[2]*b23[2]);

  return -180.0*atan2(arg2,arg1)/PI;
}

void crossprod(float res[], float vec1[], float vec2[]) {
  res[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  res[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  res[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}
