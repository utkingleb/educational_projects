#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>
#include <time.h>

int matrform(int argc, double *M, int siz, int ord, char *filename);
double mesh(int ord, int size, int i, int j);
int matr_print(double *M, int l, int m, double *R, unsigned int fl_R);
unsigned long long currentTimeNano_Pr();
unsigned long long currentTimeNano_Thr();

void R_form(double *R, const double *M, int siz);

double norm_infinit(const double *M, int siz_real, int siz_cur);

int GJ_solution(int size, double *M, double *R, const double *proto_M,
                const double *proto_R, int threads, unsigned long long *times,
                double *Buffer_Transp, double *res1, double *res2,
                double *Copy, int *Transp);

double norma(const double *a, int siz);

void error(const double *X, int siz);

int length(FILE *f);

struct shipment;
void *minus_spread(void *b);

#endif
