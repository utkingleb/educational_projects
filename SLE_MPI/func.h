#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>

unsigned long long currentTimeNano_Pr();

int maxi(int x, int y);
double norma(const double *a, int siz);
int length(FILE *f);
int matrform(int argc, double *M, int siz, int ord, char *filename, int number,
             int numer, int own_start, int own_end, int *kontakt_i);
void R_form(double *R, const double *M, int siz, int border);
double mesh(int ord, int size, int i, int j);
void matr_print(const double *M, const double *R, double *Buffer,
                double *Buffer_R, int l, int m, int number, int numer,
                int own_start, int own_end, int ord);

void R_form(double *R, const double *M, int siz, int border);

double norm_infinit(const double *Part, const int border, double *kontakt_d,
                    const int number, const int numer, const int size);
int GJ_solution(double *Part, double *X, int size, int own_start, int own_end,
                const double epsi, int *Transp, int *kontakt_i, const int numer,
                const int number, double *Line);
void residual(const double *Part, const int number, const int numer,
              double *Answers, const int own_start, const int own_end,
              const int size, const double *X, double *Send, const double *R);
void error(const double *X, int siz, const int number, const int numer,
           double *Line, const int own_start, const int own_end);

int length(FILE *f);

#endif
