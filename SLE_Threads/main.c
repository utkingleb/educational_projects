#include <errno.h>
#include <fenv.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "func.h"

#define deci 10
#define five 5
#define six 6
#define eigh 8
#define seven 7
#define thous 1000
#define GIGA_MODIFIER 1e9
#define NANO_MODIFIER 1e-9
#define eps 1e-10

int main(int argc, char *argv[])
{
    int siz;
    int m;
    int i;
    int ord;
    int thr;
    double *M;
    double *copy_M;
    double *R;
    double *X;
    int *Transp;
    double *Buffer_Transp;
    double *res1;
    double *res2;
    double *Copy;
    char *endptr;
    char *filename;
    unsigned long long *times;
    
    //exception checking

    if ((argc < five) || (argc > six)) {
        printf("Invalid number of arguments\n");
        return -1;
    }
    thr = (int)strtol(argv[1], &endptr, deci);
    if (thr <= 0) {
        printf("Incorrect value for argument threads\n");
        return -1;
    }
    siz = (int)strtol(argv[2], &endptr, deci);
    if (siz <= 0 || (errno == 0 && argv[2] && *endptr != 0)) {
        printf("Incorrect value for argument n\n");
        return -1;
    }
    m = (int)strtol(argv[3], &endptr, deci);
    if (m <= 0 || m > siz || (errno == 0 && argv[2] && *endptr != 0)) {
        printf("Incorrect value for argument m\n");
        return -1;
    }
    errno = 0;
    ord = (int)strtol(argv[4], &endptr, deci);
    if ((errno == 0 && argv[2] && *endptr != 0) || (argv[4] == endptr) ||
        ord < 0 || ord > 4 || (ord == 0 && argc == 4)) {
        printf("Incorrect value for argument k\n");
        return -1;
    }

    if (ord == 0) {
        if (argv[five] == NULL) {
            return -1;
        }
        filename = argv[five];
    } else {
        filename = NULL;
    }

    //memory allocation

    if (!(M = (double *)malloc(siz * siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(copy_M = (double *)malloc(siz * siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(M);
        return -2;
    }

    if (!(R = (double *)malloc(siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(M);
        free(copy_M);
        return -2;
    }

    //creating matrices of the linear system
    int er_code = matrform(argc, M, siz, ord, filename);
    if (er_code == 0) {
        R_form(R, M, siz);
        matr_print(M, siz, m, R, 1);
    } else {
        free(M);
        free(R);
        free(copy_M);
        return er_code;
    }

    //continuation of memory allocation
    
    if (!(X = (double *)malloc(siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(M);
        free(R);
        free(copy_M);
        return -2;
    }

    if (!(times = (unsigned long long *)malloc((thr * eigh + 1) *
                                               sizeof(unsigned long long)))) {
        printf("Failed to allocate memory\n");
        free(M);
        free(copy_M);
        free(R);
        free(X);
        return -2;
    }

    if (!(Transp = (int *)malloc(siz * thr * sizeof(int)))) {
        printf("Failed to allocate memory\n");
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(times);
        return -2;
    }

    if (!(Buffer_Transp = (double *)malloc(siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(Transp);
        free(times);
        return -2;
    }
    if (!(res1 = (double *)malloc(sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(Transp);
        free(times);
        free(Buffer_Transp);
        return -2;
    }
    if (!(res2 = (double *)malloc(sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(Transp);
        free(times);
        free(Buffer_Transp);
        free(res1);
        return -2;
    }
    
    if (!(Copy = (double *)malloc(thr * (siz + 1) * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(times);
        free(Transp);
        free(Buffer_Transp);
        free(res1);
        free(res2);
        return -2;
    }
    
    times[0] = 0;
    times[0] = currentTimeNano_Pr();
    
    //solution of the linear system; MAIN FUNCTION CALL
    
    er_code = GJ_solution(siz, copy_M, X, M, R, thr, times,
                          Buffer_Transp, res1, res2, Copy, Transp);
    if (er_code == -4) {
        free(X);
        free(M);
        free(R);
        free(copy_M);
        free(times);
        free(Transp);
        free(Buffer_Transp);
        free(res1);
        free(res2);
        free(Copy);
        return -4;
    }
    
    //results printing
    times[0] = currentTimeNano_Pr() - times[0];
    printf("Time: %.2Lf\n", (long double)(times[0]) / GIGA_MODIFIER);
    for (i = 0; i < thr; i++) {
        printf("Time of thread %d: %.2Lf\n", (i + 1),
               (long double)(times[eigh * i + 1]) / GIGA_MODIFIER);
    }
    printf("Solution:\n");
    matr_print(X, 1, m, R, 0);
    printf("Residual: %10.3e\n", res1[0]);
    error(X, siz);
    
    //freeing up memory
    free(M);
    free(R);
    free(X);
    free(copy_M);
    free(times);
    free(Transp);
    free(Buffer_Transp);
    free(res1);
    free(res2);
    free(Copy);
    return 0;
}

