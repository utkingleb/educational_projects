#define _CRT_SECURE_NO_WARNINGS

#define five 5
#define six 6
#define two 2.0
#define eigh 8
#define seven 7
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <string.h>

#define GIGA_MODIFIER 1e9
#define eps 1e-10

//Wrapper function for clock_gettime
unsigned long long currentTimeNano_Pr()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

int maxi(int x, int y)
{
    return (x > y) ? x : y;
}

//Euclidean distance
double norma(const double *a, int siz)
{
    double norm = 0;
    int i;
    for (i = 0; i < siz; i++) {
        norm += a[i] * a[i];
    }
    return sqrt(norm);
}

//length of input file
int length(FILE *f)
{
    int k = 0;
    double g;
    while (fscanf(f, "%lf ", &g) == 1) {
        k++;
    }
    rewind(f);
    //printf("length :%d\n", k);
    return k;
}

//matrix cell value
double mesh(int ord, int size, int i, int j)
{
    double res = 0;
    switch (ord) {
    case 1: {
        res = (double)(size - maxi(i, j) + 1);
        break;
    }
    case 2: {
        res = (double)(maxi(i, j));
        break;
    }
    case 3: {
        res = (double)(abs(i - j));
        break;
    }
    case 4: {
        res = 1. / (double)(i + j - 1);
        break;
    }
    }
    return res;
}

//main matrix of linear system formation
int matrform(int argc, double *M, int siz, int ord, char *filename, int number,
             int numer, int own_start, int own_end, int *kontakt_i)
{
    int i;
    int j;

    switch (argc) {
    //reading matrix from input file
    case five: {
        if (numer == number - 1) {
            for (i = 0; i < (siz / number) * siz; i++) {
                M[i] = 0.;
            }
            FILE *fp;
            kontakt_i[0] = 1;
            fp = fopen(filename, "r");
            if (!fp) {
                printf("The specified file does not exist\n");
                kontakt_i[0] = 0;
                for (i = 0; i < number - 1; i++) {
                    MPI_Send(M, (siz / number) * siz, MPI_DOUBLE, i, 1,
                             MPI_COMM_WORLD);
                    MPI_Send(kontakt_i, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                }
                return -3;
            }
            fseek(fp, 0, SEEK_END);
            long pos = ftell(fp);
            rewind(fp);
            long leng = length(fp);
            if ((pos <= 0) || (leng == 0)) {
                printf("File is empty\n");
                fclose(fp);
                kontakt_i[0] = 0;
                for (i = 0; i < number - 1; i++) {
                    MPI_Send(M, (siz / number) * siz, MPI_DOUBLE, i, 1,
                             MPI_COMM_WORLD);
                    MPI_Send(kontakt_i, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                }
                return -3;
            }
            if (leng != siz * siz) {
                printf("Wrong symbols in file\n");
                fclose(fp);
                kontakt_i[0] = 0;
                for (i = 0; i < number - 1; i++) {
                    MPI_Send(M, (siz / number) * siz, MPI_DOUBLE, i, 1,
                             MPI_COMM_WORLD);
                    MPI_Send(kontakt_i, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                }
                return -3;
            }

            int period;
            for (period = 0; period < number - 1; period++) {

                for (i = 0; i < siz / number; i++) {
                    for (j = 0; j < siz; j++) {
                        if (fscanf(fp, "%lf", &M[i * siz + j]) <= 0) {
                            printf("Failed to read matrix from file\n");
                            fclose(fp);
                            kontakt_i[0] = 0;
                            for (i = period; i < number - 1; i++) {
                                MPI_Send(M, (siz / number) * siz, MPI_DOUBLE, i,
                                         1, MPI_COMM_WORLD);
                            }
                            for (i = 0; i < number - 1; i++) {
                                MPI_Send(kontakt_i, 1, MPI_INT, i, 2,
                                         MPI_COMM_WORLD);
                            }
                            return -3;
                        }
                    }
                }
                MPI_Send(M, (siz / number) * siz, MPI_DOUBLE, period, 1,
                         MPI_COMM_WORLD);
            }
            for (i = 0; i < own_end - own_start; i++) {
                for (j = 0; j < siz; j++) {
                    if (fscanf(fp, "%lf", &M[i * siz + j]) <= 0) {
                        printf("Failed to read matrix from file\n");
                        fclose(fp);
                        // MPI_Abort(MPI_COMM_WORLD, -3);
                        kontakt_i[0] = 0;
                        for (i = 0; i < number - 1; i++) {
                            MPI_Send(kontakt_i, 1, MPI_INT, i, 2,
                                     MPI_COMM_WORLD);
                        }
                        return -3;
                    }
                }
            }
            for (i = 0; i < number - 1; i++) {
                MPI_Send(kontakt_i, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
            }
            fclose(fp);
        } else {
            MPI_Status status;
            MPI_Recv(M, (own_end - own_start) * siz, MPI_DOUBLE, number - 1, 1,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(kontakt_i, 1, MPI_INT, number - 1, 2, MPI_COMM_WORLD,
                     &status);
            if (kontakt_i[0] == 0) {
                return -3;
            }
        }
        break;
    }
    //matrix generation
    case 4: {
        for (i = 0; i < own_end - own_start; i++) {
            for (j = 0; j < siz; j++) {
                M[i * siz + j] = mesh(ord, siz, i + own_start + 1, j + 1);
            }
        }
        break;
    }
    }
    return 0;
}

//right column of linear system formation
void R_form(double *R, const double *M, int siz, int border)
{
    int i;
    int j;
    for (i = 0; i < border; i++) {
        R[i] = 0.;
        for (j = 0; j <= (siz - 1) / 2; j++) {
            R[i] += M[i * siz + (2 * j)];
        }
    }
}

//printing of linear system
void matr_print(const double *M, const double *R, double *Buffer,
                double *Buffer_R, int l, int m, int number, int numer,
                int own_start, int own_end, int ord)
{
    if (numer == number - 1) {
        int k = m;
        int i;
        int j;
        int start;
        int end;
        int period;
        MPI_Status status;
        if (k > l) {
            k = l;
        }

        for (period = 0; period < number - 1; period++) {
            start = (l / number) * period;
            end = (l / number) * (period + 1);
            if (end > m) {
                end = m;
            }
            if ((ord == -1) || (ord == 0)) {
                MPI_Recv(Buffer, (end - start) * l, MPI_DOUBLE, period, 1,
                         MPI_COMM_WORLD, &status);
            }
            if ((ord == 1) || (ord == 0)) {
                MPI_Recv(Buffer_R, (end - start), MPI_DOUBLE, period, 2,
                         MPI_COMM_WORLD, &status);
            }
            for (i = 0; i < end - start; i++) {
                if ((ord == -1) || (ord == 0)) {
                    for (j = 0; j < k; j++) {
                        printf("%10.3e", Buffer[i * l + j]);
                    }
                }
                if ((ord == 1) || (ord == 0)) {
                    printf("%10.3e", Buffer_R[i]);
                }
                printf("\n");
            }
            if ((l / number) * (period + 1) > m) {
                period = number;
            }
        }
        if (own_start < m) {
            if (own_end > m) {
                own_end = m;
            }
            for (i = 0; i < own_end - own_start; i++) {
                if ((ord == -1) || (ord == 0)) {
                    for (j = 0; j < k; j++) {
                        printf("%10.3e", M[i * l + j]);
                    }
                }
                if ((ord == 1) || (ord == 0)) {
                    printf("%10.3e", R[i]);
                }
                printf("\n");
            }
        }
        printf("\n");
    } else {
        if (own_start < m) {
            if (own_end > m) {
                own_end = m;
            }
            if ((ord == -1) || (ord == 0)) {
                MPI_Send(M, (own_end - own_start) * l, MPI_DOUBLE, number - 1,
                         1, MPI_COMM_WORLD);
            }
            if ((ord == 1) || (ord == 0)) {
                MPI_Send(R, (own_end - own_start), MPI_DOUBLE, number - 1, 2,
                         MPI_COMM_WORLD);
            }
        }
    }
}

//infinity norm for answer vector
double norm_infinit(const double *Part, const int border, double *kontakt_d,
                    const int number, const int numer, const int size)
{
    double max = 0.;
    MPI_Status status;
    int i;
    int j;
    for (i = 0; i < border; i++) {
        for (j = 0; j < size; j++) {
            if (fabs(Part[i * size + j]) > max) {
                max = fabs(Part[i * size + j]);
            }
        }
    }
    if (numer == number - 1) {
        for (i = 0; i < number - 1; i++) {
            MPI_Recv(kontakt_d, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);
            if (max < kontakt_d[0]) {
                max = kontakt_d[0];
            }
        }
        for (i = 0; i < number - 1; i++) {
            kontakt_d[0] = max;
            MPI_Send(kontakt_d, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    } else {
        kontakt_d[0] = max;
        MPI_Send(kontakt_d, 1, MPI_DOUBLE, number - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(kontakt_d, 1, MPI_DOUBLE, number - 1, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
    }
    return kontakt_d[0];
}

//solution of linear system by Gauss-Jordan elimination with selection of main element by row
int GJ_solution(double *Part, double *X, int size, int own_start, int own_end,
                const double epsi, int *Transp, int *kontakt_i, const int numer,
                const int number, double *Line)
{
    MPI_Status status;
    int iterator;
    int inner_iter;
    int t;

    double maxi;
    int maxi_p;
    int sender;
    double divid;
    double koef;

    for (iterator = 0; iterator < size; ++iterator) {
        for (sender = 0; sender < number - 1; sender++) {
            if (((size / number) * sender <= iterator) &&
                ((size / number) * (sender + 1) > iterator)) {
                break;
            }
        }

        if (sender == numer) {
            maxi = 0.;
            maxi_p = 0;
            for (inner_iter = iterator; inner_iter < size; inner_iter++) {
                if (fabs(Part[(iterator - own_start) * size +
                              Transp[inner_iter]]) > fabs(maxi)) {
                    maxi = Part[(iterator - own_start) * size +
                                Transp[inner_iter]];
                    maxi_p = inner_iter;
                }
            }

            if (fabs(maxi) < epsi) {
                kontakt_i[0] = 0;
                for (inner_iter = 0; inner_iter < number; inner_iter++) {
                    if (inner_iter == numer) {
                        continue;
                    }
                    MPI_Send(kontakt_i, 1, MPI_INT, inner_iter, 0,
                             MPI_COMM_WORLD);
                }
                return -4;
            } else {
                kontakt_i[0] = 1;
                for (inner_iter = 0; inner_iter < number; inner_iter++) {
                    if (inner_iter == numer) {
                        continue;
                    }
                    MPI_Send(kontakt_i, 1, MPI_INT, inner_iter, 0,
                             MPI_COMM_WORLD);
                }
            }

            if (maxi_p != iterator) {
                t = Transp[iterator];
                Transp[iterator] = Transp[maxi_p];
                Transp[maxi_p] = t;
            }

            kontakt_i[0] = maxi_p;
            for (inner_iter = 0; inner_iter < number; inner_iter++) {
                if (inner_iter == numer) {
                    continue;
                }
                MPI_Send(kontakt_i, 1, MPI_INT, inner_iter, 0, MPI_COMM_WORLD);
            }
            memcpy(Line, &Part[(iterator - own_start) * size],
                   size * sizeof(double));
            Line[size] = X[iterator - own_start];
            for (inner_iter = 0; inner_iter < number; inner_iter++) {
                if (inner_iter == numer) {
                    continue;
                }
                MPI_Send(Line, size + 1, MPI_DOUBLE, inner_iter, 0,
                         MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(kontakt_i, 1, MPI_INT, sender, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);
            if (kontakt_i[0] == 0) {
                return -4;
            }
            MPI_Recv(kontakt_i, 1, MPI_INT, sender, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);
            maxi_p = kontakt_i[0];
            if (maxi_p != iterator) {
                t = Transp[iterator];
                Transp[iterator] = Transp[maxi_p];
                Transp[maxi_p] = t;
            }
            MPI_Recv(Line, size + 1, MPI_DOUBLE, sender, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
        }

        divid = 1 / Line[Transp[iterator]];

        for (inner_iter = 0; inner_iter < own_end - own_start; inner_iter++) {
            if (inner_iter == iterator - own_start) {
                continue;
            }
            koef = -1 * Part[inner_iter * size + Transp[iterator]] * divid;
            for (t = iterator; t < size; t++) {
                Part[inner_iter * size + Transp[t]] += koef * Line[Transp[t]];
            }
            X[inner_iter] += koef * Line[size];
        }
    }

    for (inner_iter = 0; inner_iter < own_end - own_start; inner_iter++) {
        X[inner_iter] =
            X[inner_iter] /
            Part[inner_iter * size + Transp[own_start + inner_iter]];
    }
    return 0;
}

//calculation of residual of the solution
void residual(const double *Part, const int number, const int numer,
              double *Answer, const int own_start, const int own_end,
              const int size, const double *X, double *Send, const double *R)
{
    int i;
	int j;
	int period;
    double ret;
    MPI_Status status;
    if (numer == number - 1) {
        for (period = 0; period < number - 1; period++) {
            MPI_Recv(&Answer[(size / number) * period], size / number,
                     MPI_DOUBLE, period, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
        memcpy(&Answer[(size / number) * (number - 1)], X,
               (own_end - own_start) * sizeof(double));
        for (period = 0; period < number - 1; period++) {
            MPI_Send(Answer, size, MPI_DOUBLE, period, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Send(X, size / number, MPI_DOUBLE, number - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(Answer, size, MPI_DOUBLE, number - 1, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
    }

    for (i = own_start; i < own_end; i++) {
        Send[i] = 0.;
        for (j = 0; j < size; j++) {
            Send[i] += Part[(i - own_start) * size + j] * Answer[j];
            if (j == size - 1) {
                Send[i] -= R[i - own_start];
            }
        }
    }

    if (numer == number - 1) {
        for (period = 0; period < number - 1; period++) {
            MPI_Recv(&Send[(size / number) * period], size / number, MPI_DOUBLE,
                     period, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&Answer[(size / number) * period], size / number,
                     MPI_DOUBLE, period, 2, MPI_COMM_WORLD, &status);
        }
        ret = norma(Send, size) / norma(Answer, size);
        printf("Residual: %10.3e\n", ret);
    } else {
        MPI_Send(&Send[own_start], size / number, MPI_DOUBLE, number - 1, 1,
                 MPI_COMM_WORLD);
        MPI_Send(R, size / number, MPI_DOUBLE, number - 1, 2, MPI_COMM_WORLD);
    }
}

//calculation of error of the solution
void error(const double *X, int siz, const int number, const int numer,
           double *Line, const int own_start, const int own_end)
{
    int i;
    for (i = own_start; i < own_end; i++) {
        if (i % 2 == 0) {
            Line[i] = X[i - own_start] - 1.;
        } else {
            Line[i] = X[i - own_start];
        }
    }

    if (numer == number - 1) {
        MPI_Status status;
        for (i = 0; i < number - 1; i++) {
            MPI_Recv(&Line[(siz / number) * i], siz / number, MPI_DOUBLE, i,
                     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
        printf("Error: %10.3e\n", norma(Line, siz));
    } else {
        MPI_Send(&Line[own_start], siz / number, MPI_DOUBLE, number - 1, 0,
                 MPI_COMM_WORLD);
    }
}

