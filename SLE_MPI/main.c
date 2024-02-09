#include "mpi.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main(int argc, char **argv)
{
    int numer;
    int number;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &numer);
    MPI_Comm_size(MPI_COMM_WORLD, &number);

    int siz;
    int m;
    int ord;
    int ord_cont;
    int i;
    double epsi;
    double *Part = NULL;
    double *copy_Part = NULL;
    double *R = NULL;
    double *X = NULL;
    double *Buffer = NULL;
    double *Buffer_R = NULL;
    char *endptr;
    char *filename = NULL;

    int *Transp = NULL;
    double *Line = NULL;

    int start;
    int end;
    double *kontakt_d = NULL;
    int *kontakt_i = NULL;

    unsigned long long times;

    MPI_Status status;

    //exception checking

    if ((argc < 4) || (argc > five)) {
        printf("Invalid number of arguments\n");
        MPI_Finalize();
        return -1;
    }
    siz = (int)strtol(argv[1], &endptr, deci);
    if (siz <= 0 || (errno == 0 && argv[1] && *endptr != 0)) {
        printf("Incorrect value for argument n\n");
        MPI_Finalize();
        return -1;
    }
    m = (int)strtol(argv[2], &endptr, deci);
    if (m <= 0 || m > siz || (errno == 0 && argv[1] && *endptr != 0)) {
        printf("Incorrect value for argument m\n");
        MPI_Finalize();
        return -1;
    }
    errno = 0;
    ord = (int)strtol(argv[3], &endptr, deci);
    if ((errno == 0 && argv[1] && *endptr != 0) || (argv[3] == endptr) ||
        ord < 0 || ord > 4 || (ord == 0 && argc == 4)) {
        printf("Incorrect value for argument k\n");
        MPI_Finalize();
        return -1;
    }

    if ((numer == number - 1) && (ord == 0)) {
        if (argv[4] == NULL) {
            MPI_Abort(MPI_COMM_WORLD, -1);
            return -1;
        }
        filename = argv[4];
    }

    start = (siz / number) * numer;
    if (numer == number - 1) {
        end = siz;
    } else {
        end = (siz / number) * (numer + 1);
    }

    //memory allocation

    if (!(kontakt_i = (int *)malloc(sizeof(int)))) {
        printf("Failed to allocate memory\n");
        MPI_Abort(MPI_COMM_WORLD, -2);
        return -2;
    }

    kontakt_i[0] = 1;

    if (!(Part = (double *)malloc((end - start) * siz * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        Part = NULL;
        kontakt_i[0] = 0;
        if (numer != number - 1) {
            MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
        }
    }
    if (kontakt_i[0] == 1) {
        if (!(copy_Part =
                  (double *)malloc((end - start) * siz * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            copy_Part = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (kontakt_i[0] == 1) {
        if (!(R = (double *)malloc((end - start) * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            R = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (kontakt_i[0] == 1) {
        if (!(X = (double *)malloc((end - start) * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            X = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (kontakt_i[0] == 1) {
        if (!(Transp = (int *)malloc(siz * sizeof(int)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            Transp = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (kontakt_i[0] == 1) {
        for (i = 0; i < siz; i++) {
            Transp[i] = i;
        }
        if (!(Line = (double *)malloc((siz + 1) * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            Line = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (kontakt_i[0] == 1) {
        if (!(kontakt_d = (double *)malloc(sizeof(double)))) {
            printf("Failed to allocate memory\n");
            kontakt_i[0] = 0;
            kontakt_d = NULL;
            if (numer != number - 1) {
                MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    if ((kontakt_i[0] == 1) && (numer == number - 1)) {
        if (!(Buffer =
                  (double *)malloc((siz / number) * siz * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            Buffer = NULL;
            kontakt_i[0] = 0;
        }
    }
    if ((kontakt_i[0] == 1) && (numer == number - 1)) {
        if (!(Buffer_R = (double *)malloc((siz / number) * sizeof(double)))) {
            printf("Failed to allocate memory\n");
            Buffer_R = NULL;
            kontakt_i[0] = 0;
        }
    }
    
    //MPI initialisation

    if (numer != number - 1) {
        MPI_Send(kontakt_i, 1, MPI_INT, number - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(kontakt_i, 1, MPI_INT, number - 1, MPI_ANY_TAG, MPI_COMM_WORLD,
                 &status);
    } else {
        ord_cont = kontakt_i[0];
        for (i = 0; i < number - 1; i++) {
            MPI_Status status;
            MPI_Recv(kontakt_i, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);
            if ((ord_cont == 1) && (kontakt_i[0] == 0)) {
                ord_cont = 0;
            }
        }

        kontakt_i[0] = ord_cont;
        for (i = 0; i < number - 1; i++) {
            MPI_Send(kontakt_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }

    if (kontakt_i[0] == 0) {
        if (numer == number - 1) {
            if (Buffer != NULL) {
                free(Buffer);
            }
            if (Buffer_R != NULL) {
                free(Buffer_R);
            }
        }
        free(kontakt_i);
        if (kontakt_d != NULL) {
            free(kontakt_d);
        }
        if (Part != NULL) {
            free(Part);
        }
        if (copy_Part != NULL) {
            free(copy_Part);
        }
        if (R != NULL) {
            free(R);
        }
        if (X != NULL) {
            free(X);
        }
        if (Transp != NULL) {
            free(Transp);
        }
        if (Line != NULL) {
            free(Line);
        }
        MPI_Finalize();
        return -2;
    }

    //creating matrices of the linear system
    int ret_code;
    ret_code = matrform(argc, Part, siz, ord, filename, number, numer, start,
                        end, kontakt_i);

    if (ret_code == -3) {
        if (numer == number - 1) {

            free(Buffer);

            free(Buffer_R);
        }
        free(kontakt_i);

        free(kontakt_d);

        free(Part);

        free(copy_Part);

        free(R);

        free(X);

        free(Transp);

        free(Line);

        MPI_Finalize();
        return -3;
    }
    R_form(R, Part, siz, end - start);
    memcpy(copy_Part, Part, (end - start) * siz * sizeof(double));
    memcpy(X, R, (end - start) * sizeof(double));

    epsi =
        eps * norm_infinit(Part, end - start, kontakt_d, number, numer, siz);

    
    matr_print(Part, R, Buffer, Buffer_R, siz, m, number, numer, start, end,
               0);

    times = currentTimeNano_Pr();
    
    //solution of the linear system; MAIN FUNCTION CALL
    ret_code = GJ_solution(copy_Part, X, siz, start, end, epsi, Transp,
                           kontakt_i, numer, number, Line);
    times = currentTimeNano_Pr() - times;

    if (ret_code == -4) {
        if (numer == number - 1) {

            free(Buffer);

            free(Buffer_R);
        }
        free(kontakt_i);

        free(kontakt_d);

        free(Part);

        free(copy_Part);

        free(R);

        free(X);

        free(Transp);

        free(Line);

        MPI_Finalize();
        return -4;
    }
 
    //results printing
    if (numer == number - 1) {
        printf("Time: %.2Lf\n", (long double)(times) / GIGA_MODIFIER);
        printf("Solution:\n");
    }
    matr_print(copy_Part, X, Buffer, Buffer_R, siz, m, number, numer, start,
               end, 1);
    printf("Time of process %d: %.2Lf\n", numer,
           (long double)(times) / GIGA_MODIFIER);
    residual(Part, number, numer, Line, start, end, siz, X, copy_Part, R);
    error(X, siz, number, numer, Line, start, end);

    free(Part);
    if (numer == number - 1) {
        free(Buffer);
    }
    MPI_Finalize();
    return 0;
}

