#define five 5
#define six 6
#define two 2.0
#define eigh 8
#define seven 7
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#define GIGA_MODIFIER 1e9
#define eps 1e-10

//struct for shipment between threads and main process
struct shipment {
    double *M;
    double *R;
    const double *proto_M;
    const double *proto_R;
    int *Transp;
    double *Buffer_Transp;
    unsigned long long *times;
    double *residual1;
    double *residual2;
    double *Copy;
    //double *Max_linie;
	//int *Max_p_linie;
	
    int size;
    int threads;
    int numer;
    double epsi;
    
    int ende;
    
    double *global_maxim;
    pthread_barrier_t *bar;
};

int max(int x, int y)
{
    return (x > y) ? x : y;
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
    printf("length :%d\n", k);
    return k;
}

//euclidian distance
double norma(const double *a, int siz)
{
    double norm = 0;
    int i;
    for (i = 0; i < siz; i++) {
        norm += a[i] * a[i];
    }
    return sqrt(norm);
}

//value in matrix cell 
double mesh(int ord, int size, int i, int j)
{
    double res = 0;
    switch (ord) {
    case 1: {
        res = (double)(size - max(i, j) + 1);
        break;
    }
    case 2: {
        res = (double)(max(i, j));
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
int matrform(int argc, double *M, int siz, int ord, char *filename)
{
    int i;
    int j;
    switch (argc) {
    case six: {
        FILE *fp;
        fp = fopen(filename, "r");
        if (!fp) {
            printf("The specified file does not exist\n");
            return -3;
        }

        fseek(fp, 0, SEEK_END);
        long pos = ftell(fp);
        rewind(fp);
        long leng = length(fp);
        if ((pos <= 0) || (leng == 0)) {
            printf("File is empty\n");
            fclose(fp);
            return -3;
        }
        if (leng != siz * siz) {
            printf("Wrong symbols in file\n");
            fclose(fp);
            return -3;
        }

        for (i = 0; i < siz; i++) {
            for (j = 0; j < siz; j++) {
                if (fscanf(fp, "%lf", &M[i * siz + j]) <= 0) {
                    printf("Failed to read matrix from file\n");
                    fclose(fp);
                    return -3;
                }
            }
        }
        fclose(fp);
        break;
    }
    case five: {
        for (i = 0; i < siz; i++) {
            for (j = 0; j < siz; j++) {
                M[i * siz + j] = mesh(ord, siz, i + 1, j + 1);
            }
        }
        break;
    }
    }
    return 0;
}

//printing of linear system
int matr_print(double *M, int l, int m, double *R, unsigned int fl_R)
{
    int k = m;
    int i;
    int j;
    if (k > l) {
        k = l;
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < k; j++) {
            printf("%10.3e", M[i * l + j]);
        }
        if (fl_R == 1) {
            printf("%10.3e", R[i]);
        }
        printf("\n");
    }
    printf("\n");
    return 0;
}

//wrapper function for clock_gettime for main process time
unsigned long long currentTimeNano_Pr()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

//wrapper function for clock_gettime for threads time
unsigned long long currentTimeNano_Thr()
{
    struct timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return (long long)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

//right column of linear system formation
void R_form(double *R, const double *M, int siz)
{
    int i;
    int j;
    for (i = 0; i < siz; i++) {
        R[i] = 0;
        for (j = 0; j <= (siz - 1) / 2; j++) {
            R[i] += M[i * siz + (2 * j)];
        }
    }
}

//infinity norm for answer vector
double norm_infinit(const double *M, int siz_real, int siz_cur)
{
    double max;
    int i;
    int j;
    max = 0.;
    for (i = 0; i < siz_cur; i++) {
        for (j = 0; j < siz_cur; j++) {
            if (fabs(M[i * siz_real + j]) > max) {
                max = fabs(M[i * siz_real + j]);
            }
        }
    }
    return max;
}

//Main thread function, used for distributed search of minimal element in row, distributed row changing and distributed residual computation
void *minus_spread(void *b)
{
	//receiving the shipment
    struct shipment *ship = (struct shipment *)(b);
    int size = ship->size;
    int threads = ship->threads;
    double *M = ship->M;
    double *R = ship->R;
    const double *proto_M = ship->proto_M;
    const double *proto_R = ship->proto_R;
    unsigned long long *times = ship->times;
    int *Transp = ship->Transp;
    double *Buffer_Transp = ship->Buffer_Transp;
    double epsi = ship->epsi;
    int iterator;
    int inner_iter;
    
    int numer = ship->numer;
    double koef;
    int t;

    double divid;

    double *Copy = ship->Copy;
    double *global_maxim = ship->global_maxim;
    pthread_barrier_t *bar = ship->bar;
    
    double Max_linie;
    int Max_p_linie;
    
    times[eigh * numer + 1] = currentTimeNano_Thr();
    times[eigh * numer + 2] = 0;
    times[eigh * numer + 3] = 0;
    times[eigh * numer + 4] = 0;
    times[eigh * numer + five] = 0;
    times[eigh * numer + six] = 0;
    times[eigh * numer + seven] = 0;
    times[eigh * numer + eigh] = 0;

    for (iterator = 0; iterator < size; ++iterator) {

        if (numer == 0) {
            Max_linie =  0.;
            Max_p_linie = 0;
            for(inner_iter = iterator; inner_iter < size; inner_iter++) {
                if (fabs(M[iterator * size + Transp[inner_iter]]) >
                    fabs(Max_linie)) {
                    Max_linie =
                        M[iterator * size + Transp[inner_iter]];
                    Max_p_linie = inner_iter;
                }
            }
            global_maxim[0] = Max_linie;
            if (Max_p_linie != iterator) {
            	t = Transp[iterator];
        	    Transp[iterator] = Transp[Max_p_linie];
        	    Transp[Max_p_linie] = t;
            }
        } 
        pthread_barrier_wait(&(bar[0]));
        
        if (fabs(global_maxim[0]) < epsi) {
        	if (numer == 0) {
        		ship->ende = -4;
        	}
            return NULL;
        }
            
        if (numer == 0) {
            memcpy(Copy, &M[iterator * size], size * sizeof(double));
            Copy[size] = R[iterator];
            for (inner_iter = size + 1; inner_iter < threads * (size + 1); inner_iter += size + 1) {
            	memcpy(&Copy[inner_iter], Copy, (size + 1) * sizeof(double));
            }
            
            for (inner_iter = size; inner_iter < threads * size; inner_iter += size) {
            	memcpy(&Transp[inner_iter], Transp, size * sizeof(int));
            }
        }

        pthread_barrier_wait(&(bar[0]));
        
        divid = 1 / Copy[numer * (size + 1) + Transp[numer * size + iterator]];

        for (inner_iter = numer; inner_iter < size; inner_iter += threads) {
            if (inner_iter == iterator) {
            	continue;
            }
            koef =
                -1 * M[inner_iter * size +
                    ship->Transp[numer * size + iterator]] * divid;
            for (t = iterator; t < size; t++) {
                M[inner_iter * size + Transp[numer * size + t]] +=
                    koef * Copy[numer * (size + 1) + Transp[numer * size + t]];
            }
            R[inner_iter] +=
                koef * Copy[numer * (size + 1) + size];
        }

        pthread_barrier_wait(&(bar[0]));

    }

    for (inner_iter = numer; inner_iter < size; inner_iter += threads) {
        R[inner_iter] = R[inner_iter] / M[inner_iter * size + Transp[numer * size + inner_iter]];
    }
    pthread_barrier_wait(&(bar[0]));

    for (inner_iter = numer; inner_iter < size; inner_iter += threads) {
        Buffer_Transp[inner_iter] = 0.;
        for (t = 0; t < size; t++) {
            Buffer_Transp[inner_iter] +=
                proto_M[inner_iter * size + t] * R[t];
            if (t == size - 1) {
                Buffer_Transp[inner_iter] -=
                    proto_R[inner_iter];
            }
        }
    }
    pthread_barrier_wait(&(bar[0]));
    //residual computation
    if (threads == 1) {
        (ship->residual1)[0] =
            norma(Buffer_Transp, size) / norma(proto_R, size);
        
    } else {
        if (numer == 0) {
            (ship->residual1)[0] = norma(Buffer_Transp, size);
        } else if (numer == 1) {
            (ship->residual2)[0] = norma(proto_R, size);
        }
        pthread_barrier_wait(&(bar[0]));
        if (numer == 0) {
            (ship->residual1)[0] = (ship->residual1)[0] / (ship->residual2)[0];
        }
    }

    if (numer == 0) {
        ship->ende = 0;
    }
    times[eigh * numer + 1] = currentTimeNano_Thr() - times[eigh * numer + 1];
    return NULL;
}

//solution of linear system by Gauss-Jordan elimination with selection of main element by row
int GJ_solution(int size, double *M, double *R, const double *proto_M,
                const double *proto_R, int threads, unsigned long long *times,
                double *Buffer_Transp, double *res1, double *res2,
                double *Copy, int* Transp)
{
    int i;
    int j;
    pthread_barrier_t bar;
    
    for (i = 0; i < size * size; i++) {
        M[i] = proto_M[i];
    }
    for (j = 0; j < size; j++) {
        R[j] = proto_R[j];
    }

    //formation of array of shipments
    struct shipment ship[threads];
    
    if (threads > 1) {
    	i = threads - 1;
    } else {
    	i = 1;
    }
    pthread_t thr[i];
    
    for (i = 0; i < threads; i++) {
    	ship[i].M = M;
        ship[i].R = R;
        ship[i].proto_M = proto_M;
        ship[i].proto_R = proto_R;
        ship[i].Transp = Transp;
        ship[i].Buffer_Transp = Buffer_Transp;
        ship[i].times = times;
        ship[i].residual1 = res1;
        ship[i].residual2 = res2;
        ship[i].Copy = Copy;
	    
	    ship[i].size = size;
        ship[i].threads = threads;
        ship[i].numer = i;
        ship[i].epsi = eps * norm_infinit(M, size, size);
        ship[i].global_maxim = res1;
        ship[i].bar = &bar;
    }
    
    pthread_barrier_init(&bar, NULL, threads);

    for (j = 0; j < size; j++) {
        Transp[j] = j;
    }

    for (j = 0; j < threads - 1; j++) {
        pthread_create(&thr[j], NULL, minus_spread, &(ship[j + 1]));
    }
    
    minus_spread(&ship[0]);

    for (j = 0; j < threads - 1; j++) {
        pthread_join(thr[j], NULL);
    }
    pthread_barrier_destroy (&bar);
    if (ship[0].ende == -4) {
        return -4;
    }
    return 0;
}

//error computation in main process
void error(const double *X, int siz)
{
    double *y;
    int i;
    y = (double *)malloc(siz * sizeof(double));
    for (i = 0; i < siz; i++) {
        if (i % 2 == 0) {
            y[i] = X[i] - 1;
        } else {
            y[i] = X[i];
        }
    }
    printf("Error: %10.3e\n", norma(y, siz));
    free(y);
}

