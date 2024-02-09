#define _CRT_SECURE_NO_WARNINGS

#include "func.h"
#include <errno.h>
#include <fenv.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define GIGA_MODIFIER 1e9
#define NANO_MODIFIER 1e-6
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main()
{
    int ord_f = 6;
    int ord_alpha = ord_f;
    int N1, N2;
    double h1, h2, d1, D1, d2, D2, m1, m2;
    double* M_prog;
    double* R_prog;
    double* X_prog;
    double* U_prog;
    double* V_prog;

    //int n;
    double a, nu, nu_s, epsi, epsi_s, b, t, r, s, q;
    double delta, hi;
    //double* W;
    //double* Mu;

    double* A1;
    double A2_1, A2_2;
    double* F;
    double* Y;
    double* X;
    double* V;

    double* M_prog2;
    double* R_prog2;
    double* X1_prog2;
    double* X2_prog2;
    double* U_prog2;
    double* V_prog2;

    double y0;

    int i, gr;
    double iter;
    
    //setting approximation grid parameters

    printf("N1 = ");
    scanf("%d", &N1);
    h1 = 1. / N1;

    printf("N2 = ");
    scanf("%d", &N2);
    h2 = (2. * M_PI) / N2;
    
    //memory allocation

    if (!(M_prog = (double*)malloc((N1 - 1) * 3 * sizeof(double)))) {
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(R_prog = (double*)malloc((N1 - 1) * sizeof(double)))) {
        free(M_prog);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(X_prog = (double*)malloc((N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(U_prog = (double*)malloc((N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(V_prog = (double*)malloc((N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        printf("Failed to allocate memory\n");
        return -2;
    }

    //calculation of bounds of operators
    d2 = 0.; D2 = 4. / pow(h2, 2);

    //lower bound d1
    i = 1;
    iter = h1;
    M_prog[1] = (-2.) * pow(iter, 2.);
    M_prog[2] = iter * (iter + 0.5 * h1);
    R_prog[0] = (-1.) * pow(h1, 2.);
    i++;
    iter += h1;

    for (; i <= N1 - 2; i++, iter += h1) {
        M_prog[(i - 1) * 3] = iter * (iter - 0.5 * h1);
        M_prog[(i - 1) * 3 + 1] = (-2.) * pow(iter, 2.);
        M_prog[(i - 1) * 3 + 2] = iter * (iter + 0.5 * h1);
        R_prog[i - 1] = (-1.) * pow(h1, 2.);
    }

    M_prog[(i - 1) * 3] = iter * (iter - 0.5 * h1);
    M_prog[(i - 1) * 3 + 1] = (-2.) * pow(iter, 2.);
    R_prog[i - 1] = (-1.) * pow(h1, 2.);

    if (thomas_algo(N1 - 2, M_prog, R_prog, X_prog, U_prog, V_prog) == -4) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        return -4;
    }

    d1 = 1. / fmax(vmax(N1 - 2, X_prog), 0.);

    //upper bound D1

    i = 1;
    iter = h1;
    M_prog[1] = -(2. * pow(iter, 2.) + pow(h1, 2.));
    M_prog[2] = iter * (iter + 0.5 * h1);
    R_prog[0] = 0.;
    i++;
    iter += h1;

    for (; i <= N1 - 2; i++, iter += h1) {
        M_prog[(i - 1) * 3] = iter * (iter - 0.5 * h1);
        M_prog[(i - 1) * 3 + 1] = -(2. * pow(iter, 2.) + pow(h1, 2.));
        M_prog[(i - 1) * 3 + 2] = iter * (iter + 0.5 * h1);
        R_prog[i - 1] = 0.;
    }

    M_prog[(i - 1) * 3] = iter * (iter - 0.5 * h1);
    M_prog[(i - 1) * 3 + 1] = -(2. * pow(iter, 2.) + pow(h1, 2.));
    R_prog[i - 1] = 0.;

    if (thomas_algo(N1 - 2, M_prog, R_prog, X_prog, U_prog, V_prog) == -4) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        return -4;
    }

    m1 = fmax(0., vmax(N1 - 2, X_prog));
    m2 = fmax(((2. - h1)*(4. - h1)) / (4. * pow(h1, 2.)), fmax(0.5, (4. * pow(1 - h1, 2.)) / pow(h1, 2.)));
    D1 = m1 + m2 * (1. + m1);

    //iteration parameters

    a = sqrt(( (D1 - d1) * (D2 - d2) ) / ( (D1 + d2) * (D2 + d1) ));
    nu = (1. - a) / (1. + a);

    b = (a * (D2 + d1)) / (D1 - d1);
    t = (1. - b) / (1. + b);
    r = (D2 + D1 * b) / (1. + b);
    s = (D2 - D1 * b) / (1. + b);

    //2nd memory allocation

    if (!(A1 = (double*)malloc(3 * (N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(Y = (double*)malloc(N2 * (N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(X = (double*)malloc(N2 * (N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(V = (double*)malloc(N2 * (N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(F = (double*)malloc(N2 * (N1 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(M_prog2 = (double*)malloc(3 * (N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);

        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(R_prog2 = (double*)malloc((N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(X1_prog2 = (double*)malloc((N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        free(R_prog2);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(X2_prog2 = (double*)malloc((N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        free(R_prog2);
        free(X1_prog2);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(U_prog2 = (double*)malloc((N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        free(R_prog2);
        free(X1_prog2);
        free(X2_prog2);
        printf("Failed to allocate memory\n");
        return -2;
    }
    if (!(V_prog2 = (double*)malloc((N2 - 1) * sizeof(double)))) {
        free(M_prog);
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        free(R_prog2);
        free(X1_prog2);
        free(X2_prog2);
        free(U_prog2);
        printf("Failed to allocate memory\n");
        return -2;
    }

    double* Y0 = (double*)malloc(N2 * (N1 - 1) * sizeof(double));

    //creating starting matrices
    A1_form(A1, N1, h1);
    A2_1 = 2. / pow(h2, 2.);
    A2_2 = -1. / pow(h2, 2.);
    f_form(ord_f, ord_alpha, F, N2, N1, h1, h2);
    Y0_form(Y, N1, N2);

    double W1, W2;

    //constants W1, W2 from common case

    W1 = (r * sqrt(nu) + s) / (1. + t * sqrt(nu));
    W2 = (r * sqrt(nu) - s) / (1. - t * sqrt(nu));

    int puff;
    int j;

    //solution by the method of variable directions

    printf("Starting norm = %f\n\n", norma(Y, ord_f, N1, N2, h1, h2, 0.));

    int it;
    it = 0;
    epsi = NANO_MODIFIER;
    double stop;

    do {

        A2_mult(A2_1, A2_2, Y, X, N2, N1);
        vect_koef(Y, W1, N2, N1);
        vect_summ(Y, X, N2, N1);
        vect_summ(Y, F, N2, N1);

        first_system(A1, W1, Y, X, N2, N1, M_prog, R_prog, X_prog, V_prog, U_prog);

        A1_mult(A1, X, Y, N2, N1);
        vect_koef(X, W2, N2, N1);
        vect_summ(X, Y, N2, N1);
        vect_summ(X, F, N2, N1);

        second_system(A2_1, A2_2, W2, X, Y, N2, N1, M_prog2, R_prog2, X1_prog2, X2_prog2, V_prog2, U_prog2); 

        it++;

        A1_mult(A1, Y, X, N2, N1);
        A2_mult(A2_1, A2_2, Y, Y0, N2, N1);
        vect_summ(X, Y0, N2, N1);
        vect_koef(F, -1., N2, N1);
        stop = vect_norm(X, F, N2, N1);

        printf("SLU-Delta of %d = %.11f\n", it, stop);
        vect_koef(F, -1., N2, N1);
    } while (stop >= epsi);

    A1_mult(A1, Y, X, N2, N1);
    A2_mult(A2_1, A2_2, Y, Y0, N2, N1);
    vect_summ(X, Y0, N2, N1);
    vect_koef(F, -1., N2, N1);
    printf("SLU-Delta = %.11f\n", vect_norm(X, F, N2, N1));

    //calculation of value in zero point

    i = 1;
    iter = h1;
    M_prog[1] = -2. * iter;
    M_prog[2] = iter + 0.5 * h1;
    R_prog[0] = -0.5 * h1;
    i++;
    iter += h1;

    for (; i <= N1 - 2; i++, iter += h1) {
        M_prog[(i - 1) * 3] = iter - 0.5 * h1;
        M_prog[(i - 1) * 3 + 1] = -2. * iter;
        M_prog[(i - 1) * 3 + 2] = iter + 0.5 * h1;
        R_prog[i - 1] = 0.;
    }

    M_prog[(i - 1) * 3] = iter - 0.5 * h1;
    M_prog[(i - 1) * 3 + 1] = -2. * iter;
    R_prog[i - 1] = 0.;

    if (thomas_algo(N1 - 2, M_prog, R_prog, X_prog, U_prog, V_prog) == -4) {
        free(R_prog);
        free(X_prog);
        free(U_prog);
        free(V_prog);
        free(A1);
        free(Y);
        free(X);
        free(V);
        free(F);
        free(M_prog2);
        free(R_prog2);
        free(X1_prog2);
        free(X2_prog2);
        free(U_prog2);
        free(V_prog2);
        return -4;
    }
    
    printf("Norma0 = %.11f\n", norma(Y, ord_f, N1, N2, h1, h2, answ(ord_f, 0., 0.)));

    double schet = 0.;
    for (i = 0; i <= N2 - 1; i++) {
        schet += Y[i];
    }

    y0 = -(M_PI * pow(h1, 2.) * f(ord_f, 0., 0.) + 2. * h2 * schet) / (2. * h2 * N2 * (X_prog[0] - 1.) );

    printf("Zero-norm = %f\n\n", fabs(y0 - answ(ord_f, 0., 0.)));

    //final approximation for point (0, 0)
    for (i = 1; i <= N1 - 1; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            Y[(i - 1) * N2 + j] = Y[(i - 1) * N2 + j] + y0 * X_prog[i - 1];
        }
    }

    printf("Norma = %.11f\n", norma(Y, ord_f, N1, N2, h1, h2, y0));

    //memory freeing up
    free(M_prog);
    free(R_prog);
    free(X_prog);
    free(U_prog);
    free(V_prog); 
    //free(W);
    free(A1);
    free(Y);
    free(X);
    free(V);
    free(F);
    free(M_prog2);
    free(R_prog2);
    free(X1_prog2);
    free(X2_prog2);
    free(U_prog2);
    free(V_prog2);
    return 0;
}
