#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#define five 5
#include "func.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define GIGA_MODIFIER 1e9
#define eps 1e-10

//filling a matrix with zeros
void Y0_form(double *Y, int N1, int N2)
{
    int i, j;

    double h = 1. / N1;
    double iter;
 
    for (i = 1, iter = h; i <= N1 - 1; i++, iter += h) {
        for (j = 0; j <= N2 - 1; j++) {
            Y[(i - 1) * N2 + j] = 0.;
        }
    }
}

//value of laplacian of choosen function at the point
double answ(int ord, double rad, double grad)
{
    switch (ord)
    {
    case(3):
        return (-3. * pow(rad, 5.) * cos(5. * grad));
    case(4):
        return pow(rad, 2.);
    case(5):
        return (pow(rad, 2.) * sin(grad));
    case(6):
        return (pow(rad, 3.) * sin(grad) + 9.);
    case(7):
        return (pow(rad, 4.) * pow(cos(grad), 2.));
    case(9):
        return ((pow(rad, 2.) - pow(rad, 3.)) * cos(grad) + 9.);
    case(10):
        return ((pow(rad, 2.) - pow(rad, 3.)) * cos(grad));
    case(13):
        return (exp( pow(rad, 4.) + 7. * pow(rad, 3.) + 2.));
    }
}

//value of choosen boundary condition at the point
double alpha(int ord, double grad)
{
    double x = cos(grad);
    double y = sin(grad);

    switch (ord)
    {
    case(1):
        return 7.;
    case(2):
        return (pow(x, 2.) + pow(y, 2.));
    case(3):
        return ((-3) * cos(5. * grad));
    case(4):
        return 1.;
    case(5):
        return sin(grad);
    case(6):
        return sin(grad) + 9.;
    case(7):
        return (4. * pow(cos(grad), 2.));
    case(9):
        return 9.;
    case(10):
        return 0.;
    case(13):
        return exp(10.);
    }
}

//value of choosen function at the point
double f(int ord, double rad, double grad)
{
    double x = rad * cos(grad);
    double y = rad * sin(grad);
    
    switch (ord)
    {
    case(1):
        return 4.;
    case(2):
        return (pow(x, 2.) + pow(y, 2.));
    case(3):
        return 0.;
    case(4):
        return -4.;
    case(5):
        return (-3. * sin(grad));
    case(6):
        return (-8. * rad * sin(grad));
    case(7):
        return (-2 * pow(rad, 2.) * (7 * pow(cos(grad), 2.) + pow(sin(grad), 2.)));
    case(9):
        return ( (9. * rad - 4.) * cos(grad) + (1. - rad) * sin(grad) );
    case(10):
        return ((9. * rad - 4.) * cos(grad) + (1. - rad) * sin(grad));
    case(13):
        return -(answ(ord, rad, grad) * (rad * pow(4. * pow(rad, 3.) + 21. * pow(rad, 2.), 2.) + 16. * pow(rad, 3.) + 63. * pow(rad, 2.)));
    }
}

//infinite measure between real function and approximate values
double norma(double* Y, int ord, int N1, int N2, double h1, double h2, double y0)
{
    double iter1, iter2;
    int i, j;
    double puff;
    double ret = fabs(answ(ord, 0., 0.) - y0);
    double m_rad = 0.;
    double m_grad = 0.;
    int i_m = 0;
    int j_m = 0;

    for (i = 1, iter1 = h1; i <= N1 - 1; i++, iter1 += h1) {
        for (j = 0, iter2 = 0.; j <= N2 - 1; j++, iter2 += h2) {
            puff = fabs(answ(ord, iter1, iter2) - Y[(i - 1) * N2 + j]);
            if (puff > ret) {
                ret = puff;
                m_rad = iter1;
                m_grad = iter2;
                i_m = i - 1;
                j_m = j;
            }
        }
    }
    printf("max = %f, m_rad = %f, m_grad = %f, i_m = %d, j_m = %d\n\n", ret, m_rad, m_grad, i_m, j_m);
    return ret;
}

//formation of F matrix
void f_form(int ord_f, int ord_alpha, double* V, int N2, int N1, double h1, double h2)
{
    double iter_1;
    double iter_2;
    int i, j;

    for (i = 1, iter_1 = h1; i <= N1 - 2; i++, iter_1 += h1) {
        for (j = 0, iter_2 = 0.; j <= N2 - 1; j++, iter_2 += h2) {
            V[(i - 1) * N2 + j] = pow(iter_1, 2.) * f(ord_f, iter_1, iter_2);
            //printf("ro^2 = %f, F = %f\n", pow(iter_1, 2.), V[(i - 1) * N2 + j]);
        }
    }

    for (j = 0, iter_2 = 0.; j <= N2 - 1; j++, iter_2 += h2) {
        V[(i - 1) * N2 + j] = pow(iter_1, 2.) * f(ord_f, iter_1, iter_2) + (iter_1 * (iter_1 + 0.5 * h1) * alpha(ord_alpha, iter_2)) / pow(h1, 2.);
        //printf("iter_1 = %f, iter_2 = %f, F = %f\n", iter_1, iter_2, V[(i - 1) * N2 + j]);
    }
}

void A1_form(double* A1, int N1, double h1)
{
    int i = 1;
    double iter = h1;
    A1[1] = (2. * pow(iter, 2.)) / pow(h1, 2.);
    A1[2] = ((-1.) * iter * (iter + 0.5 * h1)) / pow(h1, 2.);
    i++;
    iter += h1;

    for (; i <= N1 - 2; i++, iter += h1) {
        A1[(i - 1) * 3] = ((-1.) * iter * (iter - 0.5 * h1)) / pow(h1, 2.);
        A1[(i - 1) * 3 + 1] = (2. * pow(iter, 2.)) / pow(h1, 2.);
        A1[(i - 1) * 3 + 2] = ((-1.) * iter * (iter + 0.5 * h1)) / pow(h1, 2.);
    }

    A1[(i - 1) * 3] = ((-1.) * iter * (iter - 0.5 * h1)) / pow(h1, 2.);
    A1[(i - 1) * 3 + 1] = (2. * pow(iter, 2.)) / pow(h1, 2.);

}

//multiplying A1 by V to X
void A1_mult(double* A1, double* V, double* X, int N2, int N1)
{
    int i, j;

    for (j = 0; j <= N2 - 1; j++) {
        X[j] = -(A1[1] * V[j] + A1[2] * V[N2 + j]);
        //printf("Ord_1: X = %f\n", X[j]);
    }
    for (i = 2; i <= N1 - 2; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            X[(i - 1) * N2 + j] = -(A1[(i - 1) * 3] * V[(i - 2) * N2 + j] + A1[(i - 1) * 3 + 1] * V[(i - 1) * N2 + j] + A1[(i - 1) * 3 + 2] * V[i * N2 + j]);
            //printf("Ord_2: X = %f\n", X[(i - 1) * N2 + j]);
        }
    }
    for (j = 0; j <= N2 - 1; j++) {
        X[(i - 1) * N2 + j] = -(A1[(i - 1) * 3] * V[(i - 2) * N2 + j] + A1[(i - 1) * 3 + 1] * V[(i - 1) * N2 + j]);
        //printf("Ord_3: X = %f\n", X[(i - 1) * N2 + j]);
    }
}

//multiplying A2 by V to X
void A2_mult(double A2_1, double A2_2, double* V, double* X, int N2, int N1)
{
    int i, j;
    for (i = 1; i <= N1 - 1; i++) {
        X[(i - 1) * N2] = -(A2_2 * (V[(i - 1) * N2 + N2 - 1] + V[(i - 1) * N2 + 1]) + A2_1 * V[(i - 1) * N2]);
        for (j = 1; j <= N2 - 2; j++) {
            X[(i - 1) * N2 + j] = -(A2_2 * (V[(i - 1) * N2 + j - 1] + V[(i - 1) * N2 + j + 1]) + A2_1 * V[(i - 1) * N2 + j]);
        }
        X[(i - 1) * N2 + N2 - 1] = -(A2_2 * (V[(i - 1) * N2 + N2 - 2] + V[(i - 1) * N2]) + A2_1 * V[(i - 1) * N2 + N2 - 1]);
    }
}

//coordinate sum of vectors
void vect_summ(double* X, double* V, int N2, int N1)
{
    int i, j;
    for (i = 1; i <= N1 - 1; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            X[(i - 1) * N2 + j] += V[(i - 1) * N2 + j];
        }
    }
}

//multiplying vector X by number k
void vect_koef(double* X, double k, int N2, int N1)
{
    int i, j;
    for (i = 1; i <= N1 - 1; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            X[(i - 1) * N2 + j] *= k;
        }
    }
}

//copying vector Y0 to Y
void vect_copy(double* Y, double* Y0, int N2, int N1)
{
    int i, j;
    for (i = 1; i <= N1 - 1; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            Y[(i - 1) * N2 + j] = Y0[(i - 1) * N2 + j];
        }
    }
}

//euclid measure between vectors
double vect_norm(double* X, double* Y, int N2, int N1) {
    int i, j;
    double puff = 0.;

    for (i = 0; i <= N1 - 2; i++) {
        for (j = 0; j <= N2 - 1; j++) {
            puff += pow(X[i * N2 + j] - Y[i * N2 + j], 2.);
        }
    }
    return (sqrt(puff));
}

//max coordinate of vector
double vmax(int h_gr, double* V)
{
    double res = V[0];
    for (int i = 1; i <= h_gr; i++) {
        if (V[i] > res) {
            res = V[i];
        }
    }
    return res;
}

int first_system(double* A1, double k, double* R, double* X, int N2, int N1, double* M_prog, double* R_prog, double* X_prog, double* V_prog, double* U_prog)
{
    int i, j;

    M_prog[1] = A1[1] + k;
    M_prog[2] = A1[2];

    for (i = 2; i <= N1 - 2; i++) {
        M_prog[(i - 1) * 3] = A1[(i - 1) * 3];
        M_prog[(i - 1) * 3 + 1] = A1[(i - 1) * 3 + 1] + k;
        M_prog[(i - 1) * 3 + 2] = A1[(i - 1) * 3 + 2];
    }

    M_prog[(i - 1) * 3] = A1[(i - 1) * 3];
    M_prog[(i - 1) * 3 + 1] = A1[(i - 1) * 3 + 1] + k;

    for (j = 0; j <= N2 - 1; j++) {
        R_prog[0] = R[j];

        for (i = 2; i <= N1 - 2; i++) {
            R_prog[i - 1] = R[(i - 1) * N2 + j];
        }

        R_prog[i - 1] = R[(i - 1) * N2 + j];

        if (thomas_algo(N1 - 2, M_prog, R_prog, X_prog, U_prog, V_prog) == -4) {
            return -4;
        }
        for (i = 1; i <= N1 - 1; i++) {
            X[(i - 1) * N2 + j] = X_prog[i - 1];
        }
    }
    return 0;
}

int second_system(double A2_1, double A2_2, double k, double* R, double* X, int N2, int N1, double* M_prog, double* R_prog, double* X1_prog, double* X2_prog, double* V_prog, double* U_prog)
{
    int i, j;
    double y0;

    M_prog[1] = A2_1 + k;
    M_prog[2] = A2_2;

    for (j = 2; j <= N2 - 2; j++) {
        M_prog[(j - 1) * 3] = A2_2;
        M_prog[(j - 1) * 3 + 1] = A2_1 + k;
        M_prog[(j - 1) * 3 + 2] = A2_2;
    }
    M_prog[(j - 1) * 3] = A2_2;
    M_prog[(j - 1) * 3 + 1] = A2_1 + k;

    for (i = 1; i <= N1 - 1; i++) {
    
        //R_prog[0] = R[(i - 1) * N2 + 1];

        for (j = 1; j <= N2 - 1; j++) {
            R_prog[j - 1] = R[(i - 1) * N2 + j];
        }
        //R_prog[j - 1] = R[(i - 1) * N2 + j];

        if (thomas_algo(N2 - 2, M_prog, R_prog, X1_prog, U_prog, V_prog) == -4) {
            return -4;
        }

        R_prog[0] = (-1.) * A2_2;

        for (j = 2; j <= N2 - 2; j++) {
            R_prog[j - 1] = 0.;
        }
        
        R_prog[j - 1] = (-1.) * A2_2;

        if (thomas_algo(N2 - 2, M_prog, R_prog, X2_prog, U_prog, V_prog) == -4) {
            return -4;
        }

        y0 = (R[(i - 1) * N2] - A2_2 * (X1_prog[0] + X1_prog[N2 - 2])) / (A2_1 + k + A2_2 * (X2_prog[0] + X2_prog[N2 - 2]));
        
        X[(i - 1) * N2] = y0;
        for (j = 1; j <= N2 - 1; j++) {
            X[(i - 1) * N2 + j] = X1_prog[j - 1] + y0 * X2_prog[j - 1];
        }
    }
    return 0;
}

//tridiagonal matrix algorithm
int thomas_algo(int h_gr, double* M, double* R, double* X, double* U, double* V)
{
    double puffer;
    int i;

    if (fabs(M[1]) < eps) {
        printf("%d -> %f\n", 0, M[1]);
        return -4;
    }
    V[0] = (-1.) * M[2] / M[1];
    U[0] = R[0] / M[1];

    for (i = 1; i < h_gr; i++) {
        puffer = (-1.) * M[i * 3 + 1] - M[i * 3 + 0] * V[i - 1];
        if (fabs(puffer) < eps) {
            printf("%d -> %f\n", i, puffer);
            return -4;
        }
        V[i] = M[i * 3 + 2] / puffer;
        U[i] = (M[i * 3 + 0] * U[i - 1] - R[i]) / puffer;
    }

    puffer = (-1.) * M[h_gr * 3 + 1] - M[h_gr * 3 + 0] * V[h_gr - 1];
    if (fabs(puffer) < eps) {
        printf("%d -> %f\n", h_gr, puffer);
        return -4;
    }
    V[h_gr] = 0.;
    U[h_gr] = (M[h_gr * 3 + 0] * U[h_gr - 1] - R[h_gr]) / puffer;


    X[h_gr] = U[h_gr];
    for (i = h_gr - 1; i >= 0; i--) {
        X[i] = V[i] * X[i + 1] + U[i];
    }

    return 0;
}

int zircle_thomas_algo(int h_gr, double* M, double* R, double* X, double* X1, double* X2, double* U, double* V, double* W)
{
    int i;
    double puff, yN;

    puff = M[1];
    V[0] = M[2] / puff;
    U[0] = R[0] / puff;
    W[0] = M[0] / puff;

    for (i = 1; i <= h_gr; i++) {
        puff = M[i * 3 + 1] - V[i - 1] * M[i * 3];
        V[i] = M[i * 3 + 2] / puff;
        U[i] = (R[i] + M[i * 3] * U[i - 1]) / puff;
        W[i] = (M[i * 3] * W[i - 1]) / puff;

    }

    X1[h_gr - 1] = U[h_gr - 1];
    X2[h_gr - 1] = V[h_gr - 1] + W[h_gr - 1];
    printf("%d -> X1: %f, X2: %f\n", h_gr - 1, X1[h_gr - 1], X2[h_gr - 1]);
    for (i = h_gr - 2; i >= 0; i--) {
        X1[i] = V[i + 1] * X1[i] + U[i];
        X2[i] = V[i + 1] * X2[i] + W[i];
        printf("%d -> X1: %f, X2: %f\n", i, X1[i], X2[i]);
    }
    return 0;
}
