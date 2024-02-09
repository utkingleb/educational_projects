#ifndef FUNC_H
#define FUNC_H

int thomas_algo(int h_gr, double* M, double* R, double* X, double* U, double* V);
int zircle_thomas_algo(int h_gr, double* M, double* R, double* X, double* X1, double* X2, double* U, double* V, double* W);

double vmax(int h_gr, double* V);

double answ(int ord, double rad, double grad);
double alpha(int ord, double grad);
double f(int ord, double rad, double grad);
double norma(double* Y, int ord, int N1, int N2, double h1, double h2, double y0);
void f_form(int ord_f, int ord_alpha, double* V, int N2, int N1, double h1, double h2);

void A1_form(double* A1, int N1, double h1);
void A1_mult(double* A1, double* V, double* X, int N2, int N1);
void A2_mult(double A2_1, double A2_2, double* V, double* X, int N2, int N1);
void vect_summ(double* X, double* V, int N2, int N1);
void vect_koef(double* X, double k, int N2, int N1);
double vect_norm(double* X, double* Y, int N2, int N1);
void vect_copy(double* Y, double* Y0, int N2, int N1);

int first_system(double* A1, double k, double* R, double* X, int N2, int N1, double* M_prog, double* R_prog, double* X_prog, double* V_prog, double* U_prog);
int second_system(double A2_1, double A2_2, double k, double* R, double* X, int N2, int N1, double* M_prog, double* R_prog, double* X1_prog, double* X2_prog, double* V_prog, double* U_prog);

void Y0_form(double* Y, int N1, int N2);

#endif