/* vecmatd.h - Routinen zur Vektor- und Matrixrechnung */
/*             double precision                        */
/* Arbeiten mit Arraydarstellung aus Numerical Recipes */
/* Volker Ahlers, November 1998                        */

void print_vec(double *x,int dim);
void print_mat(double **a,int dim1,int dim2);
void set_mat_id(double **a,int dim);
void copy_vec(double *y,double *x,int dim);
void copy_mat(double **b,double **a,int dim1,int dim2);
double skalprod(double *x,double *y,int dim);
double norm_vec(double *x,int dim);
double renorm_vec(double *x,int dim);
void mult_mat_vec(double *y,double **a,double *x,int dim1,int dim2);
void mult_mat_mat(double **c,double **a,double **b,int dim1,int dim2);
void cgs(double **q,double *qtemp,double *r,int m,int n);
void mgsr(double **q,double *r,int m,int n);
void mgsr_sum(double **q,double *r,int m,int n);
