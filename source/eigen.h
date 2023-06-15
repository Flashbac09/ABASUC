#ifndef _EIGEN_H_
#define _EIGEN_H_
#include "global.h"
#include "grid.h"
// check if you have lapack compiled and installed. //
#include <lapacke.h>
class eigen{
public:
    double* H;
    int n_fixed;
    double w[MAX_SPHERE];//for eigenvalues
    double integral_as_sum(grid& G1,int i,int j);
    void Hsolver();
    void output_finalize(grid&G1);
    void Hconstructor(grid& G1);
    double distance(double x1,double y1,double z1,double x2,double y2,double z2);
    double calc_diagonal_element(int i);
    double erdf(grid& G1,int i);
    void spline_interpolation(grid &G1);
    int  mesh;
    double *y;
    double* m;//for spline interpolation parameter
    double drdf(double x,double h);//rdf from distancce
};


#endif
