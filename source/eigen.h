#ifndef _EIGEN_H_
#define _EIGEN_H_
#include "global.h"
#include "grid.h"
class eigen{
public:
    double* H;
    double integral_as_sum(grid& G1,int i,int j);
    void Hsolver();
    void output_finalize(grid&G1);
    void Hconstructor(grid& G1);
    double distance(double x1,double y1,double z1,double x2,double y2,double z2);
    double calc_diagonal_element(int i);
    double rdf(grid& G1,int i);
};


#endif