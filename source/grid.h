#ifndef _GRID_H_
#define _GRID_H_
#include "global.h"
#include "input.h"
class grid
{
public:
    int lx, ly, lz;
    int nx;
    int ny;
    int nz;
    int n_fixed;
    double d;                                     // side length of blocks, if it's strict uniform(cube).
    double d1,d2,d3;//uniform on one dimension, cuboid for meshes.
    double x[55] = {0}, y[55] = {0}, z[55] = {0}; // fixed points
    double rdf_cutoff;
    double dr;  // rdf Discretization distance
    int mesh;   // rdf Discretization number
    int l;      // Azimuthal quantum number,1 by default.
    double *gd; // 3D grid
    double *rdf;
    bool good_grid_condition();
    bool good_fixed_points(int n_fixed);
    //void generate_uniform_grid(ofstream &fout);
    void read_fixed_points(const string &file_points, ofstream &fout);
    void read_venergy_and_construct(input& I1 ,ofstream &fout);
    void read_distribution(const string &file_distribution, ofstream &fout);
    void integral_as_sum();
    void rdf_para_init();
    bool rdf_para_check();
};

#endif