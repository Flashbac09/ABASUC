#include "global.h"
#include "input.h"
#include "grid.h"
#include "eigen.h"
int main()
{
    preset_overall();
    input I1;
    grid G1;
    eigen E1;
    I1.input_initialize("INPUT.txt", fout);
    G1.read_venergy_and_construct(I1, fout);
    G1.read_fixed_points(I1.points_path, fout);
    G1.read_distribution(I1.distribution_path, fout);
    E1.spline_interpolation(G1);
    E1.Hconstructor(G1);
    E1.Hsolver();
    E1.output_finalize(G1);

    //new
    // 1. LCAO intro
    // 2. V/rdf/points figure
    // 3. test list
    // 4. final report
}