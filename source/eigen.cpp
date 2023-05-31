#include "eigen.h"
inline double eigen::rdf(grid& G1,int i)
{
    return 1;
}
inline double eigen::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
}
inline double eigen::integral_as_sum(grid &G1, int i, int j) // inline/parallel? //
{
    /*
    if(i==j)
    {
        return calc_diagonal_element(int i);
    }
    */
    double dis = sqrt(pow((G1.x[i] - G1.x[j]), 2) + pow((G1.y[i] - G1.y[j]), 2) + pow((G1.z[i] - G1.z[j]), 2));
    if (dis > 2 * G1.rdf_cutoff)
        return 0;
    cout<<"dis: "<< dis<<endl;
    // 1. determine overlap center
    double Cx = 0.5 * (G1.x[i] + G1.x[j]), Cy = 0.5 * (G1.y[i] + G1.y[j]), Cz = 0.5 * (G1.z[i] + G1.z[j]);
    // 2. determine overlap cuboid x,y,z
    double Superior_x, Inferior_x, Superior_y, Inferior_y, Superior_z, Inferior_z;
    Superior_x = max(G1.x[i], G1.x[j]);
    Inferior_x = min(G1.x[i], G1.x[j]);
    Superior_y = max(G1.y[i], G1.y[j]);
    Inferior_y = min(G1.y[i], G1.y[j]);
    Superior_z = max(G1.z[i], G1.z[j]);
    Inferior_z = min(G1.z[i], G1.z[j]);
    cout<<"Superior/Inferior x,y,z for "<<i<<"  "<<j<<endl;
    cout<<Superior_x<<"  "<<Inferior_x<<endl;
    cout<<Superior_y<<"  "<<Inferior_y<<endl;
    cout<<Superior_z<<"  "<<Inferior_z<<endl;
    /*
    if (Superior_x - Inferior_x >= 2 * G1.rdf_cutoff || Superior_y - Inferior_y >= 2 * G1.rdf_cutoff || Superior_z - Inferior_z >= 2 * G1.rdf_cutoff)
        return 0; // for stability
    */
    double forward_x, forward_y, forward_z, back_x, back_y, back_z;
    forward_x = min(Inferior_x + G1.rdf_cutoff, (double)G1.lx); // boundary
    back_x = max(Superior_x - G1.rdf_cutoff, 0.0);
    forward_y = min(Inferior_y + G1.rdf_cutoff, (double)G1.ly);
    back_y = max(Superior_y - G1.rdf_cutoff, 0.0);
    forward_z = min(Inferior_z + G1.rdf_cutoff, (double)G1.lz);
    back_z = max(Superior_z - G1.rdf_cutoff, 0.0);
    double Overlap_x, Overlap_y, Overlap_z;
    Overlap_x = forward_x - back_x;
    Overlap_y = forward_y - back_y;
    Overlap_z = forward_z - back_z;
    cout<<"Back/Forward x,y,z for "<<i<<"  "<<j<<endl;
    cout<<back_x<<"  "<<forward_x<<endl;
    cout<<back_y<<"  "<<forward_y<<endl;
    cout<<back_z<<"  "<<forward_z<<endl;
    /*
    if (Overlap_x <= 0 || Overlap_y <= 0 || Overlap_z <= 0)
        return 0;//for stability
    */
    // 3. go through grid and sum
    int x_start,y_start,z_start,x_end,y_end,z_end;
    //here check if points's positions are special(marginal part)
    x_start=(int)(back_x/G1.d1);
    y_start=(int)(back_y/G1.d2);
    z_start=(int)(back_z/G1.d3);
    x_end=min(((int)(forward_x/G1.d1)+1),G1.nx-1);
    y_end=min(((int)(forward_y/G1.d2)+1),G1.ny-1);
    z_end=min(((int)(forward_z/G1.d3)+1),G1.nz-1);
    cout<<"start/end of x,y,z for "<<i<<"  "<<j<<endl;
    cout<<x_start<<"  "<<x_end<<endl;
    cout<<y_start<<"  "<<y_end<<endl;
    cout<<z_start<<"  "<<z_end<<endl;
    double res=0;
    for(int m=x_start;m<=x_end;++m)
    {
        for(int n=y_start;n<=y_end;++n)
        {
            for(int k=z_start;k<=z_end;++k)
            {
                //cout<<m<<"  "<<n<<"  "<<k<<"  "<<G1.gd[m*G1.ny*G1.nz+n*G1.nz+k]<<endl;
                res+=G1.gd[m*G1.ny*G1.nz+n*G1.nz+k]*rdf(G1,i)*rdf(G1,j)*G1.d1*G1.d2*G1.d3;
            }
        }
    }
    return res;
    /*
    double A[6] = {0}, B[6] = {0}, mid[3] = {0}, LA[3] = {0}, LB[3] = {0};
    A[0] = G1.x[i] - G1.rdf_cutoff;
    A[1] = G1.y[i] - G1.rdf_cutoff;
    A[2] = G1.z[i] - G1.rdf_cutoff;
    A[3] = G1.x[i] + G1.rdf_cutoff;
    A[4] = G1.y[i] + G1.rdf_cutoff;
    A[5] = G1.z[i] + G1.rdf_cutoff;
    B[0] = G1.x[j] - G1.rdf_cutoff;
    B[1] = G1.y[j] - G1.rdf_cutoff;
    B[2] = G1.z[j] - G1.rdf_cutoff;
    B[3] = G1.x[j] + G1.rdf_cutoff;
    B[4] = G1.y[j] + G1.rdf_cutoff;
    B[5] = G1.z[j] + G1.rdf_cutoff;
    if (distance(A[0], A[1], A[2], B[3], B[4], B[5]) < distance(B[0], B[1], B[2], A[3], A[4], A[5]))
    {
        B[0] = B[3]; // A--- . B+++
        B[1] = B[4];
        B[2] = B[5];
    }
    else
    {
        A[0] = A[3]; // A+++ . B---
        A[1] = A[4];
        A[2] = A[5];
    }
    */
}
void eigen::Hsolver()
{
    // link
}
void eigen::output_finalize(grid &G1)
{
    delete[] H;
    delete[] G1.rdf;
    delete[] G1.gd;
}
void eigen::Hconstructor(grid &G1)
{
    cout << cu_steps << "/" << ov_steps << ": Calculate Integral for H elements.("<<G1.n_fixed<<"*"<<G1.n_fixed<<")\n";
    fout << cu_steps << "/" << ov_steps << ": Calculate Integral for H elements.("<<G1.n_fixed<<"*"<<G1.n_fixed<<")\n";
    cu_steps+=1;
    cout << cu_steps << "/" << ov_steps << ": Construct H.";
    fout << cu_steps << "/" << ov_steps << ": Construct H.";
    H = new double[G1.n_fixed * G1.n_fixed]();
    for (int i = 0; i < G1.n_fixed; ++i)
    {
        for (int j = i; j < G1.n_fixed; ++j) // symmetry
        {
            H[i * G1.n_fixed + j] = integral_as_sum(G1, i, j);
        }
    }
    
    for (int i = 0; i < G1.n_fixed; ++i)
    {
        for (int j = i; j < G1.n_fixed; ++j) // symmetry
        {
            cout<<"result of element for "<<i<<"  "<<j<<" ";
            cout<<H[i * G1.n_fixed + j]<<endl;
        }
    }
    
    cout << "   -->Done." << endl;
    fout << "   -->Done." << endl;
    cu_steps += 1;
}
