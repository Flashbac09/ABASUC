#include "eigen.h"

inline double eigen::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
}

inline double eigen::drdf(double x, double h) // rdf from distancce
{
    int i;
    double xi, xip;
    i = (int)(x / h);
    xi = i * h;
    xip = xi + h;
    double res = 0;
    res = (1 + 2 * (x - xi) / h) * ((x - xip) / h) * ((x - xip) / h) * y[i] +
          (1 - 2 * (x - xip) / h) * ((x - xi) / h) * ((x - xi) / h) * y[i + 1] +
          (x - xi) * ((x - xip) / h) * ((x - xip) / h) * m[i] +
          (x - xip) * ((x - xi) / h) * ((x - xi) / h) * m[i + 1];
    return res;
}
void eigen::spline_interpolation(grid &G1)
{
    cout << cu_steps << "/" << ov_steps << ": Cubic spline interpolation for RDF.";
    fout << cu_steps << "/" << ov_steps << ": Cubic spline interpolation for RDF.";

    mesh = G1.mesh;
    double *al = new double[mesh]();
    double *be = new double[mesh]();
    double *A = new double[mesh]();
    double *B = new double[mesh]();
    int n = mesh - 1;
    double h = G1.dr;
    y = new double[mesh]();
    m = new double[mesh]();
    for (int i = 0; i < mesh; ++i)
    {
        y[i] = G1.rdf[i];
    }

    // spline, hermite //
    al[0] = 1;
    al[n] = 0;
    be[0] = 3 * (y[1] - y[0]) / h;
    be[n] = 3 * (y[n] - y[n - 1]) / h;
    A[0] = -0.5 * al[0];
    B[0] = 0.5 * be[0];
    for (int i = 1; i <= n - 1; i++)
    {
        al[i] = 0.5;
        be[i] = 3 * ((1 - al[i]) * (y[i] - y[i - 1]) / h + al[i] * (y[i + 1] - y[i]) / h);
    }
    for (int i = 1; i <= n - 1; i++)
    {
        A[i] = (-al[i]) / (2 + (1 - al[i]) * A[i - 1]);
        B[i] = (be[i] - (1 - al[i]) * B[i - 1]) / (2 + (1 - al[i]) * A[i - 1]);
    }
    B[n] = (be[n] - (1 - al[n]) * B[n - 1]) / (2 + (1 - al[n]) * A[n - 1]);
    m[n] = B[n];

    for (int i = n - 1; i >= 0; i--)
    {
        m[i] = A[i] * m[i + 1] + B[i];
    }
    // spline, hermite //
    cout << "   -->Done." << endl;
    fout << "   -->Done." << endl;
    cu_steps += 1;
    delete[] al;
    delete[] be;
    delete[] A;
    delete[] B;
    delete[] y;
    delete[] m;
}
inline double eigen::integral_as_sum(grid &G1, int i, int j) // inline/parallel? //
{
    double dis = sqrt(pow((G1.x[i] - G1.x[j]), 2) + pow((G1.y[i] - G1.y[j]), 2) + pow((G1.z[i] - G1.z[j]), 2));
    if (dis > 2 * G1.rdf_cutoff)
        return 0;
    // 1. determine overlap cuboid x,y,z
    double Superior_x, Inferior_x, Superior_y, Inferior_y, Superior_z, Inferior_z;
    Superior_x = max(G1.x[i], G1.x[j]);
    Inferior_x = min(G1.x[i], G1.x[j]);
    Superior_y = max(G1.y[i], G1.y[j]);
    Inferior_y = min(G1.y[i], G1.y[j]);
    Superior_z = max(G1.z[i], G1.z[j]);
    Inferior_z = min(G1.z[i], G1.z[j]);
    double forward_x, forward_y, forward_z, back_x, back_y, back_z;
    forward_x = min(Inferior_x + G1.rdf_cutoff, (double)G1.lx); // do not exceed the boundary
    back_x = max(Superior_x - G1.rdf_cutoff, 0.0);
    forward_y = min(Inferior_y + G1.rdf_cutoff, (double)G1.ly);
    back_y = max(Superior_y - G1.rdf_cutoff, 0.0);
    forward_z = min(Inferior_z + G1.rdf_cutoff, (double)G1.lz);
    back_z = max(Superior_z - G1.rdf_cutoff, 0.0);
    // 2. go through grid and sum
    int x_start, y_start, z_start, x_end, y_end, z_end;
    // here notice if points's positions are special(marginal part)
    x_start = (int)(back_x / G1.d1);
    y_start = (int)(back_y / G1.d2);
    z_start = (int)(back_z / G1.d3);
    x_end = min(((int)(forward_x / G1.d1) + 1), G1.nx - 1); // do not exceed the n-grid.
    y_end = min(((int)(forward_y / G1.d2) + 1), G1.ny - 1);
    z_end = min(((int)(forward_z / G1.d3) + 1), G1.nz - 1);
    double res = 0;
    for (int m = x_start; m <= x_end; ++m)
    {
        for (int n = y_start; n <= y_end; ++n)
        {
            for (int k = z_start; k <= z_end; ++k)
            {
                // every point : x=m*d1,y=n*d2,z=k*d3, check the both sphere within this cuboid.
                if (i == j)
                {
                    if (sqrt(pow((m * G1.d1 - G1.x[i]), 2) + pow((n * G1.d2 - G1.y[i]), 2) + pow((k * G1.d3 - G1.z[i]), 2)) < G1.rdf_cutoff)
                        res += pow(drdf(distance(m * G1.d1, n * G1.d2, k * G1.d3, G1.x[i], G1.y[i], G1.z[i]), G1.dr), 2) * G1.gd[m * G1.ny * G1.nz + n * G1.nz + k] * G1.d1 * G1.d2 * G1.d3;
                }
                else if (sqrt(pow((m * G1.d1 - G1.x[i]), 2) + pow((n * G1.d2 - G1.y[i]), 2) + pow((k * G1.d3 - G1.z[i]), 2)) < G1.rdf_cutoff && sqrt(pow((m * G1.d1 - G1.x[j]), 2) + pow((n * G1.d2 - G1.y[j]), 2) + pow((k * G1.d3 - G1.z[j]), 2)) < G1.rdf_cutoff)
                    res += drdf(distance(m * G1.d1, n * G1.d2, k * G1.d3, G1.x[i], G1.y[i], G1.z[i]), G1.dr) * drdf(distance(m * G1.d1, n * G1.d2, k * G1.d3, G1.x[j], G1.y[j], G1.z[j]), G1.dr) * G1.gd[m * G1.ny * G1.nz + n * G1.nz + k] * G1.d1 * G1.d2 * G1.d3;
            }
        }
    }
    return res;
}
/*
void mpiHsolver(double *A)
{
    MPI_Init(NULL, NULL);

   
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // my scalapack context
    int num_procs_per_row = sqrt(size);
    int num_procs_per_col = size / num_procs_per_row;
    int context, my_row, my_col;
    Cblacs_pinfo(&rank, &size);
    Cblacs_get(-1, 0, &context);
    Cblacs_gridinit(&context, "Row-major", num_procs_per_row, num_procs_per_col);
    Cblacs_gridinfo(context, &num_procs_per_row, &num_procs_per_col, &my_row, &my_col);

    int block_size = N / num_procs_per_row;
    double *local_A = (double *)malloc(block_size * block_size * sizeof(double));
    int descA[9];
    int info;
    int na = N, nb = N, lda = block_size, ldb = block_size, rsrc = 0, csrc = 0;
    int lld_A = lda + 2;
    descinit_(descA, &na, &nb, &block_size, &block_size, &rsrc, &csrc, &context, &lld_A, &info);
    pdistribute_(A, descA, local_A, descA, &my_row, &my_col, &context);

    double *W = (double *)malloc(block_size * sizeof(double));
    double *Z = (double *)malloc(block_size * block_size * sizeof(double));
    int lwork = -1, liwork = -1;
    double wkopt, work[1];
    int iwkopt, iwork[1];
    pdsyev_("V", "U", &block_size, local_A, &rsrc, &csrc, descA, W, Z, &lld_A, &wkopt, &lwork, &iwkopt, &liwork, &info);
    lwork = (int)wkopt;
    liwork = iwkopt;
    double *work2 = (double *)malloc(lwork * sizeof(double));
    int *iwork2 = (int *)malloc(liwork * sizeof(int));
    pdsyev_("V", "U", &block_size, local_A, &rsrc, &csrc, descA, W, Z, &lld_A, work2, &lwork, iwork2, &liwork, &info);

    if (rank == 0)
    {
        printf("Eigenvalues:\n");
        for (int i = 0; i < N; ++i)
        {
            printf("%.6f\n", W[i]);
        }
        printf("\nEigenVectors:\n");
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                printf("%.6f ", Z[i * N + j]);
            }
            printf("\n");
        }
    }

    if (rank == 0)
    {
        free(A);
    }
    free(local_A);
    free(W);
    free(Z);
    MPI_Finalize();
}
*/
void eigen::Hsolver()
{
    // link
    // LAPACKE
    cout << cu_steps << "/" << ov_steps << ": Solve eigenpairs of H.";
    fout << cu_steps << "/" << ov_steps << ": Solve eigenpairs of H.";
    t = new double[n_fixed];
    H_copy=new double [n_fixed*n_fixed]();
    for (int i = 0; i < n_fixed; ++i)
    {
        for (int j = 0; j < n_fixed; ++j)
        {
            H_copy[i * n_fixed + j]= H[i * n_fixed + j] ;
        }
    }
    int calc = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n_fixed, H, n_fixed, t);
    // major/jobz/uplo/row/data/lda/eigenvalue/
    // eigenvector in every column of H. H[i*n_fixed+j]
    //pdsyev("")
    if (calc != 0)
    {
        cout << "     -->Error.\nLAPACKE didn't run successfully.\n";
        fout << "     -->Error.\nLAPACKE didn't run successfully.\n";
        exit(1);
    }
    cu_steps += 1;
    cout << "   -->Done." << endl;
    fout << "   -->Done." << endl;
}
void eigen::output_finalize(grid &G1)
{
    cout << cu_steps << "/" << ov_steps << ": Output and finalize.\n\n";
    fout << cu_steps << "/" << ov_steps << ": Output and finalize.\n\n";
    cu_steps += 1;
    cout << "-----Check complete output in folder \"output\"-----" << endl;
    fout << "-----Check complete output in folder \"output\"-----" << endl;
    fout.close();
    fout.open("output/eigenvalues.log");
    for (int i = 0; i < n_fixed; i++)
    {
        fout << t[i] << endl;
    }
    fout.close();
    fout.open("output/eigenvectors.log");
    for (int i = 0; i < n_fixed; i++)
    {
        for (int j = 0; j < n_fixed; j++)
        {
            fout << H[i * n_fixed + j] << " ";
        }
        fout << endl;
    }
    fout.close();
    fout.open("output/H.log");
    fout<<"-----Only Upper part of H is valid, for H is a symmetric matrix.-----"<<endl;
    for (int i = 0; i < n_fixed; i++)
    {
        for (int j = 0; j < n_fixed; j++)
        {
            fout << H_copy[i * n_fixed + j] << " ";
        }
        fout << endl;
    }
    delete[] t;
    delete[] H;
    delete[] G1.rdf;
    delete[] G1.gd;
    delete[] H_copy;
}
void eigen::Hconstructor(grid &G1)
{
    cout << cu_steps << "/" << ov_steps << ": Calculate Integral for H elements(" << G1.n_fixed << "*" << G1.n_fixed << ").\n";
    fout << cu_steps << "/" << ov_steps << ": Calculate Integral for H elements(" << G1.n_fixed << "*" << G1.n_fixed << ").\n";
    cu_steps += 1;
    cout << cu_steps << "/" << ov_steps << ": Construct H.";
    fout << cu_steps << "/" << ov_steps << ": Construct H.";
    H = new double[G1.n_fixed * G1.n_fixed]();
    n_fixed = G1.n_fixed;
    for (int i = 0; i < G1.n_fixed; ++i)
    {
        for (int j = i; j < G1.n_fixed; ++j) // symmetry upper
        {
            H[i * G1.n_fixed + j] = integral_as_sum(G1, i, j);
        }
    }
    cout << "   -->Done." << endl;
    fout << "   -->Done." << endl;
    cu_steps += 1;
}

// draft:

/*
    n_fixed = 5;
    double t[5] = {0};
    double a[25] = {
        1.96,
        -6.49,
        -0.47,
        -7.20,
        -0.65,
        -6.49,
        3.80,
        -6.39,
        1.50,
        -6.34,
        -0.47,
        -6.39,
        4.17,
        -1.51,
        2.67,
        -7.20,
        1.50,
        -1.51,
        5.70,
        1.80,
        -0.65,
        -6.34,
        2.67,
        1.80,
        -7.10,
    };
    */
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

/*

inline double eigen::integral_as_sum(grid &G1, int i, int j) // inline/parallel? //
{
    double dis = sqrt(pow((G1.x[i] - G1.x[j]), 2) + pow((G1.y[i] - G1.y[j]), 2) + pow((G1.z[i] - G1.z[j]), 2));
    if (dis > 2 * G1.rdf_cutoff)
        return 0;
    cout << "dis: " << dis << endl;
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
    cout << "Superior/Inferior x,y,z for " << i << "  " << j << endl;
    cout << Superior_x << "  " << Inferior_x << endl;
    cout << Superior_y << "  " << Inferior_y << endl;
    cout << Superior_z << "  " << Inferior_z << endl;

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
    cout << "Back/Forward x,y,z for " << i << "  " << j << endl;
    cout << back_x << "  " << forward_x << endl;
    cout << back_y << "  " << forward_y << endl;
    cout << back_z << "  " << forward_z << endl;
    //
    if (Overlap_x <= 0 || Overlap_y <= 0 || Overlap_z <= 0)
        return 0;//for stability
    //
    // 3. go through grid and sum
    int x_start, y_start, z_start, x_end, y_end, z_end;
    // here check if points's positions are special(marginal part)
    x_start = (int)(back_x / G1.d1);
    y_start = (int)(back_y / G1.d2);
    z_start = (int)(back_z / G1.d3);
    x_end = min(((int)(forward_x / G1.d1) + 1), G1.nx - 1); // do not exceed the grid.
    y_end = min(((int)(forward_y / G1.d2) + 1), G1.ny - 1);
    z_end = min(((int)(forward_z / G1.d3) + 1), G1.nz - 1);
    cout << "start/end of x,y,z for " << i << "  " << j << endl;
    cout << x_start << "  " << x_end << endl;
    cout << y_start << "  " << y_end << endl;
    cout << z_start << "  " << z_end << endl;
    double res = 0;
    for (int m = x_start; m <= x_end; ++m)
    {
        for (int n = y_start; n <= y_end; ++n)
        {
            for (int k = z_start; k <= z_end; ++k)
            {
                // cout<<m<<"  "<<n<<"  "<<k<<"  "<<G1.gd[m*G1.ny*G1.nz+n*G1.nz+k]<<endl;
                // cout<<res<<endl;
                //every point : x=m*d1,y=n*d2,z=k*d3
                if (sqrt(pow((m*G1.d1-G1.x[i]),2)+pow((n*G1.d2-G1.y[i]),2)+pow((k*G1.d3-G1.z[i]),2))<G1.rdf_cutoff&&sqrt(pow((m*G1.d1-G1.x[j]),2)+pow((n*G1.d2-G1.y[j]),2)+pow((k*G1.d3-G1.z[j]),2))<G1.rdf_cutoff)
                    res += G1.gd[m * G1.ny * G1.nz + n * G1.nz + k] * rdf(G1, i) * rdf(G1, j) * G1.d1 * G1.d2 * G1.d3;
            }
        }
    }
    return res;
}

*/
