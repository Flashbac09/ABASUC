#include "grid.h"
bool grid::good_grid_condition()
{
    /*
    if (((double)lx) / nx != ((double)ly) / ny)
        return 0;
    else if (((double)lx) / nx != ((double)lz) / nz)
        return 0;
        */
       //else 
    {
        d1 = (double)lx / nx;
        d2=(double)ly / ny;
        d3=(double)lz / nz;
        return 1;
    }
}
bool grid::good_fixed_points(int n_fixed)
{
    for(int i=0;i<n_fixed;++i)
    {
        if(x[i]<0||x[i]>lx||y[i]<0||y[i]>ly||z[i]<0||z[i]>lz)
        return 0;
    }
    return 1;
}
void grid::read_fixed_points(const string &file_points, ofstream &fout)
{
    cout << cu_steps << "/" << ov_steps << ": Read in fixed points.";
    fout << cu_steps << "/" << ov_steps << ": Read in fixed points.";
    ifstream fin(file_points);
    if (!fin.is_open())
    {
        cout << "   -->Error.\ncan't open input file: " << file_points << endl;
        fout << "   -->Error.\ncan't open input file: " << file_points << endl;
        exit(1);
    }
    fin.clear();
    fin.seekg(0);
    fin.rdstate();
    char temp[100], bit_f, bit_b, bit_d1, bit_d2, bit_n;
    n_fixed = 0;
    while (fin.good())
    {
        fin >> bit_f >> x[n_fixed] >> bit_d1 >> y[n_fixed] >> bit_d2 >> z[n_fixed] >> bit_b;
        fin.ignore(100, '\n'); // crucial '\n' //
        if (bit_f != '(' || bit_b != ')' || bit_d1 != ',' || bit_d2 != ',')
        {
            cout << "   -->Error.\nFormat of points is not compatible." << endl;
            fout << "   -->Error.\nFormat of points is not compatible." << endl;
            exit(1);
        }
        else
        {
            n_fixed += 1;
        }
        fin.rdstate();
        if (fin.eof() != 0)
        {
            break;
        }
        else if (fin.fail() != 0)
        {
            cout << "   -->Error, in reading points " << endl;
            fout << "   -->Error, in reading points " << endl;
            exit(1);
        }
        else if (fin.good() == 0)
        {
            break;
        }
    }
    fin.close();
    if(!good_fixed_points(n_fixed))
    {
        cout<<"   -->Error, points have conflicts with cuboid model."<<endl;
        fout<<"   -->Error, points have conflicts with cuboid model."<<endl;
        exit(1);
    }
    /*
    for(int i=0;i<n_fixed;++i)
    {
        cout<<x[i]<<"  "<<y[i]<<"  "<<z[i]<<endl;
    }
    */
    cout << "   -->Success." << endl;
    fout << "   -->Success." << endl;
    cu_steps += 1;
}
void grid::read_venergy_and_construct(input &I1, ofstream &fout)
{
    // 1:grid
    cout << cu_steps << "/" << ov_steps << ": Generate uniform grid.";
    fout << cu_steps << "/" << ov_steps << ": Generate uniform grid.";
    ifstream fin(I1.venergy_path);
    if (!fin.is_open())
    {
        cout << "   -->Error, can't open input file: " << I1.venergy_path << endl;
        fout << "   -->Error, can't open input file: " << I1.venergy_path << endl;
        exit(0);
    }
    fin.clear();
    fin.seekg(0);
    fin.rdstate();
    char temp[100], ntemp[100];
    lx = I1.lx;
    ly = I1.ly;
    lz = I1.lz;
    nx = 0;
    ny = 0;
    nz = 0;
    while (fin.good() && !fin.eof())
    {
        fin >> ntemp;
        if (strcmp(ntemp, "nx") == 0)
            input::read_values(fin, nx);
        else if (strcmp(ntemp, "ny") == 0)
            input::read_values(fin, ny);
        else if (strcmp(ntemp, "nz") == 0)
            input::read_values(fin, nz);
        else if (strcmp(ntemp, "V:") == 0)
        {
            fin.ignore(100, '\n');
            break;
        }
    }
    if (!(nx > 0 && ny > 0 && nz > 0))
    {
        cout << "   -->Error, in reading grid parameter nx/ny/nz." << endl;
        fout << "   -->Error, in reading grid parameter nx/ny/nz." << endl;
        exit(1);
    }
    if (!good_grid_condition())
    {
        cout << "   -->Error, non-uniform grid." << endl;
        fout << "   -->Error, non-uniform grid." << endl;
        exit(1);
    }
    // vector? //
    gd = new double[(nx) * (ny) * (nz) + 1]();
    cout << "   -->Success." << endl;
    fout << "   -->Success." << endl;
    cu_steps += 1;
    // 2:V energy
    cout << cu_steps << "/" << ov_steps << ": Read in V energy.";
    fout << cu_steps << "/" << ov_steps << ": Read in V energy.";
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                fin >> gd[i * ny * nz + j * nz + k];
            }
        }
    }
    fin.rdstate();
    if (fin.fail() != 0)
    {
        cout << "   -->Error, in reading V energy. " << endl;
        fout << "   -->Error, in reading V energy. " << endl;
        exit(1);
    }
    fin.close();
    cout << "   -->Success." << endl;
    fout << "   -->Success." << endl;
    cu_steps += 1;
}
void grid::rdf_para_init()
{
    mesh = 0;
    dr = 0;
    rdf_cutoff = 0;
    l = 0;
}
bool grid::rdf_para_check()
{
    if (mesh == 0 || dr == 0 || rdf_cutoff == 0 || l == 0)
        return 0;
    else
        return 1;
}
void grid::read_distribution(const string &file_distribution, ofstream &fout)
{
    cout << cu_steps << "/" << ov_steps << ": Read in Radius Distribution Function(RDF).";
    fout << cu_steps << "/" << ov_steps << ": Read in Radius Distribution Function(RDF).";
    ifstream fin(file_distribution);
    if (!fin.is_open())
    {
        cout << "   -->Error.\ncan't open input file: " << file_distribution << endl;
        fout << "   -->Error.\ncan't open input file: " << file_distribution << endl;
        exit(1);
    }
    fin.clear();
    fin.seekg(0);
    fin.rdstate();
    char ntemp[100], bit_d;
    rdf_para_init();
    while (fin.good())
    {
        fin >> ntemp;
        if (strcmp(ntemp, "cutoff") == 0)
            input::read_values(fin, rdf_cutoff);
        else if (strcmp(ntemp, "mesh") == 0)
            input::read_values(fin, mesh);
        else if (strcmp(ntemp, "dr") == 0)
            input::read_values(fin, dr);
        else if (strcmp(ntemp, "l") == 0)
            input::read_values(fin, l);
        else if (strcmp(ntemp, "f:") == 0)
        {
            fin.ignore(100, '\n');
            if (mesh <= 0)
            {
                cout << "   -->Error, mesh is not set." << endl;
                fout << "   -->Error, mesh is not set." << endl;
            }
            rdf = new double[mesh]();
            for (int i = 0; i < mesh; ++i)
            {
                fin >> rdf[i] >> bit_d;
            }
            ntemp[0] = '\0'; // reset
        }
        fin.rdstate();
        if (fin.eof() != 0)
        {
            break;
        }
        else if (fin.fail() != 0)
        {
            cout << "   -->Error, in reading RDF " << endl;
            fout << "   -->Error, in reading RDF " << endl;
            exit(1);
        }
        else if (fin.good() == 0)
        {
            break;
        }
    }
    fin.close();
    if (!rdf_para_check())
    {
        cout << "   ->Error, lack necessary RDF parameter." << endl;
        fout << "   ->Error, lack necessary RDF parameter." << endl;
        exit(1);
    }
    cout << "   -->Success." << endl;
    fout << "   -->Success." << endl;
    cu_steps += 1;
}