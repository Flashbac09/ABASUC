#ifndef _INPUT_H_
#define _INPUT_H_
#include "global.h"
class input
{
public:
    // 1:essential
    int lx;
    int ly;
    int lz;
    string diago_lib;
    string points_path;
    string venergy_path;
    string distribution_path;
    // 2:alternative
    int support_SH = 0;
    int support_Periodic_Boundary = 0;
    int multi_parallel_strategies = 0;
    int isHexahedral = 0;
    double thetaxy = 0;
    double thetayz = 0;
    double thetaxz = 0;
    // 3:function
    static void str_tolower(char t[], char nt[]);
    void input_initialize(const string &file_input, ofstream &fout);
    void input_distribution(const string &file_distribution, ofstream &fout);
    void input_points(const string &file_points, ofstream &fout);
    void input_venergy(const string &file_venergy, ofstream &fout);
    void generate_uniform_grid(const string &file_venergy, ofstream &fout);
    void essential_para_init();
    bool essential_para_check();
    //
    template <typename T>
    static void read_values(ifstream &fin, T &var);
};

#endif