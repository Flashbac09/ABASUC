#include "input.h"

void input::str_tolower(char t[], char nt[])
{
    char a;
    int l = strlen(t);
    for (int i = 0; i < l; i++)
    {
        a = t[i];
        nt[i] = tolower(a);
    }
    nt[l] = '\0';
}

template<typename T>
void input::read_values(ifstream &fin, T &var)
{
    fin >> var;
    fin.ignore(100, '\n');
}
//Note: for g++ O2/O3, we need to declare all possbile functions from template function.
//In this case, we only need to give instances used in grid.cpp, not in input.cpp itself.
//If not,there will be 'undefined reference' error.
template void input::read_values<int>(ifstream &fin,int &var);
template void input::read_values<double>(ifstream &fin,double &var);
void input::essential_para_init()
{
    lx = 0;
    ly = 0;
    lz = 0;
    diago_lib = "";
    points_path = "";
    venergy_path = "";
    distribution_path = "";
}
bool input::essential_para_check()
{
    if (lx <= 0 || ly <= 0 || lz <= 0)
    {
        return 0;
    }
    else if (strcmp(diago_lib.c_str(), "") == 0)
        return 0;
    else if (strcmp(points_path.c_str(), "") == 0)
        return 0;
    else if (strcmp(venergy_path.c_str(), "") == 0)
        return 0;
    else if (strcmp(distribution_path.c_str(), "") == 0)
        return 0;
    else
        return 1;
}
void input::input_initialize(const string &file_input, ofstream &fout)
{
    cout << cu_steps << "/" << ov_steps << ": Read in parameters & generate cuboid model.";
    fout << cu_steps << "/" << ov_steps << ": Read in parameters & generate cuboid model.";
    ifstream fin(file_input);
    if (!fin.is_open())
    {
        cout << "   -->Error.\ncan't open input file: " << file_input << endl;
        fout << "   -->Error.\ncan't open input file: " << file_input << endl;
        exit(1);
    }
    fin.clear();
    fin.seekg(0);
    fin.rdstate();
    char temp[100], ntemp[100];
    essential_para_init();
    while (fin.good())
    {
        fin >> ntemp;
        // str_tolower(temp,ntemp);// use concise input keys
        if (strcmp(ntemp, "isHexahedral") == 0)
            read_values(fin, isHexahedral);
        else if (strcmp(ntemp, "lx") == 0)
            read_values(fin, lx);
        else if (strcmp(ntemp, "ly") == 0)
            read_values(fin, ly);
        else if (strcmp(ntemp, "lz") == 0)
            read_values(fin, lz);
        else if (strcmp(ntemp, "thetaxy") == 0)
            read_values(fin, thetaxy);
        else if (strcmp(ntemp, "thetayz") == 0)
            read_values(fin, thetayz);
        else if (strcmp(ntemp, "thetaxz") == 0)
            read_values(fin, thetaxz);
        else if (strcmp(ntemp, "support_SH") == 0)
            read_values(fin, support_SH);
        else if (strcmp(ntemp, "diago_lib") == 0)
            read_values(fin, diago_lib);
        else if (strcmp(ntemp, "support_Periodic_Boundary") == 0)
            read_values(fin, support_Periodic_Boundary);
        else if (strcmp(ntemp, "multi_parallel_strategies") == 0)
            read_values(fin, multi_parallel_strategies);
        else if (strcmp(ntemp, "points_path") == 0)
            read_values(fin, points_path);
        else if (strcmp(ntemp, "venergy_path") == 0||strcmp(ntemp,"v_path")==0)
            read_values(fin, venergy_path);
        else if (strcmp(ntemp, "distribution_path") == 0)
            read_values(fin, distribution_path);
        else
        {
            cout << "\nWarning: the parameter " << ntemp << " is not preset." << endl;
            fout << "\nWarning: the parameter " << ntemp << " is not preset." << endl;
            fin.ignore(100, '\n');
        }
        fin.rdstate();
        if (fin.eof() != 0)
        {
            break;
        }
        else if (fin.fail() != 0)
        {
            cout << "   -->Error, in reading parameter " << ntemp << endl;
            fout << "   -->Error, in reading parameter " << ntemp << endl;
            exit(1);
        }
        else if (fin.good() == 0)
        {
            break;
        }
    }
    fin.close();
    if (!essential_para_check())
    {
        cout << "   -->Error.\nsome parameters don't have default value and you have to set it." << endl;
        fout << "   -->Error.\nsome parameters don't have default value and you have to set it." << endl;
        exit(1);
    }
    cout << "   -->Success." << endl;
    fout << "   -->Success." << endl;
    cu_steps += 1;
}
