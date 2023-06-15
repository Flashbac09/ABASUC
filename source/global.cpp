#include "global.h"

// Global //

int cu_steps;
int ov_steps;
ofstream fout;
void preset_overall()
{
    //1:
    cu_steps = 1;
    ov_steps = 10;
    //2:
    string folder = "output";
    stringstream ss;
    ss << " test -d " << folder << " || mkdir " << folder;
    if(system(ss.str().c_str())) //return 0 if the command runs successfullly.
    {
        cout<<"   -->Error.Can't create output folder.\n";
        fout<<"   -->Error.Can't create output folder.\n";
        exit(1);
    }
    fout.open(folder+"/OUTPUT.txt");
}

// Global //