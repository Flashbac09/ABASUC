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
    system(ss.str().c_str());
    fout.open(folder+"/OUTPUT.txt");
}

// Global //