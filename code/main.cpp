#include "ObsFile.h"
#include "NavFile.h"
#include "SPP.h"


int  main()
{
	
    NavFile navfile = NavFile("brdc0300.17n");
    ObsFile obsfile = ObsFile("bjnm0300.17o");



    Adjust_Scheme scheme = { 1,15.0 ,-1,0,1,0.5};

    SPP::propos(navfile, obsfile,scheme);
    

    return 0;

}
