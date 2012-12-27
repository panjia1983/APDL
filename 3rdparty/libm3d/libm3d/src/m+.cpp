#include "m+.h"

namespace libm3d 
{

bool parseArg(int argc, char ** argv)
{
	if(argc!=4) return false;
	D=(float)atof(argv[1]);
	P.build(argv[2]);
	Q.build(argv[3]);
}

void printUsage()
{
	cerr<<"usage: m+ res *.obj *.obj"<<endl;
}

int main(int argc, char ** argv)
{
	if(parseArg(argc,argv)){
		printUsage();
		return 1;
	}

	{
		float box[6];
		bbox_point_minkowski_sum(P,Q,box);
		G.build(box,D);
	}
	
	std::list<mathtool::Point3d> M; //minkowski sum boundary

	{
		G.registering(P,Q);
		G.first_wave();
		G.second_wave();
		G.getBoundaryPts(M);
		G.destroy();
	}
}

}