#include "libm3d.h"

int main()
{
	std::string file_name1 = "../data/grate1.obj";
	std::string file_name2 = "../data/grate2.obj";

	LibM3dParam param;
	param.PmaxD = param.QmaxD = 10;
	param.threadsize = 4;

	model P, Q;
	if(!P.build(file_name1, false)) return 1;
	if(!Q.build(file_name2, true)) return 1;

	double current_rot[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

	Point3d current_radian;

	if(param.perturb != 0)
	{
		current_radian.set(param.perturb*drand48(),param.perturb*drand48(),param.perturb*drand48());
		Matrix3x3 tmpM;
		computeRotationMatrix(current_radian,tmpM);
		tmpM.get(current_rot);
		P.rotate(tmpM);
		cout<<"- Perturb P using: ("<<current_radian<<")"<<endl;
	}

	P.sample(param.PmaxD, param.escaleP);
	Q.sample(param.QmaxD, param.escaleQ);

	Q.negate();

	P.build_collision_detection(true);
	Q.build_collision_detection(false);

	MKPTS allBds;
	ComputePointMKSum(P, Q, current_rot, param, allBds);

	std::cout << allBds.size() << std::endl;

	return 0;
}