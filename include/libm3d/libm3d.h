#ifndef LIBM3D_H
#define LIBM3D_H

#include "model.h"
#include "minkowski.h"

namespace libm3d 
{

struct LibM3dParam
{
	bool bOctfilter;
	bool bCDfilter;
	bool bNormalfilter;

	float PmaxD;
	float QmaxD;

	int escaleP;
	int escaleQ;

	float perturb;

	int threadsize;

	LibM3dParam()
	{
		bOctfilter = false;
		bCDfilter = true;
		bNormalfilter = true;

		PmaxD = 1e10;
		QmaxD = 1e10;

		escaleP = 1;
		escaleQ = 1;

		perturb = 0;

		threadsize = 1;
	}
};



void ComputePointMKSum(model& P, model& Q, 
					   double R[3][3],
					   const LibM3dParam& param, std::vector<std::vector<float> >& points);

}

#endif