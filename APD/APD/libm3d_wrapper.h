#ifndef LIBM3D_WRAPPER_H
#define LIBM3D_WRAPPER_H


#include "data_vector.h"
#include "mesh_io.h"

#include <libm3d.h>

namespace APDL
{
	inline std::vector<DataVector> computePointMKDiff(libm3d::model& P, libm3d::model& Q, double R[3][3], double d, double shift)
	{
		mathtool::Matrix3x3 m(
			R[0][0], R[0][1], R[0][2], 
			R[1][0], R[1][1], R[1][2], 
			R[2][0], R[2][1], R[2][2]);

		Q.rotate(m);

		libm3d::LibM3dParam param;
		param.PmaxD = param.QmaxD = d;
		param.threadsize = 4;

		double current_rot[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

		P.sample(param.PmaxD, param.escaleP);
		Q.sample(param.QmaxD, param.escaleQ);

		Q.negate();

		P.build_collision_detection(true);
		Q.build_collision_detection(false);

		std::vector<std::vector<float> > points;

		ComputePointMKSum(P, Q, current_rot, param, points);

		std::vector<DataVector> results;
		DataVector p(3), n(3);
		for(std::size_t i = 0; i < points.size(); ++i)
		{
			p[0] = points[i][0]; p[1] = points[i][1]; p[2] = points[i][2];
			n[0] = points[i][3]; n[1] = points[i][4]; n[2] = points[i][5];

			results.push_back(p + shift * n);
			results.push_back(p - shift * n);
		}

		return results;

	}
}

#endif