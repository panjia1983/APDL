#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/math_utility.h>

namespace APDL
{	
	void test()
	{
		{
			Quaternion q(1, 3, 4, 5);
			for(int i = 0; i < 4; ++i)
				std::cout << q[i] << std::endl;
			double R[3][3];
			Quat2Rot(R, q);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}
			Rot2Quat(q, R);
			for(int i = 0; i < 4; ++i)
				std::cout << q[i] << std::endl;
			Quat2Rot(R, q);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}

			return;

		}
		{
			double a = 1, b = 2, c = 3;
			double R[3][3];
			Euler2Rot(R, a, b, c);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}
			Rot2Euler(a, b, c, R);
			std::cout << a << " " << b << " " << c << std::endl;

			Euler2Rot(R, a, b, c);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}
		}

		{
			Quaternion q(1, 2, 3, 4);
			double a, b, c;
			Quat2Euler(a, b, c, q);
			double R[3][3];
			Quat2Rot(R, q);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}

			Euler2Rot(R, a, b, c);
			for(int i = 0; i < 3; ++i)
			{
				for(int j = 0; j < 3; ++j)
					std::cout << R[i][j] << " ";
				std::cout << std::endl;
			}
		}

		{
			double R[3][3];
			R[0][0] = 0.0;
			R[0][1] = -1.0;
			R[0][2] = 0.0;
			R[1][0] = 0.984808;
			R[1][1] = 0.0;
			R[1][2] = 0.173648;
			R[2][0] = -0.173648;
			R[2][1] = 0.0;
			R[2][2] = 0.984808;

			double a, b, c;
			Rot2Euler(a, b, c, R);
			Quaternion q;
			Euler2Quat(q, a, b, c);

			Euler2Rot(R, a, b, c);
			for(int j = 0; j < 3; ++j)
			{
				for(int k = 0; k < 3; ++k)
					std::cout << R[j][k] << " ";
				std::cout << std::endl;
			}

			Quat2Rot(R, q);
			for(int j = 0; j < 3; ++j)
			{
				for(int k = 0; k < 3; ++k)
					std::cout << R[j][k] << " ";
				std::cout << std::endl;
			}
		}
	}
}

void main()
{

}

