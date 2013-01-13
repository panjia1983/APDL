#ifndef MATH_UTILITY_H
#define MATH_UTILITY_H

#include <boost/math/constants/constants.hpp>
#include "data_vector.h"

namespace APDL
{

	class Quaternion
	{
		double data[4];
	public:
		Quaternion(double a, double b, double c, double d)
		{
			double inv_len = 1.0 / sqrt(a * a + b * b + c * c + d * d);
			data[0] = a * inv_len;
			data[1] = b * inv_len;
			data[2] = c * inv_len;
			data[3] = d * inv_len;
		}
		
		Quaternion()
		{
			data[0] = 1;
			data[1] = 0;
			data[2] = 0;
			data[3] = 0;
		}
		
		Quaternion& inverse()
		{
			data[1] = -data[1];
			data[2] = -data[2];
			data[3] = -data[3];
			return *this;
		}
		
		Quaternion operator * (const Quaternion& other) const
		{
			return Quaternion(data[0] * other.data[0] - data[1] * other.data[1] - data[2] * other.data[2] - data[3] * other.data[3],
			                  data[0] * other.data[1] + data[1] * other.data[0] + data[2] * other.data[3] - data[3] * other.data[2],
			                  data[0] * other.data[2] - data[1] * other.data[3] + data[2] * other.data[0] + data[3] * other.data[1],
			                  data[0] * other.data[3] + data[1] * other.data[2] - data[2] * other.data[1] + data[3] * other.data[0]);
		}
		
		double operator [](std::size_t i) const
		{
			return data[i];
		}
		double& operator [](std::size_t i)
		{
			return data[i];
		}
	};

	inline Quaternion inverse(const Quaternion& q)
	{
		return Quaternion(q[0], -q[1], -q[2], -q[3]);
	}
	
	void Quat2Euler(double& a, double& b, double& c, const Quaternion& quat);
	
	void Quat2Rot(double R[3][3], const Quaternion& quat);
	
	void Euler2Rot(double R[3][3], double a, double b, double c);
	
	void Euler2Quat(Quaternion& quat, double a, double b, double c);

	void Rot2Quat(Quaternion& quat, double R[3][3]);

	void Rot2Euler(double& a, double& b, double& c, double R[3][3]);
	
	double angleTruncate(double angle);
	
	DataVector InterpConfig2D(const DataVector& q1, const DataVector& q3, double t);
	
	DataVector InterpConfigEuler(const DataVector& q1, const DataVector& q2, double t);
	
	DataVector InterpConfigQuat(const DataVector& q1, const DataVector& q2, double t);
	
	DataVector ConfigQuat2Euler(const DataVector& q);
	
	DataVector ConfigEuler2Quat(const DataVector& q);
	
	
	
	
	
	void ConfigEuler2RotTran(double R[3][3], double T[3], const DataVector& q);
	
	void ConfigQuat2RotTrans(double R[3][3], double T[3], const DataVector& q);
	
	void IdentityRotTrans(double R[3][3], double T[3]);

	DataVector relative2D(const DataVector& q1, const DataVector& q2);

	DataVector relative3D(const DataVector& q1, const DataVector& q2);
	
}

#endif