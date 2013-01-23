#include "math_utility.h"
#include <iostream>

namespace APDL
{

	void Quat2Rot(double R[3][3], const Quaternion& quat)
	{
		double twoX  = 2.0 * quat[1];
		double twoY  = 2.0 * quat[2];
		double twoZ  = 2.0 * quat[3];
		double twoWX = twoX * quat[0];
		double twoWY = twoY * quat[0];
		double twoWZ = twoZ * quat[0];
		double twoXX = twoX * quat[1];
		double twoXY = twoY * quat[1];
		double twoXZ = twoZ * quat[1];
		double twoYY = twoY * quat[2];
		double twoYZ = twoZ * quat[2];
		double twoZZ = twoZ * quat[3];
		
		R[0][0] = 1.0 - (twoYY + twoZZ);
		R[0][1] = twoXY - twoWZ;
		R[0][2] = twoXZ + twoWY;
		R[1][0] = twoXY + twoWZ;
		R[1][1] = 1.0 - (twoXX + twoZZ);
		R[1][2] = twoYZ - twoWX;
		R[2][0] = twoXZ - twoWY;
		R[2][1] = twoYZ + twoWX;
		R[2][2] = 1.0 - (twoXX + twoYY);
	}
	
	void Quat2Euler(double& a, double& b, double& c, const Quaternion& quat)
	{
		double R[3][3];
		
		Quat2Rot(R, quat);
		//
		//if(R[0][2] < 1.0)
		//{
		//	if(R[0][2] > -1.0)
		//	{
		//		b = asin(R[0][2]);
		//		a = atan2(-R[1][2], R[2][2]);
		//		c = atan2(-R[0][1], R[0][0]);
		//	}
		//	else
		//	{
		//		b = -0.5 * boost::math::constants::pi<double>();
		//		a = -atan2(R[1][0], R[1][1]);
		//		c = 0.0;
		//	}
		//}
		//else
		//{
		//	b = boost::math::constants::pi<double>();
		//	a = atan2(R[1][0], R[1][1]);
		//	c = 0.0;
		//}
		a = atan2(R[1][0], R[0][0]);
		b = asin(-R[2][0]);
		c = atan2(R[2][1], R[2][2]);

		if(b == boost::math::constants::pi<double>() * 0.5)
		{
			if(a > 0)
				a -= boost::math::constants::pi<double>();
			else 
				a += boost::math::constants::pi<double>();

			if(c > 0)
				c -= boost::math::constants::pi<double>();
			else
				c += boost::math::constants::pi<double>();
		}

	}
	
	void Euler2Quat(Quaternion& quat, double a, double b, double c)
	{
		double R[3][3];
		Euler2Rot(R, a, b, c);
		
		// R to quaternion
		Rot2Quat(quat, R);
	}
	
	void Rot2Quat(Quaternion& quat, double R[3][3])
	{
		double a, b, c, d;
		double trace = R[0][0] + R[1][1] + R[2][2];
		double root;
		if(trace > 0.0)
		{
			root = sqrt(trace + 1.0);
			a = 0.5 * root;
			root = 0.5 / root;
			b = (R[2][1] - R[1][2]) * root;
			c = (R[0][2] - R[2][0]) * root;
			d = (R[1][0] - R[0][1]) * root;
		}
		else
		{
			int i = 0;
			if (R[1][1] > R[0][0])
				i = 1;
			if (R[2][2] > R[i][i])
				i = 2;

			int j, k;
			if (i == 0)
			{
				j = 1;
				k = 2;
			}
			else if (i == 1)
			{
				j = 2;
				k = 0;
			}
			else if (i == 2)
			{
				j = 0;
				k = 1;
			}

			root = sqrt(R[i][i] - R[j][j] - R[k][k] + 1);

			if (i == 0)
			{
				b = 0.5 * root;
				root = 0.5 / root;
				a = (R[k][j] - R[j][k]) * root;
				c = (R[j][i] + R[i][j]) * root;
				d = (R[k][i] + R[i][k]) * root;
			}
			else if (i == 1)
			{
				c = 0.5 * root;
				root = 0.5 / root;
				a = (R[k][j] - R[j][k]) * root;
				d = (R[j][i] + R[i][j]) * root;
				b = (R[k][i] + R[i][k]) * root;
			}
			else if (i == 2)
			{
				d = 0.5 * root;
				root = 0.5 / root;
				a = (R[k][j] - R[j][k]) * root;
				b = (R[j][i] + R[i][j]) * root;
				c = (R[k][i] + R[i][k]) * root;
			}
		}

		quat = Quaternion(a, b, c, d);
	}

	void Rot2Euler(double& a, double& b, double& c, double R[3][3])
	{
		Quaternion q;
		Rot2Quat(q, R);
		Quat2Euler(a, b, c, q);
	}

	double angleTruncate(double angle)
	{
		if(angle > 0)
		{
			if(angle > boost::math::constants::pi<double>())
			{
				int n = (int)std::floor(angle / (2 * boost::math::constants::pi<double>()));
				angle -= n * 2 * boost::math::constants::pi<double>();
				
				if(angle > boost::math::constants::pi<double>())
					angle -= 2 * boost::math::constants::pi<double>();
			}
		}
		else if(angle < 0)
		{
			if(angle < -boost::math::constants::pi<double>())
			{
				int n = (int)std::floor(-angle / (2 * boost::math::constants::pi<double>()));
				angle += n * 2 * boost::math::constants::pi<double>();
				
				if(angle < -boost::math::constants::pi<double>())
					angle += 2 * boost::math::constants::pi<double>();
			}
		}
		
		return angle;
	}
	
	
	DataVector InterpConfig2D(const DataVector& q1, const DataVector& q2, double t)
	{
		DataVector q(3);
		q[0] = (1 - t) * q1[0] + t * q2[0];
		q[1] = (1 - t) * q1[1] + t * q2[1];
		
		double a0 = angleTruncate(q1[2]);
		double a1 = angleTruncate(q2[2]);
		if(std::fabs(a0 - a1) <= boost::math::constants::pi<double>())
			q[2] = (1 - t) * a0 + t * a1;
		else if(a0 >= 0)
			q[2] = angleTruncate((1 - t) * a0 + t * (a1 + 2 * boost::math::constants::pi<double>()));
		else
			q[2] = angleTruncate((1 - t) * (a0 + 2 * boost::math::constants::pi<double>()) + t * a1);
			
		return q;
	}
	
	
	DataVector InterpConfigEuler(const DataVector& q1, const DataVector& q2, double t)
	{
		DataVector q(q1.dim());
		q[0] = (1 - t) * q1[0] + t * q2[0];
		q[1] = (1 - t) * q1[1] + t * q2[1];
		q[2] = (1 - t) * q1[2] + t * q2[2];
		
		double a0 = angleTruncate(q1[3]);
		double a1 = angleTruncate(q2[3]);
		if(std::fabs(a0 - a1) <= boost::math::constants::pi<double>())
			q[3] = (1 - t) * a0 + t * a1;
		else if(a0 >= 0)
			q[3] = angleTruncate((1 - t) * a0 + t * (a1 + 2 * boost::math::constants::pi<double>()));
		else
			q[3] = angleTruncate((1 - t) * (a0 + 2 * boost::math::constants::pi<double>()) + t * a1);
			
		a0 = angleTruncate(q1[4]);
		a1 = angleTruncate(q2[4]);
		if(std::fabs(a0 - a1) <= boost::math::constants::pi<double>())
			q[4] = (1 - t) * a0 + t * a1;
		else if(a0 >= 0)
			q[4] = angleTruncate((1 - t) * a0 + t * (a1 + 2 * boost::math::constants::pi<double>()));
		else
			q[4] = angleTruncate((1 - t) * (a0 + 2 * boost::math::constants::pi<double>()) + t * a1);
			
		a0 = angleTruncate(q1[5]);
		a1 = angleTruncate(q2[5]);
		if(std::fabs(a0 - a1) <= boost::math::constants::pi<double>())
			q[5] = (1 - t) * a0 + t * a1;
		else if(a0 >= 0)
			q[5] = angleTruncate((1 - t) * a0 + t * (a1 + 2 * boost::math::constants::pi<double>()));
		else
			q[5] = angleTruncate((1 - t) * (a0 + 2 * boost::math::constants::pi<double>()) + t * a1);
			
		return q;
	}
	
	DataVector InterpConfigQuat(const DataVector& q1, const DataVector& q2, double t)
	{
		DataVector q(q1.dim());
		q[0] = (1 - t) * q1[0] + t * q2[0];
		q[1] = (1 - t) * q1[1] + t * q2[1];
		q[2] = (1 - t) * q1[2] + t * q2[2];
		
		double cs = q1[3] * q2[3] + q1[4] * q2[4] + q1[5] * q2[5] + q1[6] * q2[6];
		double angle = std::acos(cs);
		
		if(std::fabs(angle) > 0.0)
		{
			double sn = sin(angle);
			double inv_sn = 1.0 / sn;
			double t_angle = t * angle;
			double coeff0 = sin(angle - t_angle) * inv_sn;
			double coeff1 = sin(t_angle) * inv_sn;
			
			q[3] = coeff0 * q1[3] + coeff1 * q2[3];
			q[4] = coeff0 * q1[4] + coeff1 * q2[4];
			q[5] = coeff0 * q1[5] + coeff1 * q2[5];
			q[6] = coeff0 * q1[6] + coeff1 * q2[6];
		}
		else
		{
			q[3] = q1[3];
			q[4] = q1[4];
			q[5] = q1[5];
			q[6] = q1[6];
		}
		
		return q;
	}
	
	DataVector ConfigQuat2Euler(const DataVector& q)
	{
		DataVector qout(6);
		qout[0] = q[0];
		qout[1] = q[1];
		qout[2] = q[2];
		
		double a, b, c;
		Quat2Euler(a, b, c, Quaternion(q[3], q[4], q[5], q[6]));
		
		qout[3] = a;
		qout[4] = b;
		qout[5] = c;
		
		return qout;
	}
	
	
	
	DataVector ConfigEuler2Quat(const DataVector& q)
	{
		DataVector qout(7);
		qout[0] = q[0];
		qout[1] = q[1];
		qout[2] = q[2];
		
		Quaternion quat;
		Euler2Quat(quat, q[3], q[4], q[5]);
		qout[3] = quat[0];
		qout[4] = quat[1];
		qout[5] = quat[2];
		qout[6] = quat[3];
		
		return qout;
	}
	
	
	void Euler2Rot(double R[3][3], double a, double b, double c)
	{
		Euler2RotZYX(R, c, b, a);
	}

	void Euler2RotZYX(double R[3][3], double a, double b, double c)
	{

		double ci(cos(a));
		double cj(cos(b));
		double ch(cos(c));
		double si(sin(a));
		double sj(sin(b));
		double sh(sin(c));
		double cc = ci * ch;
		double cs = ci * sh;
		double sc = si * ch;
		double ss = si * sh;
		
		R[0][0] = cj * ch;
		R[0][1] = sj * sc - cs;
		R[0][2] = sj * cc + ss;
		R[1][0] = cj * sh;
		R[1][1] = sj * ss + cc;
		R[1][2] = sj * cs - sc;
		R[2][0] = -sj;
		R[2][1] = cj * si;
		R[2][2] = cj * ci;
	}
	
	
	void ConfigEuler2RotTran(double R[3][3], double T[3], const DataVector& q)
	{
		Euler2Rot(R, q[3], q[4], q[5]);
		T[0] = q[0];
		T[1] = q[1];
		T[2] = q[2];
	}
	
	
	
	void ConfigQuat2RotTrans(double R[3][3], double T[3], const DataVector& q)
	{
		Quat2Rot(R, Quaternion(q[3], q[4], q[5], q[6]));
		
		T[0] = q[0];
		T[1] = q[1];
		T[2] = q[2];
	}
	
	void IdentityRotTrans(double R[3][3], double T[3])
	{
		R[0][0] = 1;
		R[0][1] = 0;
		R[0][2] = 0;
		R[1][0] = 0;
		R[1][1] = 1;
		R[1][2] = 0;
		R[2][0] = 0;
		R[2][1] = 0;
		R[2][2] = 1;
		
		T[0] = 0;
		T[1] = 0;
		T[2] = 0;
	}

	DataVector relative2D(const DataVector& q1, const DataVector& q2)
	{
		//DataVector q(3);
		//q[0] = q2[0] - q1[0];
		//q[1] = q2[1] - q1[1];
		//q[2] = angleTruncate(q2[2] - q1[2]);
		//return q;

		DataVector q(3);

		Mat2D m(-q1[2]);
		Vec2D v(q2[0] - q1[0], q2[1] - q1[1]);
		v = m * v;
		q[0] = v.x;
		q[1] = v.y;
		q[2] = angleTruncate(q2[2] - q1[2]);
		
		return q;
	}

	DataVector unrelative3D(const DataVector& q1, const DataVector q)
	{
		if(q1.dim() == 6)
		{
			DataVector q2(6);

			double t1 = q1[0];
			double t2 = q1[1];
			double t3 = q1[2];

			Quaternion quat1, quat;
			Euler2Quat(quat1, q1[3], q1[4], q1[5]);
			Euler2Quat(quat, q[3], q[4], q[5]);

			Quaternion r = quat1 * Quaternion(0, t1, t2, t3, false) * inverse(quat1);
			q2[0] = r[1] + q[0];
			q2[1] = r[2] + q[1];
			q2[2] = r[3] + q[2];

			Quaternion quat2 = quat1 * quat;

			double a, b, c;
			Quat2Euler(a, b, c, quat2);
			q2[3] = a;
			q2[4] = b;
			q2[5] = c;

			return q2;

		}
		else if(q1.dim() == 7)
		{
			DataVector q2(7);

			double t1 = q1[0];
			double t2 = q1[1];
			double t3 = q1[2];

			Quaternion quat1(q1[3], q1[4], q1[5], q1[6]);
			Quaternion quat(q[3], q[4], q[5], q[6]);

			Quaternion r = quat1 * Quaternion(0, t1, t2, t3, false) * inverse(quat1);
			q2[0] = r[1] + q[0];
			q2[1] = r[2] + q[1];
			q2[2] = r[3] + q[2];

			Quaternion quat2 = quat1 * quat;

			q2[3] = quat2[0];
			q2[4] = quat2[1];
			q2[5] = quat2[2];
			q2[6] = quat2[3];

			return q2;
		}
	}

	DataVector relative3D(const DataVector& q1, const DataVector& q2)
	{
		if(q1.dim() == 6)
		{
			DataVector q(6);
			double t1 = q2[0] - q1[0];
			double t2 = q2[1] - q1[1];
			double t3 = q2[2] - q1[2];

			Quaternion quat1, quat2;
			Euler2Quat(quat1, q1[3], q1[4], q1[5]);
			Euler2Quat(quat2, q2[3], q2[4], q2[5]);


			Quaternion r = inverse(quat1) * Quaternion(0, t1, t2, t3, false) * quat1;

			q[0] = r[1];
			q[1] = r[2];
			q[2] = r[3];

			Quaternion quat = inverse(quat1) * quat2;

			double a, b, c;
			Quat2Euler(a, b, c, quat);
			q[3] = a;
			q[4] = b;
			q[5] = c;

			return q;
		}
		else if(q2.dim() == 7)
		{
			DataVector q(7);
			double t1 = q2[0] - q1[0];
			double t2 = q2[1] - q1[1];
			double t3 = q2[2] - q1[2];

			Quaternion quat1(q1[3], q1[4], q1[5], q1[6]);
			Quaternion quat2(q2[3], q2[4], q2[5], q2[6]);

			Quaternion r = inverse(quat1) * Quaternion(0, t1, t2, t3) * quat1;
			q[0] = r[1];
			q[1] = r[2];
			q[2] = r[3];

			Quaternion quat = inverse(quat1) * quat2;

			q[3] = quat[0];
			q[4] = quat[1];
			q[5] = quat[2];
			q[6] = quat[3];

			return q;
		}
		else
		{
			std::cout << "dim not correct" << std::endl;
			return DataVector();
		}
	}
	
}