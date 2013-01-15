#ifndef MATH_UTILITY_H
#define MATH_UTILITY_H

#include <boost/math/constants/constants.hpp>
#include "data_vector.h"
#include <cassert>

namespace APDL
{
	class Vec2D
	{
	public:
		Vec2D()
		{
			x = 0;
			y = 0;
		}
		Vec2D(double x_, double y_) : x(x_), y(y_) {}
		
		void set(double x_, double y_)
		{
			x = x_;
			y = y_;
		}
		
		Vec2D operator - () const
		{
			return Vec2D(-x, -y);
		}
		
		Vec2D& operator += (const Vec2D& v)
		{
			x += v.x;
			y += v.y;
			return *this;
		}
		
		Vec2D& operator -= (const Vec2D& v)
		{
			x -= v.x;
			y -= v.y;
			return *this;
		}
		
		Vec2D& operator *= (double d)
		{
			x *= d;
			y *= d;
			return *this;
		}
		
		Vec2D operator * (double d) const
		{
			return Vec2D(x * d, y * d);
		}
		
		Vec2D operator + (const Vec2D& v) const
		{
			return Vec2D(x + v.x, y + v.y);
		}
		
		Vec2D operator - (const Vec2D& v) const
		{
			return Vec2D(x - v.x, y - v.y);
		}
		
		double sqrLength() const
		{
			return x * x + y * y;
		}
		
		double length() const
		{
			return sqrt(x * x + y * y);
		}
		
		Vec2D& normalize()
		{
			double len = sqrt(x * x + y * y);
			if(len > 0)
			{
				x /= len;
				y /= len;
			}
			
			return *this;
		}
		
		double x, y;
	};
	
	inline Vec2D normalize(const Vec2D& v)
	{
		Vec2D v_(v);
		v_.normalize();
		return v_;
	}
	
	inline Vec2D operator * (double s, const Vec2D& v)
	{
		return Vec2D(v.x * s, v.y * s);
	}
	
	inline double cross(const Vec2D& a, const Vec2D& b)
	{
		return a.x * b.y - a.y * b.x;
	}
	
	inline Vec2D cross(const Vec2D& a, double s)
	{
		return Vec2D(s * a.y, -s * a.x);
	}
	
	inline Vec2D cross(double s, const Vec2D& a)
	{
		return Vec2D(-s * a.y, s * a.x);
	}
	
	inline double dot(const Vec2D& a, const Vec2D& b)
	{
		return a.x * b.x + a.y * b.y;
	}
	
	inline double sqrDistance(const Vec2D& a, const Vec2D& b)
	{
		return (b - a).sqrLength();
	}
	
	inline double distance(const Vec2D& a, const Vec2D& b)
	{
		return (b - a).length();
	}
	
	inline std::ostream& operator << (std::ostream& os, const Vec2D& v)
	{
		os << v.x << " " << v.y;
		return os;
	}
	
	
	struct Mat2D
	{
		Mat2D()
		{
			col[0] = Vec2D(1, 0);
			col[1] = Vec2D(0, 1);
		}
		
		Mat2D(double angle)
		{
			double c = cos(angle);
			double s = sin(angle);
			col[0] = Vec2D(c, s);
			col[1] = Vec2D(-s, c);
		}
		
		Mat2D(const Vec2D& col1, const Vec2D& col2)
		{
			col[0] = col1;
			col[1] = col2;
		}
		
		Mat2D(double a, double b, double c, double d)
		{
			col[0] = Vec2D(a, c);
			col[1] = Vec2D(b, d);
		}
		
		Mat2D& transpose()
		{
			double temp = col[0].y;
			col[0].y = col[1].x;
			col[1].x = temp;
			return *this;
		}
		
		Mat2D& inverse()
		{
			double a = col[0].x;
			double b = col[1].x;
			double c = col[0].y;
			double d = col[1].y;
			
			double det = a * d - b * c;
			assert(det != 0);
			det = 1.0 / det;
			col[0].x = det * d;
			col[0].y = -det * c;
			col[1].x = -det * b;
			col[1].y = det * a;
			
			return *this;
		}
		
		void set(double angle)
		{
			double c = cos(angle);
			double s = sin(angle);
			col[0].x = c;
			col[0].y = s;
			col[1].x = -s;
			col[1].y = c;
		}
		
		Vec2D operator * (const Vec2D& v) const
		{
			return Vec2D(col[0].x * v.x + col[1].x * v.y, col[0].y * v.x + col[1].y * v.y);
		}
		
		Vec2D col[2];
	};
	
	inline std::ostream& operator << (std::ostream& os, const Mat2D& m)
	{
		os << m.col[0].x << " " << m.col[1].x << std::endl;
		os << m.col[0].y << " " << m.col[1].y;
		return os;
	}
	
	
	inline Mat2D transpose(const Mat2D& m)
	{
		return Mat2D(m.col[0].x, m.col[0].y, m.col[1].x, m.col[1].y);
	}
	
	inline Mat2D inverse(const Mat2D& m)
	{
		Mat2D r(m);
		r.inverse();
		return r;
	}
	
	
	struct Transform2D
	{
		Mat2D R;
		Vec2D t;
		
		Transform2D()
		{
			R = Mat2D(1, 0, 0, 1);
			t = Vec2D(0, 0);
		}
		
		Transform2D(const Mat2D& R_, const Vec2D& t_) : R(R_), t(t_)
		{}
		
		Vec2D transform(const Vec2D& v) const
		{
			return R * v + t;
		}
		
		Vec2D untransform(const Vec2D& v) const
		{
			return inverse(R) * (v - t);
		}
		
		Vec2D transform_dir(const Vec2D& v) const
		{
			return R * v;
		}
		
		Vec2D untransform_dir(const Vec2D& v) const
		{
			return inverse(R) * v;
		}
	};
	

	class Quaternion
	{
		double data[4];
	public:
		Quaternion(double a, double b, double c, double d, bool normalize = true)
		{
			if(normalize)
			{
				double inv_len = 1.0 / sqrt(a * a + b * b + c * c + d * d);
				data[0] = a * inv_len;
				data[1] = b * inv_len;
				data[2] = c * inv_len;
				data[3] = d * inv_len;
			}
			else
			{
				data[0] = a;
				data[1] = b;
				data[2] = c;
				data[3] = d;
			}
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
			                  data[0] * other.data[3] + data[1] * other.data[2] - data[2] * other.data[1] + data[3] * other.data[0], false);
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

	void Euler2RotZYX(double R[3][3], double a, double b, double c);
	
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

	DataVector unrelative3D(const DataVector& q1, const DataVector q);
	
}

#endif