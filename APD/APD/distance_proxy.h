#ifndef DISTANCE_PROXY_H
#define DISTANCE_PROXY_H

#include "data_vector.h"
#include "math_utility.h"
#include <iostream>
#include <flann/flann.hpp>



namespace APDL
{
	extern double angle_weight;
	namespace FLANN_WRAPPER
	{
		struct DistanceRN
		{
			typedef bool is_kdtree_distance;

			typedef double ElementType;
			typedef double ResultType;

			template <typename Iterator1, typename Iterator2>
			ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
			{
				ResultType result = ResultType();
				ResultType diff;
				for(size_t i = 0; i < size; ++i ) 
				{
					diff = *a++ - *b++;
					result += diff*diff;
				}
				return result;
			}

			template <typename U, typename V>
			inline ResultType accum_dist(const U& a, const V& b, int) const
			{
				return (a-b)*(a-b);
			}
		};


		struct DistanceSE2
		{
			typedef double ElementType;
			typedef double ResultType;

			DistanceSE2()
			{}

			/// size must be 3
			template <typename Iterator1, typename Iterator2>
			ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
			{
				ResultType result = ResultType();
				ResultType diff;
				for(size_t i = 0; i < size - 1; ++i ) 
				{
					diff = *a++ - *b++;
					result += diff*diff;
				}
				
				diff = angleTruncate(*a - *b);
				result += angle_weight * diff * diff;

				return result;
			}
		};

		struct DistanceSE3EulerAngle
		{
			typedef double ElementType;
			typedef double ResultType;

			DistanceSE3EulerAngle()
			{}

			/// size must be 6
			template <typename Iterator1, typename Iterator2>
			ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
			{
				ResultType result = ResultType();
				ResultType diff;
				for(size_t i = 0; i < size - 3; ++i ) 
				{
					diff = *a++ - *b++;
					result += diff*diff;
				}
				
				for(size_t i = 0; i < 3; ++i)
				{
					diff = *a++ - *b++;
					diff = angleTruncate(diff);
					result += angle_weight * diff * diff;
				}

				return result;
			}
		};


		struct DistanceSE3Quat
		{
			typedef double ElementType;
			typedef double ResultType;

			DistanceSE3Quat()
			{}

			/// size must be 7
			template <typename Iterator1, typename Iterator2>
			ResultType operator()(Iterator1 a, Iterator2 b, size_t size, ResultType /*worst_dist*/ = -1) const
			{
				ResultType result = ResultType();
				ResultType diff;
				for(size_t i = 0; i < size - 4; ++i ) 
				{
					diff = *a++ - *b++;
					result += diff*diff;
				}

				Quaternion q1(a[0], a[1], a[2], a[3]);
				Quaternion q2(b[0], b[1], b[2], b[3]);

				Quaternion q = q1 * q2.inverse();

				result += angle_weight * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);

				return result;
			}
		};
	}

	class DistanceProxyBase
	{
	public:
		virtual double sqrDistance(const DataVector& v1, const DataVector& v2) const = 0;
		virtual std::size_t dim() const = 0;
	};
	
	class DistanceProxyRN : public DistanceProxyBase
	{
	public:
		DistanceProxyRN(const DataVector& weight_) : weight(weight_)
		{}

		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double result = 0;
			for(std::size_t i = 0; i < v1.dim(); ++i)
			{
				double a = v1[i] - v2[i];
				result += a * a * weight[i];
			}

			return result;
		}

		std::size_t dim() const
		{
			return weight.dim();
		}

		DataVector weight;
	};

	class DistanceProxyR2 : public DistanceProxyBase
	{
	public:
		DistanceProxyR2()
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double a = v1[0] - v2[0];
			double b = v1[1] - v2[1];			
			return a * a + b * b;
		}
		
		std::size_t dim() const
		{
			return 2;
		}
	};
	
	
	class DistanceProxySE2 : public DistanceProxyBase
	{
	public:
		DistanceProxySE2(double angle_weight_ = 1) : angle_weight(angle_weight_)
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double a = v1[0] - v2[0];
			double b = v1[1] - v2[1];
			double c = angleTruncate(v1[2] - v2[2]);
			
			return a * a + b * b + angle_weight * c * c;
		}
		
		std::size_t dim() const
		{
			return 3;
		}
		
	protected:
	
		double angle_weight;
	};
	
	class DistanceProxyR3 : public DistanceProxyBase
	{
	public:
		DistanceProxyR3()
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double a = v1[0] - v2[0];
			double b = v1[1] - v2[1];	
			double c = v1[2] - v2[2];
			return a * a + b * b + c * c;
		}
		
		std::size_t dim() const
		{
			return 3;
		}
	};
	
	class DistanceProxySE3EulerAngle : public DistanceProxyBase
	{
	public:
		DistanceProxySE3EulerAngle(double angle_weight_ = 1) : angle_weight(angle_weight_)
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double dx = v1[0] - v2[0];
			double dy = v1[1] - v2[1];
			double dz = v1[2] - v2[2];
			
			double da = angleTruncate(v1[3] - v2[3]);
			double db = angleTruncate(v1[4] - v2[4]);
			double dc = angleTruncate(v1[5] - v2[5]);
			
			return dx * dx + dy * dy + dz * dz + angle_weight * (da * da + db * db + dc * dc);
		}
		
		std::size_t dim() const
		{
			return 6;
		}
		
	protected:
	
		double angle_weight;
	};
	
	
	class DistanceProxySE3Quaternion : public DistanceProxyBase
	{
	public:
		DistanceProxySE3Quaternion(double q_weight_) : q_weight(q_weight_) {}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double dx = v1[0] - v2[0];
			double dy = v1[1] - v2[1];
			double dz = v1[2] - v2[2];
			
			Quaternion q1(v1[3], v1[4], v1[5], v1[6]);
			Quaternion q2(v2[3], v2[4], v2[5], v2[6]);
			
			Quaternion q = q1 * q2.inverse();
			
			
			double d = dx * dx + dy * dy + dz * dz + q_weight * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
			return d;
		}
		
		std::size_t dim() const
		{
			return 7;
		}
		
	protected:
	
		double q_weight;
	};
	
}

#endif