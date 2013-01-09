#ifndef DISTANCE_PROXY_H
#define DISTANCE_PROXY_H

#include "data_vector.h"
#include "math_utility.h"
#include <iostream>
#include <flann/flann.hpp>



namespace APDL
{
	extern double distance_weight[7]; // must be set for PDg

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
				for(size_t i = 0; i < size - 1; ++i) 
				{
					diff = *a++ - *b++;
					result += diff*diff * distance_weight[i];
				}
				
				diff = angleTruncate(*a - *b);
				result +=  diff * diff * distance_weight[size - 1];

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
					result += diff*diff * distance_weight[i];
				}
				
				for(size_t i = 0; i < 3; ++i)
				{
					diff = *a++ - *b++;
					diff = angleTruncate(diff);
					result += diff * diff * distance_weight[size - 3 + i];
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
				for(size_t i = 0; i < size - 4; ++i) 
				{
					diff = *a++ - *b++;
					result += diff*diff * distance_weight[i];
				}

				Quaternion q1(a[0], a[1], a[2], a[3]);
				Quaternion q2(b[0], b[1], b[2], b[3]);

				Quaternion q = q1 * q2.inverse();

				result += (q[1] * q[1] * distance_weight[size - 4] + q[2] * q[2] * distance_weight[size - 3] + q[3] * q[3] * distance_weight[size - 2]);

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
		DistanceProxySE2()
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double a = v1[0] - v2[0];
			double b = v1[1] - v2[1];
			double c = angleTruncate(v1[2] - v2[2]);
			
			return a * a * distance_weight[0] + b * b * distance_weight[1] + c * c * distance_weight[2];
		}
		
		std::size_t dim() const
		{
			return 3;
		}
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
		DistanceProxySE3EulerAngle()
		{}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double dx = v1[0] - v2[0];
			double dy = v1[1] - v2[1];
			double dz = v1[2] - v2[2];
			
			double da = angleTruncate(v1[3] - v2[3]);
			double db = angleTruncate(v1[4] - v2[4]);
			double dc = angleTruncate(v1[5] - v2[5]);
			
			return dx * dx * distance_weight[0] + dy * dy * distance_weight[1] + dz * dz * distance_weight[2] + 
				(da * da * distance_weight[3] + db * db * distance_weight[4] + dc * dc * distance_weight[5]);
		}
		
		std::size_t dim() const
		{
			return 6;
		}
	};
	
	
	class DistanceProxySE3Quaternion : public DistanceProxyBase
	{
	public:
		DistanceProxySE3Quaternion() {}
		
		double sqrDistance(const DataVector& v1, const DataVector& v2) const
		{
			double dx = v1[0] - v2[0];
			double dy = v1[1] - v2[1];
			double dz = v1[2] - v2[2];
			
			Quaternion q1(v1[3], v1[4], v1[5], v1[6]);
			Quaternion q2(v2[3], v2[4], v2[5], v2[6]);
			
			Quaternion q = q1 * q2.inverse();
			
			
			double d = dx * dx * distance_weight[0] + dy * dy * distance_weight[1] + dz * dz * distance_weight[2] + 
				(q[1] * q[1] * distance_weight[3] + q[2] * q[2] * distance_weight[4] + q[3] * q[3] * distance_weight[5]);
			return d;
		}
		
		std::size_t dim() const
		{
			return 7;
		}
	};
	
}

#endif