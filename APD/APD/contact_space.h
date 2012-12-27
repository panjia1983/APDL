#ifndef CONTACT_SPACE_H
#define CONTACT_SPACE_H

#include "collider.h"
#include "sampler.h"
#include <vector>
#include "distance_proxy.h"
#include <C2A/LinearMath.h>
#include "libm3d_wrapper.h"

namespace APDL
{
	inline std::vector<DataVector> readPM3dFile(const std::string& filename)
	{
		std::vector<DataVector> points;

		FILE* fp = fopen(filename.c_str(), "rb");
		if(fp == NULL)
		{
			std::cerr << "!Error: Cannot open file (" << filename << ")" << std::endl;
			return points;
		}

		float radian[3];
		int total_bdSize;
		fread(&total_bdSize, sizeof(int), 1, fp);
		fread(radian, sizeof(float) * 3, 1, fp);

		float c_r = cos(radian[0]);
		float s_r = sin(radian[0]);

		Matrix3x3 mx(1,0,0, 0, c_r,-s_r, 0, s_r,c_r);
		c_r = cos(radian[1]);
		s_r = sin(radian[1]);
		Matrix3x3 my(c_r, 0, s_r, 0, 1, 0, -s_r, 0, c_r);
		c_r = cos(radian[2]);
		s_r = sin(radian[2]);
		Matrix3x3 mz(c_r,-s_r, 0, s_r,c_r,0,0,0,1);
		Matrix3x3 M = mz * my * mx;

		float* bdM = new float[6*total_bdSize];
		fread(bdM, sizeof(float)*6*total_bdSize, 1, fp);

		fclose(fp);

		
		DataVector p(3);
		for(std::size_t i = 0; i < total_bdSize; ++i)
		{
			p[0] = bdM[6*i];
			p[1] = bdM[6*i+1];
			p[2] = bdM[6*i+2];
			points.push_back(p);
		}

		delete [] bdM;

		return points;
	}


	struct ContactSpaceSampleData
	{
		DataVector v;
		bool col;
		
		ContactSpaceSampleData(const DataVector& v_, bool col_) : v(v_), col(col_) {}
	};

	struct Scaler
	{
		DataVector v_min, v_max;

		std::size_t active_dim;

		double scale_mag;

		Scaler(const DataVector& v_min_, const DataVector& v_max_, std::size_t active_dim_) : v_min(v_min_), v_max(v_max_), active_dim(active_dim_)
		{
			scale_mag = 0;
			for(std::size_t i = 0; i < active_dim; ++i)
				scale_mag += (v_min[i] - v_max[i]) * (v_min[i] - v_max[i]);
			scale_mag = sqrt(scale_mag);
		}

		DataVector scale(const DataVector& v) const
		{
			DataVector res(v);
			for(std::size_t i = 0; i < active_dim; ++i)
			{
				res[i] = (v[i] - v_min[i]) / (v_max[i] - v_min[i]);
			}

			return res;
		}

		DataVector unscale(const DataVector& v) const
		{
			DataVector res(v);
			for(std::size_t i = 0; i < active_dim; ++i)
			{
				res[i] = v[i] * (v_max[i] - v_min[i]) + v_min[i];
			}

			return res;
		}

		double getScale() const
		{
			return scale_mag;
		}
	};

	inline std::ofstream& operator << (std::ofstream& os, const Scaler& scaler)
	{
		for(std::size_t i = 0; i < scaler.active_dim; ++i)
		{
			os << 1.0 / (scaler.v_max[i] - scaler.v_min[i]) << " " << - scaler.v_min[i] / (scaler.v_max[i] - scaler.v_min[i]) << std::endl;
		}
		return os;
	}

	template<typename Collider>
	bool checkStatus(const Collider& collider, std::size_t dim, const DataVector& query)
	{
		std::size_t in_dim = query.dim();

		DataVector v(dim);
		if(in_dim < dim)
		{
			for(std::size_t j = 0; j < in_dim; ++j)
				v[j] = query[j];
			for(std::size_t j = in_dim; j < dim; ++j)
				v[j] = 0;
		}
		else
		{
			for(std::size_t j = 0; j < dim; ++j)
				v[j] = query[j];
		}

		return collider.isCollide(v);
	}

	template<typename Collider>
	std::vector<bool> checkStatus(const Collider& collider, std::size_t dim, const std::vector<DataVector>& queries)
	{
		std::vector<bool> statuses(queries.size());

		for(std::size_t i = 0; i < queries.size(); ++i)
			statuses[i] = checkStatus(collider, dim, queries[i]);

		return statuses;
	}
	
	class ContactSpaceR2
	{
	public:
		typedef FLANN_WRAPPER::DistanceRN DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider2D ColliderType;

		ContactSpaceR2(const Polygon& p1, const Polygon& p2, double delta = 0) : collider(&p1, &p2)
		{
			AABB2D aabb1 = p1.getAABB();
			AABB2D aabb2 = p2.getAABB();
			
			double delta_x_min = aabb1.x_min - aabb2.x_max - delta;
			double delta_x_max = aabb1.x_max - aabb2.x_min + delta;
			double delta_y_min = aabb1.y_min - aabb2.y_max - delta;
			double delta_y_max = aabb1.y_max - aabb2.y_min + delta;
			sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max);
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			DataVector v(3);
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v_ = sampler.sample();
			
				v[0] = v_[0];
				v[1] = v_[1];
				v[2] = 0;
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			DataVector v(3);
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v_ = sampler.sample();
			
				v[0] = v_[0];
				v[1] = v_[1];
				v[2] = 0;
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 3;
		}
		
		std::size_t active_data_dim() const
		{
			return 2;
		}

		Scaler getScaler() const
		{
			double x_min, x_max, y_min, y_max;
			sampler.getBound(x_min, x_max, y_min, y_max);
			DataVector v_min(2);
			DataVector v_max(2);
			v_min[0] = x_min; v_min[1] = y_min;
			v_max[0] = x_max; v_max[1] = y_max;
			return Scaler(v_min, v_max, 2);
		}
		
		SamplerR2 sampler;
		
		Collider2D collider;
	};
	
	class ContactSpaceSE2_disk
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE2 DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider2D ColliderType;

		ContactSpaceSE2_disk(const Polygon& p1, const Polygon& p2, double delta = 0) : collider(&p1, &p2)
		{
			std::pair<Vec2D, double> c1 = p1.getCircle();
			std::pair<Vec2D, double> c2 = p2.getCircle();
			sampler.setBound(c1.first.x, c1.first.y , c1.second + 0.5 * delta,
			                 c2.first.x, c2.first.y, c2.second + 0.5 * delta);
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 3;
		}
		
		std::size_t active_data_dim() const
		{
			return 3;
		}

		Scaler getScaler() const
		{
			double x_center, y_center, r1, r2;
			sampler.getBound(x_center, y_center, r1, r2);
			DataVector v_min(3);
			DataVector v_max(3);
			v_min[0] = x_center - r2; v_min[1] = y_center - r2; v_min[2] = -boost::math::constants::pi<double>();
			v_max[0] = x_center + r2; v_max[1] = y_center + r2; v_max[2] = boost::math::constants::pi<double>();
			
			return Scaler(v_min, v_max, 3);
		}
		
		SamplerSE2_disk sampler;
		
		Collider2D collider;
	};

	class ContactSpaceSE2
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE2 DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider2D ColliderType;

		ContactSpaceSE2(const Polygon& p1, const Polygon& p2, double delta = 0) : collider(&p1, &p2)
		{
			std::pair<Vec2D, double> c1 = p1.getCircle();
			std::pair<Vec2D, double> c2 = p2.getCircle();

			SamplerSE2_disk sampler_;
			sampler_.setBound(c1.first.x, c1.first.y , c1.second + 0.5 * delta,
			                  c2.first.x, c2.first.y, c2.second + 0.5 * delta);

			double x_center, y_center, r1, r2;
			sampler_.getBound(x_center, y_center, r1, r2);

			sampler.setBound(x_center - r2, x_center + r2, y_center - r2, y_center + r2);
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 3;
		}
		
		std::size_t active_data_dim() const
		{
			return 3;
		}

		Scaler getScaler() const
		{
			double x_min, x_max, y_min, y_max;
			sampler.getBound(x_min, x_max, y_min, y_max);
			DataVector v_min(3);
			DataVector v_max(3);
			v_min[0] = x_min; v_min[1] = y_min; v_min[2] = -boost::math::constants::pi<double>();
			v_max[0] = x_max; v_max[1] = y_max; v_max[2] = boost::math::constants::pi<double>();

			return Scaler(v_min, v_max, 3);
		}
		
		SamplerSE2 sampler;
		
		Collider2D collider;
	};
	
	class ContactSpaceR3
	{
	public:
		typedef FLANN_WRAPPER::DistanceRN DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider3D ColliderType;

		ContactSpaceR3(C2A_Model* model1, C2A_Model* model2, double delta = 0) : collider(model1, model2)
		{
			AABB3D aabb1 = computeAABB(model1);
			AABB3D aabb2 = computeAABB(model2);
			
			double delta_x_min = aabb1.b_min[0] - aabb2.b_max[0] - delta;
			double delta_x_max = aabb1.b_max[0] - aabb2.b_min[0] + delta;
			
			double delta_y_min = aabb1.b_min[1] - aabb2.b_max[1] - delta;
			double delta_y_max = aabb1.b_max[1] - aabb2.b_min[1] + delta;
			
			double delta_z_min = aabb1.b_min[2] - aabb2.b_max[2] - delta;
			double delta_z_max = aabb1.b_max[2] - aabb2.b_min[2] + delta;
			
			sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max, delta_z_min, delta_z_max);
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			DataVector v(6);

			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v_ = sampler.sample();
			
				v[0] = v_[0];
				v[1] = v_[1];
				v[2] = v_[2];
				v[3] = 0;
				v[4] = 0;
				v[5] = 0;
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			DataVector v(6);

			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v_ = sampler.sample();
			
				v[0] = v_[0];
				v[1] = v_[1];
				v[2] = v_[2];
				samples.push_back(v);
			}

			return samples;
		}

		std::vector<ContactSpaceSampleData> uniform_sample(libm3d::model& P, libm3d::model& Q, double d, double shift) const
		{
			double R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
			std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

			std::vector<ContactSpaceSampleData> samples;
			DataVector v(6);
			for(std::size_t i = 0; i < samples_.size(); ++i)
			{
				v[0] = samples_[i][0];
				v[1] = samples_[i][1];
				v[2] = samples_[i][2];

				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}
			
			return samples;
		}

		std::vector<DataVector> uniform_sample0(libm3d::model& P, libm3d::model& Q, double d, double shift) const
		{
			double R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
			std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

			std::vector<DataVector> samples;
			DataVector v(6);
			for(std::size_t i = 0; i < samples_.size(); ++i)
			{
				v[0] = samples_[i][0];
				v[1] = samples_[i][1];
				v[2] = samples_[i][2];
				samples.push_back(v);
			}

			return samples;
		}
		
		std::vector<ContactSpaceSampleData> data;
		
		std::size_t data_dim() const
		{
			return 6;
		}
		
		std::size_t active_data_dim() const
		{
			return 3;
		}

		Scaler getScaler() const
		{
			double x_min, y_min, z_min, x_max, y_max, z_max;
			sampler.getBound(x_min, x_max, y_min, y_max, z_min, z_max);
			DataVector v_min(3), v_max(3);
			v_min[0] = x_min; v_min[1] = y_min; v_min[2] = z_min;
			v_max[0] = x_max; v_max[1] = y_max; v_max[2] = z_max;
			return Scaler(v_min, v_max, 3);
		}
		
		SamplerR3 sampler;
		
		Collider3D collider;
	};
	
	class ContactSpaceSE3Euler
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE3EulerAngle DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider3D ColliderType;

		ContactSpaceSE3Euler(C2A_Model* model1, C2A_Model* model2, double delta = 0) : collider(model1, model2)
		{
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			
			sampler.setBound(model1->com[0], model1->com[1], model1->com[2], model1->radius + 0.5 * delta,
			                 model2->com[0], model2->com[1], model2->com[2], model2->radius + 0.5 * delta);
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 6;
		}
		
		std::size_t active_data_dim() const
		{
			return 6;
		}

		Scaler getScaler() const
		{
			double center_x, center_y, center_z, r1, r2;
			sampler.getBound(center_x, center_y, center_z, r1, r2);
			DataVector v_min(6), v_max(6);
			v_min[0] = center_x - r2; v_min[1] = center_y - r2; v_min[2] = center_z - r2; v_min[3] = -boost::math::constants::pi<double>(); v_min[4] = -boost::math::constants::pi<double>(); v_min[5] = -boost::math::constants::pi<double>();
			v_max[0] = center_x + r2; v_max[1] = center_y + r2; v_max[2] = center_z + r2; v_max[3] = boost::math::constants::pi<double>(); v_max[4] = boost::math::constants::pi<double>(); v_max[5] = boost::math::constants::pi<double>();
			return Scaler(v_min, v_max, 6);
		}
		
		SamplerSE3Euler_ball sampler;
		
		Collider3D collider;
	};
	
	class ContactSpaceSE3Quat
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE3Quat DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider3D ColliderType;

		ContactSpaceSE3Quat(C2A_Model* model1, C2A_Model* model2, double delta = 0) : collider(model1, model2)
		{
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			
			sampler.setBound(model1->com[0], model1->com[1], model1->com[2], model1->radius + 0.5 * delta,
			                 model2->com[0], model2->com[1], model2->com[2], model2->radius + 0.5 * delta);
			                 
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			for(std::size_t i = 0; i < n; ++i)
			{
				DataVector v = sampler.sample();
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 7;
		}
		
		std::size_t active_data_dim() const
		{
			return 7;
		}

		Scaler getScaler() const
		{
			double center_x, center_y, center_z, r1, r2;
			sampler.getBound(center_x, center_y, center_z, r1, r2);
			DataVector v_min(7), v_max(7);
			v_min[0] = center_x - r2; v_min[1] = center_y - r2; v_min[2] = center_z - r2; v_min[3] = -1; v_min[4] = -1; v_min[5] = -1; v_min[6] = -1;
			v_max[0] = center_x + r2; v_max[1] = center_y + r2; v_max[2] = center_z + r2; v_max[3] = 1; v_max[4] = 1; v_max[5] = 1; v_max[6] = 1;
			return Scaler(v_min, v_max, 7);
		}
		
		SamplerSE3Quat_ball sampler;
		
		Collider3D collider;
	};
	
	
	class ContactSpaceSE3Euler2
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE3EulerAngle DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider3D ColliderType;

		ContactSpaceSE3Euler2(C2A_Model* model1, C2A_Model* model2, double delta_ = 0) : collider(model1, model2)
		{
			aabb1 = computeAABB(model1);
			aabb2 = computeAABB(model2);
			delta = delta_;
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;
			double r_[4];
			double R[3][3];
			DataVector v(6);
			double a, b, c;

			for(std::size_t i = 0; i < n; ++i)
			{
				rng.quaternion(r_);
				Quaternion r(r_[0], r_[1], r_[2], r_[3]);

				Quat2Rot(R, r);

				AABB3D new_aabb2 = rotate(aabb2, R);

				double delta_x_min = aabb1.b_min[0] - new_aabb2.b_max[0] - delta;
				double delta_x_max = aabb1.b_max[0] - new_aabb2.b_min[0] + delta;

				double delta_y_min = aabb1.b_min[1] - new_aabb2.b_max[1] - delta;
				double delta_y_max = aabb1.b_max[1] - new_aabb2.b_min[1] + delta;

				double delta_z_min = aabb1.b_min[2] - new_aabb2.b_max[2] - delta;
				double delta_z_max = aabb1.b_max[2] - new_aabb2.b_min[2] + delta;

				sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max, delta_z_min, delta_z_max);
				DataVector t = sampler.sample();

				Quat2Euler(a, b, c, r);
				v[0] = t[0];
				v[1] = t[1];
				v[2] = t[2];
				v[3] = a;
				v[4] = b;
				v[5] = c;

				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;
			double r_[4];
			double R[3][3];
			DataVector v(6);
			double a, b, c;

			for(std::size_t i = 0; i < n; ++i)
			{
				rng.quaternion(r_);
				Quaternion r(r_[0], r_[1], r_[2], r_[3]);

				Quat2Rot(R, r);

				AABB3D new_aabb2 = rotate(aabb2, R);

				double delta_x_min = aabb1.b_min[0] - new_aabb2.b_max[0] - delta;
				double delta_x_max = aabb1.b_max[0] - new_aabb2.b_min[0] + delta;

				double delta_y_min = aabb1.b_min[1] - new_aabb2.b_max[1] - delta;
				double delta_y_max = aabb1.b_max[1] - new_aabb2.b_min[1] + delta;

				double delta_z_min = aabb1.b_min[2] - new_aabb2.b_max[2] - delta;
				double delta_z_max = aabb1.b_max[2] - new_aabb2.b_min[2] + delta;

				sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max, delta_z_min, delta_z_max);
				DataVector t = sampler.sample();

				Quat2Euler(a, b, c, r);
				v[0] = t[0];
				v[1] = t[1];
				v[2] = t[2];
				v[3] = a;
				v[4] = b;
				v[5] = c;
				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 6;
		}
		
		std::size_t active_data_dim() const
		{
			return 6;
		}

		Scaler getScaler() const
		{
			double x_min, y_min, z_min, x_max, y_max, z_max;
			sampler.getBound(x_min, x_max, y_min, y_max, z_min, z_max);
			DataVector v_min(6), v_max(6);
			v_min[0] = x_min; v_min[1] = y_min; v_min[2] = z_min; v_min[3] = -boost::math::constants::pi<double>(); v_min[4] = -boost::math::constants::pi<double>(); v_min[5] = -boost::math::constants::pi<double>();
			v_max[0] = x_max; v_max[1] = y_max; v_max[2] = z_max; v_max[3] = boost::math::constants::pi<double>(); v_max[4] = boost::math::constants::pi<double>(); v_max[5] = boost::math::constants::pi<double>();
			return Scaler(v_min, v_max, 6);
		}
		
		mutable SamplerR3 sampler;
		
		mutable RNG rng;
		
		Collider3D collider;
		
		AABB3D aabb1, aabb2;
		
		double delta;
	};
	
	class ContactSpaceSE3Quat2
	{
	public:
		typedef FLANN_WRAPPER::DistanceSE3Quat DistanceType;
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistanceType;
		typedef Collider3D ColliderType;

		ContactSpaceSE3Quat2(C2A_Model* model1, C2A_Model* model2, double delta_ = 0) : collider(model1, model2)
		{
			aabb1 = computeAABB(model1);
			aabb2 = computeAABB(model2);
			delta = delta_;
		}
		
		std::vector<ContactSpaceSampleData> uniform_sample(std::size_t n) const
		{
			std::vector<ContactSpaceSampleData> samples;

			double r_[4];
			double R[3][3];
			DataVector v(7);

			for(std::size_t i = 0; i < n; ++i)
			{
				rng.quaternion(r_);

				Quaternion r(r_[0], r_[1], r_[2], r_[3]);

				Quat2Rot(R, r);

				AABB3D new_aabb2 = rotate(aabb2, R);

				double delta_x_min = aabb1.b_min[0] - new_aabb2.b_max[0] - delta;
				double delta_x_max = aabb1.b_max[0] - new_aabb2.b_min[0] + delta;

				double delta_y_min = aabb1.b_min[1] - new_aabb2.b_max[1] - delta;
				double delta_y_max = aabb1.b_max[1] - new_aabb2.b_min[1] + delta;

				double delta_z_min = aabb1.b_min[2] - new_aabb2.b_max[2] - delta;
				double delta_z_max = aabb1.b_max[2] - new_aabb2.b_min[2] + delta;

				sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max, delta_z_min, delta_z_max);
				DataVector t = sampler.sample();


				v[0] = t[0];
				v[1] = t[1];
				v[2] = t[2];
				v[3] = r_[0];
				v[4] = r_[1];
				v[5] = r_[2];
				v[6] = r_[3];

				bool col = collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			return samples;
		}

		std::vector<DataVector> uniform_sample0(std::size_t n) const
		{
			std::vector<DataVector> samples;

			double r_[4];
			double R[3][3];
			DataVector v(7);

			for(std::size_t i = 0; i < n; ++i)
			{
				rng.quaternion(r_);

				Quaternion r(r_[0], r_[1], r_[2], r_[3]);

				Quat2Rot(R, r);

				AABB3D new_aabb2 = rotate(aabb2, R);

				double delta_x_min = aabb1.b_min[0] - new_aabb2.b_max[0] - delta;
				double delta_x_max = aabb1.b_max[0] - new_aabb2.b_min[0] + delta;

				double delta_y_min = aabb1.b_min[1] - new_aabb2.b_max[1] - delta;
				double delta_y_max = aabb1.b_max[1] - new_aabb2.b_min[1] + delta;

				double delta_z_min = aabb1.b_min[2] - new_aabb2.b_max[2] - delta;
				double delta_z_max = aabb1.b_max[2] - new_aabb2.b_min[2] + delta;

				sampler.setBound(delta_x_min, delta_x_max, delta_y_min, delta_y_max, delta_z_min, delta_z_max);
				DataVector t = sampler.sample();


				v[0] = t[0];
				v[1] = t[1];
				v[2] = t[2];
				v[3] = r_[0];
				v[4] = r_[1];
				v[5] = r_[2];
				v[6] = r_[3];

				samples.push_back(v);
			}

			return samples;
		}
		
		std::size_t data_dim() const
		{
			return 7;
		}
		
		std::size_t active_data_dim() const
		{
			return 7;
		}

		Scaler getScaler() const
		{
			double x_min, y_min, z_min, x_max, y_max, z_max;
			sampler.getBound(x_min, x_max, y_min, y_max, z_min, z_max);
			DataVector v_min(7), v_max(7);
			v_min[0] = x_min; v_min[1] = y_min; v_min[2] = z_min; v_min[3] = -1; v_min[4] = -1; v_min[5] = -1; v_min[6] = -1;
			v_max[0] = x_max; v_max[1] = y_max; v_max[2] = z_max; v_min[3] = 1; v_min[4] = 1; v_min[5] = 1; v_min[6] = 1;
			return Scaler(v_min, v_max, 7);
		}
		
		mutable SamplerR3 sampler;
		
		mutable RNG rng;
		
		Collider3D collider;
		
		AABB3D aabb1, aabb2;
		
		double delta;
	};
	
	inline std::ofstream& asciiWriter(std::ofstream& os, const std::vector<ContactSpaceSampleData>& samples)
	{
		if(samples.size() == 0) return os;

		os << samples.size() << " " << samples[0].v.dim() << std::endl;
		std::size_t dim = samples[0].v.dim();
		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			os << samples[i].col << " ";
			for(std::size_t j = 0; j < dim; ++j)
			{
				os << samples[i].v[j] << " ";
			}
			os << std::endl;
		}
		
		return os;
	}
	
	inline std::ofstream& binaryWriter(std::ofstream& os, const std::vector<ContactSpaceSampleData>& samples)
	{
		if(samples.size() == 0) return os;

		std::size_t data_num = samples.size();
		std::size_t data_dim = samples[0].v.dim();
		os.write((char*)&data_num, sizeof(std::size_t));
		os.write((char*)&data_dim, sizeof(std::size_t));
		
		double* v = new double[data_dim];
		char col;
		for(std::size_t i = 0; i < data_num; ++i)
		{
			col = samples[i].col;
			for(std::size_t j = 0; j < data_dim; ++j)
			{
				v[j] = samples[i].v[j];
			}
			
			os.write(&col, sizeof(char));
			os.write((char*)v, sizeof(double) * data_dim);
		}
		
		delete [] v;
		
		return os;
	}
	

	inline std::ifstream& asciiReader(std::ifstream& is, std::vector<ContactSpaceSampleData>& samples)
	{
		std::size_t data_num, data_dim;
		is >> data_num >> data_dim;
		if(samples.size() != 0)
			assert(data_dim == samples[0].v.dim());
			
		DataVector v(data_dim);
		bool col;
		for(std::size_t i = 0; i < data_num; ++i)
		{
			is >> col;
			for(std::size_t j = 0; j < data_dim; ++j)
				is >> v[j];
				
			samples.push_back(ContactSpaceSampleData(v, col));
		}
		
		return is;
	}
	
	inline std::ifstream& binaryReader(std::ifstream& is, std::vector<ContactSpaceSampleData>& samples)
	{
		std::size_t data_num, data_dim;
		is.read((char*)&data_num, sizeof(std::size_t));
		is.read((char*)&data_dim, sizeof(std::size_t));
		if(samples.size() != 0)
			assert(data_dim == samples[0].v.dim());
			
		DataVector v(data_dim);
		bool col;
		std::size_t buffer_size = (sizeof(double) * data_dim + sizeof(char)) * data_num;
		char* buffer = new char[buffer_size];
		is.read(buffer, buffer_size);
		
		char* p = buffer;
		for(std::size_t i = 0; i < data_num; ++i)
		{
			col = (bool)(*p);
			p++;
			
			for(std::size_t j = 0; j < data_dim; ++j)
			{
				double t = 0;
				memcpy((void*)&t, (void*)p, sizeof(double));
				v[j] = t;
				p = p + sizeof(double);
			}
			samples.push_back(ContactSpaceSampleData(v, col));
		}
		
		delete [] buffer;
		
		return is;
	}


	template<typename Distance>
	void generateIndex(const std::vector<ContactSpaceSampleData>& data, std::size_t active_data_dim,
		flann::Index<Distance>*& index,
		const flann::IndexParams& params)
	{
		std::size_t num_data = data.size();
		std::size_t dim_data = active_data_dim;
		flann::Matrix<double> dataset = flann::Matrix<double>(new double[num_data * dim_data], num_data, dim_data);

		for(std::size_t i = 0; i < num_data; ++i)
		{
			for(std::size_t j = 0; j < dim_data; ++j)
			{
				dataset[i][j] = data[i].v[j];
			}
		}

		index = new flann::Index<Distance>(dataset, params);
		index->buildIndex();
	}

	template<typename Distance>
	void knnSearch(const std::vector<ContactSpaceSampleData>& queries, 
       std::size_t active_data_dim,
       flann::Index<Distance>* index,
       std::vector<std::vector<int> >& indices,
	   std::vector<std::vector<double> >& dists,
	   std::size_t knn,
	   const flann::SearchParams& params)
	{
		std::size_t num_queries = queries.size();
		std::size_t dim_data = active_data_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[num_queries * dim_data], num_queries, dim_data);

		for(std::size_t i = 0; i < num_queries; ++i)
		{
			for(std::size_t j = 0; j < dim_data; ++j)
			{
				queryset[i][j] = queries[i].v[j];
			}
		}

		index->knnSearch(queryset, indices, dists, knn, params);
	}

	template<typename Distance>
	void radiusSearch(const std::vector<ContactSpaceSampleData>& queries,
		std::size_t active_data_dim,
		flann::Index<Distance>* index,
	    std::vector<std::vector<int> >& indices,
	    std::vector<std::vector<double> >& dists,
		double radius,
		const flann::SearchParams& params)
	{
		std::size_t num_queries = queries.size();
		std::size_t dim_data = active_data_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[num_queries * dim_data], num_queries, dim_data);

		for(std::size_t i = 0; i < num_queries; ++i)
		{
			for(std::size_t j = 0; j < dim_data; ++j)
			{
				queryset[i][j] = queries[i].v[j];
			}
		}
		index->radiusSearch(queryset, indices, dists, radius, params);
	}

	template<typename Distance>
	void knnSearch(const ContactSpaceSampleData& query,
       std::size_t active_data_dim,
       flann::Index<Distance>* index,
       std::vector<int>& indices,
	   std::vector<double>& dists,
	   std::size_t knn,
	   const flann::SearchParams& params)
	{
		std::size_t dim_data = active_data_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim_data], 1, dim_data);


		for(std::size_t j = 0; j < dim_data; ++j)
		{
			queryset[0][j] = query.v[j];
		}

		std::vector<std::vector<int> > indices_;
		std::vector<std::vector<double> > dists_;

		index->knnSearch(queryset, indices_, dists_, knn, params);
		indices = indices_[0];
		dists = dists_[0];
	}

	template<typename Distance>
	void radiusSearch(const ContactSpaceSampleData& query,
		std::size_t active_data_dim,
		flann::Index<Distance>* index,
	    std::vector<int>& indices,
	    std::vector<double>& dists,
		double radius,
		const flann::SearchParams& params)
	{
		std::size_t dim_data = active_data_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim_data], 1, dim_data);

		for(std::size_t j = 0; j < dim_data; ++j)
		{
			queryset[0][j] = query.v[j];
		}

		std::vector<std::vector<int> > indices_;
		std::vector<std::vector<double> > dists_;

		index->radiusSearch(queryset, indices_, dists_, radius, params);
		indices = indices_[0];
		dists = dists_[0];
	}

	inline std::ofstream& operator << (std::ofstream& os, const std::vector<std::vector<int> >& indices)
	{
		for(std::size_t i = 0; i < indices.size(); ++i)
		{
			for(std::size_t j = 0; j < indices[i].size(); ++j)
			{
				os << indices[i][j] << " ";
			}
			os << std::endl;
		}

		return os;
	}

	inline std::ofstream& operator << (std::ofstream& os, const std::vector<std::vector<double> >& dists)
	{
		for(std::size_t i = 0; i < dists.size(); ++i)
		{
			for(std::size_t j = 0; j < dists[i].size(); ++j)
			{
				os << dists[i][j] << " ";
			}
			os << std::endl;
		}

		return os;
	}
	
}

#endif