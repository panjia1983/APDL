#include "contact_space.h"

namespace APDL
{

std::vector<ContactSpaceSampleData> uniform_sample_Euler(
	libm3d::model& P, libm3d::model& Q, double d, double shift, std::size_t nrot,
	const Collider3D& collider,
	RNG& rng)
{
	std::vector<ContactSpaceSampleData> samples;
	double r_[4];
	double R[3][3];
	DataVector v(6);
	double a, b, c;

	for(std::size_t i = 0; i < nrot; ++i)
	{
		rng.quaternion(r_);
		Quaternion r(r_[0], r_[1], r_[2], r_[3]);
		Quat2Rot(R, r);
		Quat2Euler(a, b, c, r);

		std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

		for(std::size_t j = 0; j < samples_.size(); ++j)
		{
			v[0] = samples_[j][0];
			v[1] = samples_[j][1];
			v[2] = samples_[j][2];
			v[3] = a;
			v[4] = b;
			v[5] = c;

			bool col = collider.isCollide(v);
			samples.push_back(ContactSpaceSampleData(v, col));
		}
	}

	return samples;
}

std::vector<DataVector> uniform_sample0_Euler(
	libm3d::model& P, libm3d::model& Q, double d, double shift, std::size_t nrot,
	RNG& rng)
{
	std::vector<DataVector> samples;
	double r_[4];
	double R[3][3];
	DataVector v(6);
	double a, b, c;

	for(std::size_t i = 0; i < nrot; ++i)
	{
		rng.quaternion(r_);
		Quaternion r(r_[0], r_[1], r_[2], r_[3]);
		Quat2Rot(R, r);
		Quat2Euler(a, b, c, r);

		std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

		for(std::size_t j = 0; j < samples_.size(); ++j)
		{
			v[0] = samples_[j][0];
			v[1] = samples_[j][1];
			v[2] = samples_[j][2];
			v[3] = a;
			v[4] = b;
			v[5] = c;

			samples.push_back(v);
		}
	}

	return samples;
}

std::vector<ContactSpaceSampleData> uniform_sample_Quat(
	libm3d::model& P, libm3d::model& Q, double d, double shift, std::size_t nrot,
	const Collider3D& collider,
	RNG& rng)
{
	std::vector<ContactSpaceSampleData> samples;
	double r_[4];
	double R[3][3];
	DataVector v(7);

	for(std::size_t i = 0; i < nrot; ++i)
	{
		rng.quaternion(r_);
		Quaternion r(r_[0], r_[1], r_[2], r_[3]);
		Quat2Rot(R, r);

		std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

		for(std::size_t j = 0; j < samples_.size(); ++j)
		{
			v[0] = samples_[j][0];
			v[1] = samples_[j][1];
			v[2] = samples_[j][2];
			v[3] = r_[0];
			v[4] = r_[1];
			v[5] = r_[2];
			v[6] = r_[3];

			bool col = collider.isCollide(v);
			samples.push_back(ContactSpaceSampleData(v, col));
		}
	}

	return samples;
}

std::vector<DataVector> uniform_sample0_Quat(
	libm3d::model& P, libm3d::model& Q, double d, double shift, std::size_t nrot,
	RNG& rng)
{
	std::vector<DataVector> samples;
	double r_[4];
	double R[3][3];
	DataVector v(6);

	for(std::size_t i = 0; i < nrot; ++i)
	{
		rng.quaternion(r_);
		Quaternion r(r_[0], r_[1], r_[2], r_[3]);
		Quat2Rot(R, r);

		std::vector<DataVector> samples_ = computePointMKDiff(P, Q, R, d, shift);

		for(std::size_t j = 0; j < samples_.size(); ++j)
		{
			v[0] = samples_[j][0];
			v[1] = samples_[j][1];
			v[2] = samples_[j][2];
			v[3] = r_[0];
			v[4] = r_[1];
			v[5] = r_[2];
			v[6] = r_[3];

			samples.push_back(v);
		}
	}

	return samples;
}

}