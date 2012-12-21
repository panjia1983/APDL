#ifndef DATA_VECTOR_H
#define DATA_VECTOR_H

#include <vector>

namespace APDL
{


	class DataVector
	{
	public:
		DataVector(std::size_t n)
		{
			data.resize(n, 0);
		}

		DataVector() {}
		
		std::size_t dim() const
		{
			return data.size();
		}

		void setDim(std::size_t d)
		{
			data.resize(d);
		}
		
		void setData(const std::vector<double>& data_)
		{
			if(data.size() > data_.size())
			{
				std::size_t i = 0;
				for(i = 0; i < data_.size(); ++i)
					data[i] = data_[i];
				for(; i < data.size(); ++i)
					data[i] = 0;
			}
			else
			{
				for(std::size_t i = 0; i < data.size(); ++i)
					data[i] = data_[i];
			}
		}
		
		double operator [](std::size_t i) const
		{
			return data[i];
		}
		
		double& operator [](std::size_t i)
		{
			return data[i];
		}

		DataVector operator + (const DataVector& other) const
		{
			std::size_t dim = (std::min)(data.size(), other.dim());
			DataVector res(*this);
			for(std::size_t i = 0; i < dim; ++i)
			{
				res[i] += other[i];
			}
			return res;
		}

		DataVector operator - (const DataVector& other) const
		{
			std::size_t dim = (std::min)(data.size(), other.dim());
			DataVector res(*this);
			for(std::size_t i = 0; i < dim; ++i)
			{
				res[i] -= other[i];
			}
			return res;
		}

		DataVector operator - () const
		{
			DataVector res(*this);
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				res[i] = -res[i];
			}	
			return res;
		}

		DataVector operator * (double t) const
		{
			DataVector res(*this);
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				res[i] *= t;
			}
			return res;
		}

		DataVector& operator += (const DataVector& other)
		{
			std::size_t dim = (std::min)(data.size(), other.dim());
			for(std::size_t i = 0; i < dim; ++i)
			{
				data[i] += other[i];
			}
			return *this;
		}

		DataVector& operator -= (const DataVector& other)
		{
			std::size_t dim = (std::min)(data.size(), other.dim());
			for(std::size_t i = 0; i < dim; ++i)
			{
				data[i] -= other[i];
			}
			return *this;
		}

		DataVector& operator *= (double t)
		{
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				data[i] *= t;
			}
			return *this;
		}

		void setZero()
		{
			for(std::size_t i = 0; i < data.size(); ++i)
				data[i] = 0;
		}
		
	protected:
		std::vector<double> data;
	};

	inline double dot_prod(const DataVector& a, const DataVector& b)
	{
		double sum = 0;
		std::size_t dim = (std::min)(a.dim(), b.dim());
		for(std::size_t i = 0; i < dim; ++i)
		{
			sum += a[i] * b[i];
		}
		return sum;
	}

	inline DataVector operator * (double t, const DataVector& v)
	{
		DataVector res(v);
		for(std::size_t i = 0; i < res.dim(); ++i)
		{
			res[i] *= t;
		}
		return res;
	}
	
}

#endif