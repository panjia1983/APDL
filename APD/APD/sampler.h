#ifndef SAMPLER_H
#define SAMPLER_H


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "math_utility.h"

namespace APDL
{


	/** \brief Random number generation. An instance of this class
	    cannot be used by multiple threads at once (member functions
	    are not const). However, the constructor is thread safe and
	    different instances can be used safely in any number of
	    threads. It is also guaranteed that all created instances will
	    have a different random seed. */
	class RNG
	{
	public:
	
		/** \brief Constructor. Always sets a different random seed */
		RNG(void);
		
		/** \brief Generate a random real between 0 and 1 */
		double uniform01(void)
		{
			return uni_();
		}
		
		/** \brief Generate a random real within given bounds: [\e lower_bound, \e upper_bound) */
		double uniformReal(double lower_bound, double upper_bound)
		{
			assert(lower_bound <= upper_bound);
			return (upper_bound - lower_bound) * uni_() + lower_bound;
		}
		
		/** \brief Generate a random integer within given bounds: [\e lower_bound, \e upper_bound] */
		int uniformInt(int lower_bound, int upper_bound)
		{
			int r = (int)floor(uniformReal((double)lower_bound, (double)(upper_bound) + 1.0));
			return (r > upper_bound) ? upper_bound : r;
		}
		
		/** \brief Generate a random boolean */
		bool uniformBool(void)
		{
			return uni_() <= 0.5;
		}
		
		/** \brief Generate a random real using a normal distribution with mean 0 and variance 1 */
		double gaussian01(void)
		{
			return normal_();
		}
		
		/** \brief Generate a random real using a normal distribution with given mean and variance */
		double gaussian(double mean, double stddev)
		{
			return normal_() * stddev + mean;
		}
		
		/** \brief Generate a random real using a half-normal distribution. The value is within specified bounds [\e
		    r_min, \e r_max], but with a bias towards \e r_max. The function is implemended using a Gaussian distribution with
		    mean at \e r_max - \e r_min. The distribution is 'folded' around \e r_max axis towards \e r_min.
		    The variance of the distribution is (\e r_max - \e r_min) / \e focus. The higher the focus,
		    the more probable it is that generated numbers are close to \e r_max. */
		double halfNormalReal(double r_min, double r_max, double focus = 3.0);
		
		/** \brief Generate a random integer using a half-normal
		    distribution. The value is within specified bounds ([\e r_min, \e r_max]), but
		    with a bias towards \e r_max. The function is implemented on top of halfNormalReal() */
		int    halfNormalInt(int r_min, int r_max, double focus = 3.0);
		
		/** \brief Uniform random unit quaternion sampling. The computed value has the order (x,y,z,w) */
		void   quaternion(double value[4]);
		
		/** \brief Uniform random sampling of Euler roll-pitch-yaw angles, each in the range [-pi, pi). The computed value has the order (roll, pitch, yaw) */
		void   eulerRPY(double value[3]);
		
		/** \brief Uniform random sample on a disk with radius from r_min to r_max */
		void disk(double r_min, double r_max, double& x, double& y);
		
		/** \brief Uniform random sample in a ball with radius from r_min to r_max */
		void ball(double r_min, double r_max, double& x, double& y, double& z);
		
		/** \brief Set the seed for random number generation. Use this
		    function to ensure the same sequence of random numbers is
		    generated. */
		static void setSeed(boost::uint32_t seed);
		
		/** \brief Get the seed used for random number
		    generation. Passing the returned value to setSeed() at a
		    subsequent execution of the code will ensure deterministic
		    (repeatable) behaviour. Useful for debugging. */
		static boost::uint32_t getSeed(void);
		
	private:
	
		boost::mt19937                                                           generator_;
		boost::uniform_real<>                                                    uniDist_;
		boost::normal_distribution<>                                             normalDist_;
		boost::variate_generator<boost::mt19937&, boost::uniform_real<> >        uni_;
		boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > normal_;
		
	};
	
	class SamplerR2
	{
	public:
		SamplerR2() {}
		
		SamplerR2(double x_min_, double x_max_, double y_min_, double y_max_)
		{
			setBound(x_min_, x_max_, y_min_, y_max_);
		}
		
		void setBound(double x_min_, double x_max_, double y_min_, double y_max_)
		{
			x_min = x_min_;
			x_max = x_max_;
			y_min = y_min_;
			y_max = y_max_;
		}

		void getBound(double& x_min_, double& x_max_, double& y_min_, double& y_max_) const
		{
			x_min_ = x_min;
			y_min_ = y_min;
			x_max_ = x_max;
			y_max_ = y_max;
		}
		
		DataVector sample() const
		{
			DataVector res(2);
			res[0] = rng.uniformReal(x_min, x_max);
			res[1] = rng.uniformReal(y_min, y_max);
			
			return res;
		}
		

		mutable RNG rng;

	protected:
		
		double x_min, x_max, y_min, y_max;
	};
	
	
	class SamplerSE2
	{
	public:
		SamplerSE2() {}
		
		SamplerSE2(double x_min_, double x_max_, double y_min_, double y_max_)
		{
			setBound(x_min_, x_max_, y_min_, y_max_);
		}
		
		void setBound(double x_min_, double x_max_, double y_min_, double y_max_)
		{
			x_min = x_min_;
			x_max = x_max_;
			y_min = y_min_;
			y_max = y_max_;
		}

		void getBound(double& x_min_, double& x_max_, double& y_min_, double& y_max_) const
		{
			x_min_ = x_min;
			y_min_ = y_min;
			x_max_ = x_max;
			y_max_ = y_max;
		}
		
		DataVector sample() const
		{
			DataVector res(3);
			res[0] = rng.uniformReal(x_min, x_max);
			res[1] = rng.uniformReal(y_min, y_max);
			res[2] = rng.uniformReal(-boost::math::constants::pi<double>(), boost::math::constants::pi<double>());
			
			return res;
		}
		
		mutable RNG rng;

	protected:
		
		double x_min, x_max, y_min, y_max;
	};
	
	class SamplerSE2_disk
	{
	public:
		SamplerSE2_disk() {}
		SamplerSE2_disk(double cx1_, double cy1_, double r1_,
		                double cx2_, double cy2_, double r2_)
		{
			setBound(cx1_, cy1_, r1_, cx2_, cy2_, r2_);
		}
		
		void setBound(double cx1_, double cy1_, double r1_,
		              double cx2_, double cy2_, double r2_)
		{
			c[0] = cx1_;
			c[1] = cy1_;
			double c2_len = std::sqrt(cx2_ * cx2_ + cy2_ * cy2_);
			r1 = std::max<double>(c2_len - (r1_ + r2_), 0);
			r2 = c2_len + r1_ + r2_;
		}

		void getBound(double& c1_, double& c2_, double& r1_, double& r2_) const
		{
			c1_ = c[0];
			c2_ = c[1];
			r1_ = r1;
			r2_ = r2;
		}
		
		DataVector sample() const
		{
			DataVector res(3);
			double x, y;
			rng.disk(r1, r2, x, y);
			res[0] = x + c[0];
			res[1] = y + c[1];
			res[2] = rng.uniformReal(-boost::math::constants::pi<double>(), boost::math::constants::pi<double>());
			return res;
		}
		
		mutable RNG rng;

	protected:
		
		double c[2];
		double r1, r2;
	};
	
	class SamplerR3
	{
	public:
		SamplerR3() {}
		
		SamplerR3(double x_min_, double x_max_,
		          double y_min_, double y_max_,
		          double z_min_, double z_max_)
		{
			setBound(x_min_, x_max_, y_min_, y_max_, z_min_, z_max_);
		}
		
		void setBound(double x_min_, double x_max_, double y_min_, double y_max_, double z_min_, double z_max_)
		{
			x_min = x_min_;
			x_max = x_max_;
			y_min = y_min_;
			y_max = y_max_;
			z_min = z_min_;
			z_max = z_max_;
		}

		void getBound(double& x_min_, double& x_max_, double& y_min_, double& y_max_, double& z_min_, double& z_max_) const
		{
			x_min_ = x_min;
			y_min_ = y_min;
			x_max_ = x_max;
			y_max_ = y_max;
			z_min_ = z_min;
			z_max_ = z_max;
		}
		
		DataVector sample() const
		{
			DataVector res(3);
			res[0] = rng.uniformReal(x_min, x_max);
			res[1] = rng.uniformReal(y_min, y_max);
			res[2] = rng.uniformReal(z_min, z_max);
			
			return res;
		}
		
		mutable RNG rng;

	protected:
		
		double x_min, x_max, y_min, y_max, z_min, z_max;
	};
	
	class SamplerSE3Euler
	{
	public:
		SamplerSE3Euler() {}
		
		SamplerSE3Euler(double x_min_, double x_max_,
		                double y_min_, double y_max_,
		                double z_min_, double z_max_)
		{
			setBound(x_min_, x_max_, y_min_, y_max_, z_min_, z_max_);
		}
		
		void setBound(double x_min_, double x_max_, double y_min_, double y_max_, double z_min_, double z_max_)
		{
			x_min = x_min_;
			x_max = x_max_;
			y_min = y_min_;
			y_max = y_max_;
			z_min = z_min_;
			z_max = z_max_;
		}

		void getBound(double& x_min_, double& x_max_, double& y_min_, double& y_max_, double& z_min_, double& z_max_) const
		{
			x_min_ = x_min;
			y_min_ = y_min;
			x_max_ = x_max;
			y_max_ = y_max;
			z_min_ = z_min;
			z_max_ = z_max;
		}
		
		DataVector sample() const
		{
			DataVector res(7);
			res[0] = rng.uniformReal(x_min, x_max);
			res[1] = rng.uniformReal(y_min, y_max);
			res[2] = rng.uniformReal(z_min, z_max);
			
			double quat[4];
			rng.quaternion(quat);
			
			res[3] = quat[0];
			res[4] = quat[1];
			res[5] = quat[2];
			res[6] = quat[3];
			
			return ConfigQuat2Euler(res);
		}
		
		mutable RNG rng;

	protected:
		
		double x_min, x_max, y_min, y_max, z_min, z_max;
	};
	
	
	class SamplerSE3Quat
	{
	public:
		SamplerSE3Quat() {}
		
		SamplerSE3Quat(double x_min_, double x_max_,
		               double y_min_, double y_max_,
		               double z_min_, double z_max_)
		{
			setBound(x_min_, x_max_, y_min_, y_max_, z_min_, z_max_);
		}
		
		void setBound(double x_min_, double x_max_, double y_min_, double y_max_, double z_min_, double z_max_)
		{
			x_min = x_min_;
			x_max = x_max_;
			y_min = y_min_;
			y_max = y_max_;
			z_min = z_min_;
			z_max = z_max_;
		}

		void getBound(double& x_min_, double& x_max_, double& y_min_, double& y_max_, double& z_min_, double& z_max_) const
		{
			x_min_ = x_min;
			y_min_ = y_min;
			x_max_ = x_max;
			y_max_ = y_max;
			z_min_ = z_min;
			z_max_ = z_max;
		}
		
		DataVector sample() const
		{
			DataVector res(6);
			res[0] = rng.uniformReal(x_min, x_max);
			res[1] = rng.uniformReal(y_min, y_max);
			res[2] = rng.uniformReal(z_min, z_max);
			
			double quat[4];
			rng.quaternion(quat);
			res[3] = quat[0];
			res[4] = quat[1];
			res[5] = quat[2];
			res[6] = quat[3];
			
			return res;
		}
		
		mutable RNG rng;

	protected:
		
		double x_min, x_max, y_min, y_max, z_min, z_max;
	};
	
	
	
	class SamplerSE3Euler_ball
	{
	public:
		SamplerSE3Euler_ball() {}
		
		SamplerSE3Euler_ball(double cx1_, double cy1_, double cz1_, double r1_,
		                     double cx2_, double cy2_, double cz2_, double r2_)
		{
			setBound(cx1_, cy1_, cz1_, r1_, cx2_, cy2_, cz2_, r2_);
		}
		
		void setBound(double cx1_, double cy1_, double cz1_, double r1_,
		              double cx2_, double cy2_, double cz2_, double r2_)
		{
			c[0] = cx1_;
			c[1] = cy1_;
			c[2] = cz1_;
			double c2_len = std::sqrt(cx2_ * cx2_ + cy2_ * cy2_ + cz2_ * cz2_);
			r1 = std::max<double>(c2_len - (r1_ + r2_), 0);
			r2 = c2_len + r1_ + r2_;
		}

		void getBound(double& c1_, double& c2_, double& c3_, double& r1_, double& r2_) const
		{
			c1_ = c[0];
			c2_ = c[1];
			c3_ = c[2];
			r1_ = r1;
			r2_ = r2;
		}
		
		DataVector sample() const
		{
			DataVector res(7);
			double x, y, z;
			rng.ball(r1, r2, x, y, z);
			res[0] = x;
			res[1] = y;
			res[2] = z;
			
			double quat[4];
			rng.quaternion(quat);
			
			res[3] = quat[0];
			res[4] = quat[1];
			res[5] = quat[2];
			res[6] = quat[3];
			
			return ConfigQuat2Euler(res);
		}
		
		mutable RNG rng;

	protected:
		
		double c[3];
		double r1, r2;
	};
	
	
	class SamplerSE3Quat_ball
	{
	public:
		SamplerSE3Quat_ball() {}
		
		SamplerSE3Quat_ball(double cx1_, double cy1_, double cz1_, double r1_,
		                    double cx2_, double cy2_, double cz2_, double r2_)
		{
			setBound(cx1_, cy1_, cz1_, r1_, cx2_, cy2_, cz2_, r2_);
		}
		
		void setBound(double cx1_, double cy1_, double cz1_, double r1_,
		              double cx2_, double cy2_, double cz2_, double r2_)
		{
			c[0] = cx1_;
			c[1] = cy1_;
			c[2] = cz1_;
			double c2_len = std::sqrt(cx2_ * cx2_ + cy2_ * cy2_ + cz2_ * cz2_);
			r1 = std::max<double>(c2_len - (r1_ + r2_), 0);
			r2 = c2_len + r1_ + r2_;
		}

		void getBound(double& c1_, double& c2_, double& c3_, double& r1_, double& r2_) const
		{
			c1_ = c[0];
			c2_ = c[1];
			c3_ = c[2];
			r1_ = r1;
			r2_ = r2;
		}
		
		DataVector sample() const
		{
			DataVector res(7);
			double x, y, z;
			rng.ball(r1, r2, x, y, z);
			res[0] = x;
			res[1] = y;
			res[2] = z;
			
			double quat[4];
			rng.quaternion(quat);
			res[3] = quat[0];
			res[4] = quat[1];
			res[5] = quat[2];
			res[6] = quat[3];
			
			return res;
		}
		
		mutable RNG rng;

	protected:
		
		double c[3];
		double r1, r2;
	};
	
	
}


#endif