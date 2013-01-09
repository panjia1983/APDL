#ifndef ACTIVE_LEARNING_H
#define ACTIVE_LEARNING_H

#include "contact_space_learning.h"
#include <APD/profile.h>

extern std::ofstream test_file;

namespace APDL
{
	template<typename Learner>
	double empiricalErrorRatio(std::vector<ContactSpaceSampleData>& data, const Learner& learner)
	{
		DataVector v(learner.feature_dim);
		std::size_t n_error = 0;
		for(std::size_t i = 0; i < data.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.feature_dim; ++j)
				v[j] = data[i].v[j];
			if(data[i].col != (learner.predict(v)).label)
				n_error++;
		}

		return n_error / (double)data.size();
	}

	template<typename ContactSpace>
	double errorRatioOnGrid(const ContactSpace& contactspace, const SVMLearner& learner, std::size_t grid_n, bool input_scaled = false)
	{
		std::size_t dim = learner.feature_dim;
		DataVector upper(dim);
		DataVector lower(dim);

		if(learner.scaler && learner.use_scaler & input_scaled)
		{
			for(std::size_t i = 0; i < dim; ++i)
			{
				upper[i] = 1;
				lower[i] = 0;
			}
		}
		else
		{
			upper = contactspace.getScaler().v_max;
			lower = contactspace.getScaler().v_min;
		}


		std::size_t n = 1;
		for(std::size_t i = 0; i < dim; ++i)
			n *= grid_n;

		std::size_t* ids = new std::size_t[dim];
		DataVector v_col(contactspace.data_dim()); // cache for collision
		DataVector v_pred(dim); // cache for predict
		std::size_t n_error = 0;
		for(std::size_t i = 0; i < n; ++i)
		{
			std::size_t i_ = i;
			for(std::size_t j = 0; j < dim; ++j)
			{
				ids[j] = i_ % grid_n;
				i_ = i_ / grid_n;
			}

			for(std::size_t j = 0; j < dim; ++j)
			{
				v_col[j] = (upper[j] - lower[j]) * ids[j] / (double)grid_n + lower[j];
				v_pred[j] = v_col[j];
			}

			bool exact_col = contactspace.collider.isCollide(v_col);
			bool approx_col = ((learner.predict(v_pred)).label > 0);
			if(exact_col != approx_col)
				n_error++;
		}


		delete [] ids;

		return n_error / (double)n;
	}

	template<typename ContactSpace>
	double errorRatioOnGrid(const ContactSpace& contactspace, const MulticonlitronLearner& learner, std::size_t grid_n)
	{
		const Scaler& scaler = contactspace.getScaler();

		std::size_t dim = contactspace.active_data_dim();

		DataVector upper = scaler.v_max;
		DataVector lower = scaler.v_min;

		std::size_t n = 1;
		for(std::size_t i = 0; i < dim; ++i)
			n *= grid_n;

		std::size_t* ids = new std::size_t[dim];
		DataVector v_col(contactspace.data_dim()); // cache for collision
		DataVector v_pred(dim); // cache for predict
		std::size_t n_error = 0;
		for(std::size_t i = 0; i < n; ++i)
		{
			std::size_t i_ = i;
			for(std::size_t j = 0; j < dim; ++j)
			{
				ids[j] = i_ % grid_n;
				i_ = i_ / grid_n;
			}

			for(std::size_t j = 0; j < dim; ++j)
			{
				v_col[j] = (upper[j] - lower[j]) * ids[j] / (double)grid_n + lower[j];
				v_pred[j] = v_col[j];
			}

			for(std::size_t j = dim; j < contactspace.data_dim(); ++j)
				v_col[j] = 0;

			bool exact_col = contactspace.collider.isCollide(v_col);
			bool approx_col = (learner.predict(v_pred).label == 1);
			if(exact_col != approx_col)
				n_error++;
		}

		delete [] ids;

		return n_error / (double)n;
	}

	struct ActiveLearningParam
	{
		std::size_t num_init_uniform_samples;
		std::size_t num_exploit_samples_per_iter;
		std::size_t num_explore_samples_per_iter;
		std::size_t num_iter;

		std::size_t num_grid; //for test quality

		bool debug;

		ActiveLearningParam(std::size_t num_init_uniform_samples_,
			std::size_t num_exploit_samples_per_iter_,
			std::size_t num_explore_samples_per_iter_,
			std::size_t num_iter_) 
		{
			num_init_uniform_samples = num_init_uniform_samples_;
			num_exploit_samples_per_iter = num_exploit_samples_per_iter_;
			num_explore_samples_per_iter = num_explore_samples_per_iter_;
			num_iter = num_iter_;

			num_grid = 100;

			debug = false;
			debug_os = NULL;
		}

		std::string model_name;

		std::ofstream* debug_os;
	};


	template<typename ContactSpace, typename Learner, typename DecisionBoundarySampler>
	void active_learning(const ContactSpace& contactspace, Learner& learner, const DecisionBoundarySampler& decision_boundary_sampler, const ActiveLearningParam& param)
	{
		std::size_t num_uniform_samples = param.num_init_uniform_samples;
		std::size_t num_exploit = param.num_exploit_samples_per_iter;
		std::size_t num_explore = param.num_explore_samples_per_iter;

		tools::Profiler::Begin(param.model_name + "_init");
		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(num_uniform_samples); // initial_sample
		learner.learn(samples, contactspace.active_data_dim());
		tools::Profiler::End(param.model_name + "_init");

		if(param.debug && param.debug_os)
		{
			*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
			std::string model_file_name = param.model_name + "_base.txt";
			learner.save(model_file_name);
		}


		for(std::size_t i = 0; i < param.num_iter; ++i)
		{
			tools::Profiler::Begin(param.model_name + "_iteration");
			std::vector<DataVector> exploitation_samples = decision_boundary_sampler.sample(num_exploit);
			std::vector<ContactSpaceSampleData> exploration_samples = contactspace.uniform_sample(num_explore);

			int n_exploitation_error = 0;
			int n_exploration_error = 0;

			for(std::size_t j = 0; j < exploitation_samples.size(); ++j)
			{
				const DataVector& unscaled_exploitation_sample = (learner.scaler && learner.use_scaler) ? learner.scaler->unscale(exploitation_samples[j]) : exploitation_samples[j];

				DataVector v(contactspace.data_dim());
				for(std::size_t k = 0; k < contactspace.active_data_dim(); ++k) v[k] = unscaled_exploitation_sample[k];
				bool col = contactspace.collider.isCollide(v);
				samples.push_back(ContactSpaceSampleData(v, col));
			}

			for(std::size_t j = 0; j < exploration_samples.size(); ++j)
				samples.push_back(exploration_samples[j]);

			learner.learn(samples, contactspace.active_data_dim());

			tools::Profiler::End(param.model_name + "_iteration");

			SVMLearner* plearner = dynamic_cast<SVMLearner*>(&learner);
			if(plearner) 
			{
				plearner->setC(plearner->param.C + 1);
				//plearner->setGamma(plearner->param.gamma + 2);
			}

			if(param.debug && param.debug_os)
			{
				*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = param.model_name + "_" + convert.str() + ".txt";
				learner.save(model_file_name);
			}
		}

		if(param.debug && param.debug_os)
			*(param.debug_os) << std::endl;
	}

	template<typename ContactSpace, typename Learner, typename DecisionBoundarySampler>
	void active_learning2(const ContactSpace& contactspace, Learner& learner, const DecisionBoundarySampler& decision_boundary_sampler, const ActiveLearningParam& param)
	{
		std::size_t num_uniform_samples = param.num_init_uniform_samples;
		std::size_t num_exploit = param.num_exploit_samples_per_iter;
		std::size_t num_explore = param.num_explore_samples_per_iter;
		tools::Profiler::Begin(param.model_name + "_init");
		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(num_uniform_samples); // initial_sample
		learner.learn(samples, contactspace.active_data_dim());
		tools::Profiler::End(param.model_name + "_init");

		if(param.debug && param.debug_os)
		{
			*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
			std::string model_file_name = param.model_name + "_base.txt";
			learner.save(model_file_name);
		}


		double p = 0.5;
		double p_epsilon = 0.1;
		double p_lambda = 0.5;
		RNG rng;
		for(std::size_t i = 0; i < param.num_iter; ++i)
		{
			tools::Profiler::Begin(param.model_name + "_iteration");
			double ps = rng.uniform01();
			if(ps < p)
			{
				std::vector<ContactSpaceSampleData> exploration_samples = contactspace.uniform_sample(num_explore);
				for(std::size_t j = 0; j < exploration_samples.size(); ++j)
					samples.push_back(exploration_samples[j]);

				std::size_t n_error = 0;
				for(std::size_t j = 0; j < exploration_samples.size(); ++j)
				{
					if(exploration_samples[j].col != (learner.predict(exploration_samples[j].v).label > 0))
						n_error++;
				}

				double exploration_effectiveness = n_error / (double)exploration_samples.size();
				exploration_effectiveness = 2 * exploration_effectiveness - 1;

				p = (std::max)((std::min)(p * p_lambda * exp(exploration_effectiveness), 1 - p_epsilon), p_epsilon);
			}
			else
			{
				std::vector<DataVector> exploitation_samples = decision_boundary_sampler.sample(num_exploit);

				for(std::size_t j = 0; j < exploitation_samples.size(); ++j)
				{
					const DataVector& unscaled_exploitation_sample = (learner.scaler && learner.use_scaler) ? learner.scaler->unscale(exploitation_samples[j]) : exploitation_samples[j];

					DataVector v(contactspace.data_dim());
					for(std::size_t k = 0; k < contactspace.active_data_dim(); ++k) v[k] = unscaled_exploitation_sample[k];
					bool col = contactspace.collider.isCollide(v);
					samples.push_back(ContactSpaceSampleData(v, col));
				}
			}

			learner.learn(samples, contactspace.active_data_dim());
			tools::Profiler::End(param.model_name + "_iteration");

			SVMLearner* plearner = dynamic_cast<SVMLearner*>(&learner);
			if(plearner) 
			{
				// plearner->setC(plearner->param.C + 5);
				// plearner->setGamma(plearner->param.gamma + 1);
			}

			if(param.debug && param.debug_os)
			{
				*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = param.model_name + "_" + convert.str() + ".txt";
				learner.save(model_file_name);
			}
		}

		if(param.debug && param.debug_os)
			*(param.debug_os) << std::endl;
	}


	template<typename ContactSpace, typename Learner, typename DecisionBoundarySampler>
	void active_learning_incremental(const ContactSpace& contactspace, Learner& learner, const DecisionBoundarySampler& decision_boundary_sampler, const ActiveLearningParam& param)
	{
		std::size_t num_uniform_samples = param.num_init_uniform_samples;
		std::size_t num_exploit = param.num_exploit_samples_per_iter;
		std::size_t num_explore = param.num_explore_samples_per_iter;
		tools::Profiler::Begin(param.model_name + "_init");
		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(num_uniform_samples); // initial_sample
		learner.learn(samples, contactspace.active_data_dim());
		tools::Profiler::End(param.model_name + "_init");

		if(param.debug && param.debug_os)
		{
			*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
			std::string model_file_name = param.model_name + "_base.txt";
			learner.save(model_file_name);
		}

		for(std::size_t i = 0; i < param.num_iter; ++i)
		{
			tools::Profiler::Begin(param.model_name + "_iteration");
			samples.clear();
			tools::Profiler::Begin("decision sampler");
			std::vector<DataVector> exploitation_samples = decision_boundary_sampler.sample(num_exploit);
			tools::Profiler::End("decision sampler");
			std::vector<ContactSpaceSampleData> exploration_samples = contactspace.uniform_sample(num_explore);
			for(std::size_t j = 0; j < exploitation_samples.size(); ++j)
			{
				// unscale, because these samples are generated by learner and therefore is scaled.
				const DataVector& unscaled_exploitation_sample = (learner.scaler && learner.use_scaler) ? learner.scaler->unscale(exploitation_samples[j]) : exploitation_samples[j];

				DataVector v(contactspace.data_dim());
				for(std::size_t k = 0; k < contactspace.active_data_dim(); ++k) v[k] = unscaled_exploitation_sample[k];
				bool col = contactspace.collider.isCollide(v);

				samples.push_back(ContactSpaceSampleData(v, col));
			}


			for(std::size_t j = 0; j < exploration_samples.size(); ++j)
				samples.push_back(exploration_samples[j]);


			learner.incremental_learn(samples, contactspace.active_data_dim());

			tools::Profiler::End(param.model_name + "_iteration");

			if(param.debug && param.debug_os)
			{
				*(param.debug_os) << samples.size() << " " << empiricalErrorRatio(samples, learner) << " " << errorRatioOnGrid(contactspace, learner, param.num_grid) << " ";
				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = param.model_name + "_" + convert.str() + ".txt";
				learner.save(model_file_name);
			}
		}

		if(param.debug && param.debug_os)
			*(param.debug_os) << std::endl;
	}
	
}

#endif