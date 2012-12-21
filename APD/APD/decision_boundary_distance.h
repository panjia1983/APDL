#ifndef DECISION_BOUNDARY_DISTANCE_H
#define DECISION_BOUNDARY_DISTANCE_H

#include "learning/svm.h"
#include "distance_proxy.h"
#include "contact_space_learning.h"

#include <vector>

namespace APDL
{
	void sample_decision_boundary_brute_force(const struct svm_model* model, 
											  double* upper, double* lower,
											  std::size_t dim,
											  double on_decision_threshold, 
											  std::size_t subdivision_level, 
											  std::vector<std::vector<double> >& boundary);

	void sample_decision_boundary_brute_force(const struct svm_model* model,
											  double* upper, double* lower,
											  std::size_t dim,
											  double on_decision_threshold, 
											  std::size_t subdivision_level, 
											  std::vector<std::vector<double> >& boundary,
											  flann::Index<FLANN_WRAPPER::DistanceRN>*& index);

	double dist_to_decision_boundary_brute_force(const std::vector<std::vector<double> >& boundary, const struct svm_node* x, struct svm_node* x_res = NULL);

	double dist_to_decision_boundary_brute_force(const std::vector<std::vector<double> >& boundary, flann::Index<FLANN_WRAPPER::DistanceRN>* index, const struct svm_node* x, struct svm_node* x_res = NULL);

	struct SVMDistanceToDecisionBoundary_RoughLowerBound
	{
		SVMDistanceToDecisionBoundary_RoughLowerBound(const SVMLearner& learner_) : learner(learner_)
		{
			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;

			inv_hyperw = 1 / sqrt(learner.hyperw_normsqr);
		}

		~SVMDistanceToDecisionBoundary_RoughLowerBound()
		{
			delete [] node;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& cloest_v) const;


		const SVMLearner& learner;

		double inv_hyperw;

		mutable svm_node* node;
	};

	struct SVMDistanceToDecisionBoundary_Projection
	{
		SVMDistanceToDecisionBoundary_Projection(const SVMLearner& learner_, bool use_bound_ = false) : learner(learner_), use_bound(use_bound_)
		{
			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;

			closest_node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				closest_node[i].index = i + 1;
			closest_node[learner.feature_dim].index = -1;

			upper_bound = NULL;
			lower_bound = NULL;

			if(learner.scaler)
			{
				upper_bound = new double[learner.feature_dim];
				lower_bound = new double[learner.feature_dim];

				for(std::size_t i = 0; i < learner.feature_dim; ++i)
				{
					if(learner.use_scaler)
					{
						upper_bound[i] = 2;
						lower_bound[i] = -1;
					}
					else
					{
						double delta = learner.scaler->v_max[i] - learner.scaler->v_min[i];
						upper_bound[i] = learner.scaler->v_max[i] + 0.5 * delta;
						lower_bound[i] = learner.scaler->v_min[i] - 0.5 * delta;
					}
				}
			}

			hyperw = sqrt(learner.hyperw_normsqr);

			index = learner.constructIndexOfSupportVectors();
		}

		~SVMDistanceToDecisionBoundary_Projection()
		{
			delete [] node;
			delete [] closest_node;
			delete [] upper_bound;
			delete [] lower_bound;
			delete index;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const SVMLearner& learner;

		double hyperw;

		mutable svm_node* node;
		mutable struct svm_node* closest_node;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;

		double* upper_bound;
		double* lower_bound;

		bool use_bound;
	};

	struct SVMDistanceToDecisionBoundary_Optimization
	{
		SVMDistanceToDecisionBoundary_Optimization(const SVMLearner& learner_, bool use_bound_ = false) : learner(learner_), use_bound(use_bound_)
		{
			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;

			closest_node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				closest_node[i].index = i + 1;
			closest_node[learner.feature_dim].index = -1;

			upper_bound = NULL;
			lower_bound = NULL;

			if(learner.scaler)
			{
				upper_bound = new double[learner.feature_dim];
				lower_bound = new double[learner.feature_dim];

				for(std::size_t i = 0; i < learner.feature_dim; ++i)
				{
					if(learner.use_scaler)
					{
						upper_bound[i] = 2;
						lower_bound[i] = -1;
					}
					else
					{
						double delta = learner.scaler->v_max[i] - learner.scaler->v_min[i];
						upper_bound[i] = learner.scaler->v_max[i] + 0.5 * delta;
						lower_bound[i] = learner.scaler->v_min[i] - 0.5 * delta;
					}
				}
			}

			index = learner.constructIndexOfSupportVectors();
		}

		~SVMDistanceToDecisionBoundary_Optimization()
		{
			delete [] node;
			delete [] closest_node;
			delete [] upper_bound;
			delete [] lower_bound;
			delete index;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		mutable svm_node* node;
		mutable struct svm_node* closest_node;

		double* upper_bound;
		double* lower_bound;

		const SVMLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;

		bool use_bound;
	};

	struct SVMDistanceToDecisionBoundary_OptimizationGradient
	{
		SVMDistanceToDecisionBoundary_OptimizationGradient(const SVMLearner& learner_, bool use_bound_ = false) : learner(learner_), use_bound(use_bound_)
		{
			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;

			closest_node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				closest_node[i].index = i + 1;
			closest_node[learner.feature_dim].index = -1;

			upper_bound = NULL;
			lower_bound = NULL;

			if(learner.scaler)
			{
				upper_bound = new double[learner.feature_dim];
				lower_bound = new double[learner.feature_dim];

				for(std::size_t i = 0; i < learner.feature_dim; ++i)
				{
					if(learner.use_scaler)
					{
						upper_bound[i] = 2;
						lower_bound[i] = -1;
					}
					else
					{
						double delta = learner.scaler->v_max[i] - learner.scaler->v_min[i];
						upper_bound[i] = learner.scaler->v_max[i] + 0.5 * delta;
						lower_bound[i] = learner.scaler->v_min[i] - 0.5 * delta;
					}
				}
			}

			index = learner.constructIndexOfSupportVectors();
		}

		~SVMDistanceToDecisionBoundary_OptimizationGradient()
		{
			delete [] node;
			delete [] closest_node;
			delete [] upper_bound;
			delete [] lower_bound;
			delete index;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		mutable svm_node* node;
		mutable struct svm_node* closest_node;

		double* upper_bound;
		double* lower_bound;

		const SVMLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;

		bool use_bound;
	};

	struct SVMDistanceToDecisionBoundary_Bruteforce
	{
		SVMDistanceToDecisionBoundary_Bruteforce(const SVMLearner& learner_) : learner(learner_)
		{
			if(!learner.scaler)
			{
				std::cout << "Error! Must have scaler!" << std::endl;
				return;
			}

			double* lower = new double[learner.feature_dim];
			double* upper = new double[learner.feature_dim];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
			{
				lower[i] = learner.use_scaler ? 0 : learner.scaler->v_min[i];
				upper[i] = learner.use_scaler ? 1 : learner.scaler->v_max[i];
			}

			sample_decision_boundary_brute_force(learner.model, upper, lower, learner.feature_dim, 0.005, 100, boundary, index);
			delete [] upper;
			delete [] lower;

			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;

			closest_node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				closest_node[i].index = i + 1;
			closest_node[learner.feature_dim].index = -1;
		}

		~SVMDistanceToDecisionBoundary_Bruteforce()
		{
			delete index;
			delete [] node;
			delete [] closest_node;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const SVMLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;

		std::vector<std::vector<double> > boundary;

		mutable svm_node* node;
		mutable svm_node* closest_node;
	};



	struct MulticonlitronDistanceToDecisionBoundary_BruteForce
	{
		MulticonlitronDistanceToDecisionBoundary_BruteForce(const MulticonlitronLearner& learner_) : learner(learner_)
		{
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const MulticonlitronLearner& learner;
	};

	struct MulticonlitronDistanceToDecisionBoundary_KNN
	{
		MulticonlitronDistanceToDecisionBoundary_KNN(const MulticonlitronLearner& learner_);

		~MulticonlitronDistanceToDecisionBoundary_KNN() { delete index; }

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const MulticonlitronLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;
		std::vector<const HyperPlane*> hyperplanes;
		
		std::size_t data_dim;
	};

	struct MulticonlitronDistanceToDecisionBoundary_Optimization
	{
		MulticonlitronDistanceToDecisionBoundary_Optimization(const MulticonlitronLearner& learner_);

		~MulticonlitronDistanceToDecisionBoundary_Optimization() 
		{ 
			delete index; 
			delete [] upper;
			delete [] lower;
		}

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const MulticonlitronLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;
		std::vector<const HyperPlane*> hyperplanes;

		std::size_t data_dim;

		double* lower;
		double* upper;
	};


	struct MulticonlitronDistanceToDecisionBoundary_EmbedKNN
	{
		MulticonlitronDistanceToDecisionBoundary_EmbedKNN(const MulticonlitronLearner& learner_);

		~MulticonlitronDistanceToDecisionBoundary_EmbedKNN() { delete index; }

		double distance(const DataVector& v) const;

		double distance(const DataVector& v, DataVector& closest_v) const;

		const MulticonlitronLearner& learner;

		flann::Index<FLANN_WRAPPER::DistanceRN>* index;
		std::vector<const HyperPlane*> hyperplanes;

		std::size_t original_dim;
		std::size_t embedded_dim;
	};




	int conlitron_distanceF(long int *Status, long int *n,    double x[],
		long int    *needF,  long int *neF,  double F[],
		long int    *needG,  long int *neG,  double G[],
		char       *cu,     long int *lencu,
		long int    iu[],    long int *leniu,
		double ru[],    long int *lenru );

}

#endif