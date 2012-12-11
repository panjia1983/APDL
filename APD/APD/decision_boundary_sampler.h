#ifndef DECISION_BOUNDARY_SAMPLER_H
#define DECISION_BOUNDARY_SAMPLER_H


#include "contact_space_learning.h"
#include "decision_boundary_distance.h"
#include "learning/spatial_tree.h"

namespace APDL
{
	template<typename Distancer>
	void sample_decision_boundary_hierarchial_tree(const SVMLearner& learner,
												   std::vector<DataVector>& samples)
	{
		if(!learner.scaler)
		{
			std::cout << "Learner does not have scaler!" << std::endl;
			return;
		}

		std::size_t dim = learner.feature_dim;

		DataVector upper(dim), lower(dim);
		for(std::size_t i = 0; i < dim; ++i)
		{
			upper[i] = learner.use_scaler ? 1 : learner.scaler->v_max[i];
			lower[i] = learner.use_scaler ? 0 : learner.scaler->v_min[i];
		}

		SpatialTree<SVMLearner, Distancer> tree(lower, upper, learner);

		tree.collectBoundarySamples(0.1, samples);
	}

	// difference with sample_decision_boundary_hierarchial_tree: 
	//    we evalute the corners of the node to determine (in fact guess) whether a node is completely inside one region
	template<typename SVMEvaluator>
	void sample_decision_boundary_hierarchial_tree_E(const SpatialTreeEParam& param,
													 const SVMLearner& learner, 
													 std::vector<DataVector>& samples)
	{
		if(!learner.scaler)
		{
			std::cout << "Learner does not have scaler!" << std::endl;
			return;
		}

		std::size_t dim = learner.feature_dim;

		DataVector upper(dim), lower(dim);
		for(std::size_t i = 0; i < dim; ++i)
		{
			upper[i] = learner.use_scaler ? 1 : learner.scaler->v_max[i];
			lower[i] = learner.use_scaler ? 0 : learner.scaler->v_min[i];
		}

		SpatialTreeE<SVMLearner, SVMEvaluator> tree(lower, upper, param, learner);

		tree.collectBoundarySamples(0.1, samples);
	}




	template<typename SVMEvaluator>
	void sample_decision_boundary_hierarchial_tree_E(const SpatialTreeEParam& param,
													 const MulticonlitronLearner& learner, 
													 std::vector<DataVector>& samples)
	{
		if(!learner.scaler)
		{
			std::cout << "Learner does not have scaler!" << std::endl;
			return;
		}

		std::size_t dim = learner.feature_dim;

		DataVector upper(dim), lower(dim);
		for(std::size_t i = 0; i < dim; ++i)
		{
			upper[i] = learner.scaler->v_max[i];
			lower[i] = learner.scaler->v_min[i];
		}

		SpatialTreeE<MulticonlitronLearner, SVMEvaluator> tree(lower, upper, param, learner);

		tree.collectBoundarySamples(0.1, samples);
	}


	void sample_decision_boundary_interpolation(const SVMLearner& learner,
		                                        std::vector<DataVector>& samples,
												std::size_t search_num  = 50);

	void sample_decision_boundary_interpolation(const MulticonlitronLearner& learner, 
		                                        std::vector<DataVector>& samples,
												std::size_t search_num  = 50);




	std::vector<DataVector> sampleSelectionInitialize(const std::vector<DataVector>& samples, std::size_t n);

	std::vector<DataVector> sampleSelectionKMeans(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter = 100);

	std::vector<DataVector> sampleSelectionKCentroids(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter = 100);

}

#endif