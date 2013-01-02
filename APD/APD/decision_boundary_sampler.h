#ifndef DECISION_BOUNDARY_SAMPLER_H
#define DECISION_BOUNDARY_SAMPLER_H


#include "contact_space_learning.h"
#include "decision_boundary_distance.h"
#include "learning/spatial_tree.h"
#include "decision_boundary_evaluate.h"

namespace APDL
{
	template<typename Learner, typename Distancer>
	void sample_decision_boundary_hierarchial_tree(const SpatialTreeParam& param,
												   const Learner& learner,
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

		SpatialTree<Learner, Distancer> tree(lower, upper, param, learner);

		tree.collectBoundarySamples(samples);
	}

	//template<typename Distancer>
	//void sample_decision_boundary_hierarchial_tree(const SpatialTreeParam& param,
	//											   const MulticonlitronLearner& learner,
	//	                                           std::vector<DataVector>& samples)
	//{
	//	if(!learner.scaler)
	//	{
	//		std::cout << "Learner does not have scaler!" << std::endl;
	//		return;
	//	}

	//	std::size_t dim = learner.feature_dim;

	//	DataVector upper(dim), lower(dim);
	//	for(std::size_t i = 0; i < dim; ++i)
	//	{
	//		upper[i] = learner.use_scaler ? 1 : learner.scaler->v_max[i];
	//		lower[i] = learner.use_scaler ? 1 : learner.scaler->v_min[i];
	//	}

	//	SpatialTree<MulticonlitronLearner, Distancer> tree(lower, upper, param, learner);

	//	tree.collectBoundarySamples(samples);
	//}



	// difference with sample_decision_boundary_hierarchial_tree: 
	//    we evalute the corners of the node to determine (in fact guess) whether a node is completely inside one region
	template<typename Learner, typename Evaluator>
	void sample_decision_boundary_hierarchial_tree_E(const SpatialTreeEParam& param,
													 const Learner& learner, 
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

		SpatialTreeE<Learner, Evaluator> tree(lower, upper, param, learner);

		tree.collectBoundarySamples(samples);
	}




	//template<typename Evaluator>
	//void sample_decision_boundary_hierarchial_tree_E(const SpatialTreeEParam& param,
	//												 const MulticonlitronLearner& learner, 
	//												 std::vector<DataVector>& samples)
	//{
	//	if(!learner.scaler)
	//	{
	//		std::cout << "Learner does not have scaler!" << std::endl;
	//		return;
	//	}

	//	std::size_t dim = learner.feature_dim;

	//	DataVector upper(dim), lower(dim);
	//	for(std::size_t i = 0; i < dim; ++i)
	//	{
	//		upper[i] = learner.use_scaler ? 1 : learner.scaler->v_max[i];
	//		lower[i] = learner.use_scaler ? 1 : learner.scaler->v_min[i];
	//	}

	//	SpatialTreeE<MulticonlitronLearner, Evaluator> tree(lower, upper, param, learner);

	//	tree.collectBoundarySamples(samples);
	//}


	void sample_decision_boundary_interpolation(const SVMLearner& learner,
		                                        std::vector<DataVector>& samples,
												std::size_t search_num);

	void sample_decision_boundary_interpolation2(const SVMLearner& learner,
		                                         std::vector<DataVector>& samples,
		                                         std::size_t search_num);

	void sample_decision_boundary_interpolation(const MulticonlitronLearner& learner, 
		                                        std::vector<DataVector>& samples,
												std::size_t search_num);



	struct FilterParam
	{
		// for k-means or k-centroid
		std::size_t clustering_max_iter;

		// for filter
		std::size_t filter_threshold;

		FilterParam()
		{
			clustering_max_iter = 10;
			filter_threshold = (std::numeric_limits<double>::max)();
		}
	};


	std::vector<DataVector> sampleSelectionInitialize(const std::vector<DataVector>& samples, std::size_t n);

	std::vector<DataVector> sampleSelectionKMeans(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter);

	std::vector<DataVector> sampleSelectionKCentroids(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter);



	std::vector<DataVector> filter(const SVMEvaluator& eval, std::vector<DataVector>& samples, double threshold);
	std::vector<DataVector> filter(const MulticonlitronEvaluator& eval, std::vector<DataVector>& samples, double threshold);






	template<typename Learner, typename Distancer>
	struct DecisionBoundaryHierarchialTreeSampler
	{
		DecisionBoundaryHierarchialTreeSampler(
			const SpatialTreeParam& param_,
			const FilterParam& fparam_,
			const Learner& learner_) : param(param_), fparam(fparam_), learner(learner_) 
		{}

		std::vector<DataVector> sample(std::size_t n) const
		{
			std::vector<DataVector> samples;
			sample_decision_boundary_hierarchial_tree<Learner, Distancer>(param, learner, samples);
			return sampleSelectionKCentroids(samples, n, fparam.clustering_max_iter);
		}

		SpatialTreeParam param;
		FilterParam fparam;
		const Learner& learner;
	};

	template<typename Learner, typename Evaluator>
	struct DecisionBoundaryHierarchialTreeESampler
	{
		DecisionBoundaryHierarchialTreeESampler(
			const SpatialTreeEParam& param_,
			const FilterParam& fparam_,
			const Learner& learner_) : param(param_), fparam(fparam_), learner(learner_)
		{}

		std::vector<DataVector> sample(std::size_t n) const
		{
			std::vector<DataVector> samples;
			sample_decision_boundary_hierarchial_tree_E<Learner, Evaluator>(param, learner, samples);
			return sampleSelectionKCentroids(samples, n, fparam.clustering_max_iter);
		}
		
		SpatialTreeEParam param;
		FilterParam fparam;
		const Learner& learner;
	};

	template<typename Learner, typename Evaluator>
	struct DecisionBoundaryInterpolationSampler
	{
		DecisionBoundaryInterpolationSampler(std::size_t knn_param_,
											 const FilterParam& fparam_,
											 const Evaluator& evaluator_,
											 const Learner& learner_) : knn_param(knn_param_), fparam(fparam_), evaluator(evaluator_), learner(learner_)
		{}

		std::vector<DataVector> sample(std::size_t n) const
		{
			std::vector<DataVector> samples;
			sample_decision_boundary_interpolation(learner, samples, knn_param);
			return sampleSelectionKCentroids(filter(evaluator, samples, fparam.filter_threshold), n, fparam.clustering_max_iter);
		}

		std::size_t knn_param;
		FilterParam fparam;
		const Evaluator& evaluator;
		const Learner& learner;
	};

	// only for svm
	struct DecisionBoundaryInterpolation2Sampler
	{
		DecisionBoundaryInterpolation2Sampler(std::size_t knn_param_,
											  const FilterParam& fparam_,
											  const SVMEvaluator& evaluator_,
											  const SVMLearner& learner_) : knn_param(knn_param_), fparam(fparam_), evaluator(evaluator_), learner(learner_)
		{}

		std::vector<DataVector> sample(std::size_t n) const
		{
			std::vector<DataVector> samples;
			sample_decision_boundary_interpolation2(learner, samples, knn_param);
			return sampleSelectionKCentroids(filter(evaluator, samples, fparam.filter_threshold), n, fparam.clustering_max_iter);
		}

		std::size_t knn_param;
		FilterParam fparam;
		const SVMEvaluator& evaluator;
		const SVMLearner& learner;
	};






	


}

#endif