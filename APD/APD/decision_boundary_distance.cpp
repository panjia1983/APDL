#include "decision_boundary_distance.h"

namespace APDL
{
	struct CompareDistToDecisionBoundaryLess
	{
		bool operator ()(const std::pair<std::size_t, double>& a, const std::pair<std::size_t, double>& b) const 
		{
			return std::abs(a.second) < std::abs(b.second);
		}
	};

	void sample_decision_boundary_brute_force(const struct svm_model* model, 
											  double* upper, double* lower,
											  std::size_t dim,
											  double on_decision_threshold, 
											  std::size_t subdivision_level, 
											  std::vector<std::vector<double> >& boundary)
	{
		int max_i = 1;
		for(std::size_t i = 0; i < dim; ++i) max_i *= subdivision_level;

		std::vector<std::size_t> id(dim);
		struct svm_node* node = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i)
			node[i].index = i + 1;
		node[dim].index = -1;

		std::vector<std::pair<std::size_t, double> > cached_results(max_i);

		for(std::size_t i = 0; i < max_i; ++i)
		{
			std::size_t i_ = i;
			for(std::size_t j = 0; j < dim; ++j)
			{
				id[j] = i_ % subdivision_level;
				i_ = i_ / subdivision_level;	
			}

			for(std::size_t j = 0; j < dim; ++j)
				node[j].value = (upper[j] - lower[j]) * id[j] / (double)(subdivision_level - 1) + lower[j];

			double val = svm_predict_values_twoclass(model, node);

			cached_results[i] = std::make_pair(i, val);
		}

		std::sort(cached_results.begin(), cached_results.end(), CompareDistToDecisionBoundaryLess());

		std::size_t n_threshold = ceil(on_decision_threshold * max_i);
		if(n_threshold == 0) n_threshold = 1;
		else if(n_threshold >= max_i) n_threshold = max_i;

		for(std::size_t i = 0; i < n_threshold; ++i)
		{
			std::size_t i_ = cached_results[i].first;
			for(std::size_t j = 0; j < dim; ++j)
			{
				id[j] = i_ % subdivision_level;
				i_ = i_ / subdivision_level;	
			}

			std::vector<double> p(dim);
			for(std::size_t j = 0; j < dim; ++j)
				p[j] = (upper[j] - lower[j]) * id[j] / (double)(subdivision_level - 1) + lower[j];

			// std::cout << cached_results[i].second << std::endl;
			boundary.push_back(p);
		}

		delete [] node;
	}

	void sample_decision_boundary_brute_force(const struct svm_model* model,
		double* upper, double* lower,
		std::size_t dim,
		double on_decision_threshold, 
		std::size_t subdivision_level, 
		std::vector<std::vector<double> >& boundary,
		flann::Index<FLANN_WRAPPER::DistanceRN>*& index)
	{
		APDL::sample_decision_boundary_brute_force(model, upper, lower, dim, on_decision_threshold, subdivision_level, boundary);

		std::size_t num_data = boundary.size();
		std::size_t dim_data = dim;
		flann::Matrix<double> dataset = flann::Matrix<double>(new double[num_data * dim_data], num_data, dim_data);

		for(std::size_t i = 0; i < num_data; ++i)
		{
			for(std::size_t j = 0; j < dim_data; ++j)
			{
				dataset[i][j] = boundary[i][j];
			}
		}

		index = new flann::Index<FLANN_WRAPPER::DistanceRN>(dataset, flann::KDTreeIndexParams());
		index->buildIndex();
	}


	double dist_to_decision_boundary_brute_force(const std::vector<std::vector<double> >& boundary, const struct svm_node* x, struct svm_node* x_res)
	{
		if(boundary.size() == 0) return 0;

		std::size_t dim = boundary[0].size();
		double min_dist = (std::numeric_limits<double>::max)();
		std::size_t min_id = -1;
		for(std::size_t i = 0; i < boundary.size(); ++i)
		{
			double dist = 0;
			for(std::size_t j = 0; j < dim; ++j)
			{
				double v = (x[j].value - boundary[i][j]);
				dist += v * v;
			}
			if(min_dist > dist) { min_dist = dist; min_id = i; }
		}

		if(x_res)
		{
			for(int i = 0; i < dim; ++i)
				x_res[i].value = boundary[min_id][i];
		}

		return std::sqrt(min_dist);
	}

	double dist_to_decision_boundary_brute_force(const std::vector<std::vector<double> >& boundary, flann::Index<FLANN_WRAPPER::DistanceRN>* index, const struct svm_node* x, struct svm_node* x_res)
	{
		if(boundary.size() == 0) return 0;

		std::size_t dim = boundary[0].size();

		std::size_t dim_data = dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim_data], 1, dim_data);

		for(std::size_t i = 0; i < dim_data; ++i)
		{
			queryset[0][i] = x[i].value;
		}

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		if(x_res)
		{
			for(int i = 0; i < dim; ++i)
				x_res[i].value = boundary[indices[0][0]][i];
		}

		return sqrt(dists[0][0]);
	}



	double SVMDistanceToDecisionBoundary_RoughLowerBound::distance(const DataVector& v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		double predict_label = svm_predict_values_twoclass(learner.model, node);
		double f_dist = predict_label * inv_hyperw;
		double i_dist = sqrt(-1/learner.param.gamma * log(1 - f_dist * f_dist * 0.5));

		return i_dist;
	}

	double SVMDistanceToDecisionBoundary_Projection::distance(const DataVector& v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		
		if(!use_bound)
			return dist_to_decision_boundary_constrain_free(learner.model, hyperw, node);
		else
			return dist_to_decision_boundary_constrain_free(learner.model, hyperw, node, upper_bound, lower_bound);
	}

	double SVMDistanceToDecisionBoundary_Projection::distance(const DataVector& v, DataVector& closest_v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		double d;

		if(!use_bound)
			dist_to_decision_boundary_constrain_free(learner.model, hyperw, node, NULL, NULL, closest_node);
		else 
			dist_to_decision_boundary_constrain_free(learner.model, hyperw, node, upper_bound, lower_bound, closest_node);
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			closest_v[i] = closest_node[i].value;

		return d;
	}

	double SVMDistanceToDecisionBoundary_Optimization::distance(const DataVector& v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];

		if(!use_bound)
			return dist_to_decision_boundary(learner.model, node);
		else
			return dist_to_decision_boundary(learner.model, node, upper_bound, lower_bound);
	}

	double SVMDistanceToDecisionBoundary_Optimization::distance(const DataVector& v, DataVector& closest_v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		double d;
		
		if(!use_bound)
			dist_to_decision_boundary(learner.model, node, NULL, NULL, closest_node);
		else 
			dist_to_decision_boundary(learner.model, node, upper_bound, lower_bound, closest_node);
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			closest_v[i] = closest_node[i].value;

		return d;
	}

	double SVMDistanceToDecisionBoundary_OptimizationGradient::distance(const DataVector& v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];

		if(!use_bound)
			return dist_to_decision_boundary_with_gradient(learner.model, node);
		else
			return dist_to_decision_boundary_with_gradient(learner.model, node, upper_bound, lower_bound);
	}

	double SVMDistanceToDecisionBoundary_OptimizationGradient::distance(const DataVector& v, DataVector& closest_v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		double d;
		
		if(!use_bound)
			dist_to_decision_boundary_with_gradient(learner.model, node, NULL, NULL, closest_node);
		else 
			dist_to_decision_boundary_with_gradient(learner.model, node, upper_bound, lower_bound, closest_node);
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			closest_v[i] = closest_node[i].value;

		return d;
	}

	double SVMDistanceToDecisionBoundary_Bruteforce::distance(const DataVector& v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];

		return dist_to_decision_boundary_brute_force(boundary, index, node);
	}

	double SVMDistanceToDecisionBoundary_Bruteforce::distance(const DataVector& v, DataVector& closest_v) const
	{
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			node[i].value = v[i];
		double d = dist_to_decision_boundary_brute_force(boundary, index, node, closest_node);
		for(std::size_t i = 0; i < learner.feature_dim; ++i)
			closest_v[i] = closest_node[i].value;

		return d;
	}


	double MulticonlitronDistanceToDecisionBoundary_BruteForce::distance(const DataVector& v) const
	{
		const HyperPlane* best_hp = NULL;
		double min_dist = (std::numeric_limits<double>::max)();

		for(std::size_t i = 0; i < learner.model.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.model[i].size(); ++j)
			{
				DataVector center = (learner.model[i][j].supp1 + learner.model[i][j].supp2) * 0.5;
				
				double dist = dot_prod(center - v, center - v);
				if(min_dist > dist) 
				{
					min_dist = dist;
					best_hp = &learner.model[i][j];
				}
			}
		}

		return best_hp->distance(v);
	}

	double MulticonlitronDistanceToDecisionBoundary_BruteForce::distance(const DataVector& v, DataVector& closest_v) const
	{
		const HyperPlane* best_hp = NULL;
		double min_dist = (std::numeric_limits<double>::max)();

		for(std::size_t i = 0; i < learner.model.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.model[i].size(); ++j)
			{
				DataVector center = (learner.model[i][j].supp1 + learner.model[i][j].supp2) * 0.5;
				
				double dist = dot_prod(center - v, center - v);
				if(min_dist > dist) 
				{
					min_dist = dist;
					best_hp = &learner.model[i][j];
				}
			}
		}

		closest_v = best_hp->project(v);

		return best_hp->distance(v);
	}

	MulticonlitronDistanceToDecisionBoundary_KNN::MulticonlitronDistanceToDecisionBoundary_KNN(const MulticonlitronLearner& learner_) : learner(learner_)
	{
		if(learner.model.size() == 0) return;

		data_dim = learner.model[0][0].dim();
		std::size_t data_num = learner.model.numOfHyperPlanes();

		flann::Matrix<double> dataset = flann::Matrix<double>(new double[data_dim * data_num], data_num, data_dim);

		std::size_t id = 0;
		for(std::size_t i = 0; i < learner.model.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.model[i].size(); ++j)
			{
				DataVector p = 0.5 * (learner.model[i][j].supp1 + learner.model[i][j].supp2);
				hyperplanes.push_back(&learner.model[i][j]);

				for(std::size_t k = 0; k < data_dim; ++k)
					dataset[id][k] = p[k];

				id++;
			}
		}

		index = new flann::Index<FLANN_WRAPPER::DistanceRN>(dataset, flann::KDTreeIndexParams());
		index->buildIndex();
	}

	double MulticonlitronDistanceToDecisionBoundary_KNN::distance(const DataVector& v) const
	{
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[data_dim], 1, data_dim);
		for(std::size_t i = 0; i < data_dim; ++i) queryset[0][i] = v[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		std::size_t min_id = indices[0][0];

		return hyperplanes[min_id]->distance(v);
	}

	double MulticonlitronDistanceToDecisionBoundary_KNN::distance(const DataVector& v, DataVector& closest_v) const
	{
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[data_dim], 1, data_dim);
		for(std::size_t i = 0; i < data_dim; ++i) queryset[0][i] = v[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		std::size_t min_id = indices[0][0];

		closest_v = hyperplanes[min_id]->project(v);

		return hyperplanes[min_id]->distance(v);
	}

	MulticonlitronDistanceToDecisionBoundary_EmbedKNN::MulticonlitronDistanceToDecisionBoundary_EmbedKNN(const MulticonlitronLearner& learner_) : learner(learner_)
	{
		if(learner.model.size() == 0) return;
		
		std::size_t data_dim = learner.model[0][0].dim();
		std::size_t data_num = learner.model.numOfHyperPlanes();
		original_dim = data_dim;
		embedded_dim = (data_dim + 2) * (data_dim + 3) / 2;

		std::vector<DataVector> embedded_hyperplanes;
		double max_b_sqr = 0;
		for(std::size_t i = 0; i < learner.model.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.model[i].size(); ++j)
			{
				const DataVector& w = learner.model[i][j].w;
				double b = learner.model[i][j].b;

				DataVector embedded_p(data_dim + 2);
				for(std::size_t k = 0; k < data_dim; ++k)
					embedded_p[k] = w[k];
				embedded_p[data_dim] = b;

				double b_sqr = b * b;
				if(b_sqr > max_b_sqr) max_b_sqr = b_sqr;

				embedded_hyperplanes.push_back(embedded_p);
				hyperplanes.push_back(&learner.model[i][j]);
			}
		}

		for(std::size_t i = 0; i < embedded_hyperplanes.size(); ++i)
			embedded_hyperplanes[i][data_dim + 1] = sqrt(max_b_sqr - embedded_hyperplanes[i][data_dim] * embedded_hyperplanes[i][data_dim]);


		for(std::size_t i = 0; i < embedded_hyperplanes.size(); ++i)
		{
			DataVector embedded_v(embedded_dim);
			
			int id = 0;
			for(std::size_t j = 0; j < (data_dim + 2); ++j)
			{
				for(std::size_t k = j; k < data_dim + 2; ++k)
				{
					double ratio = (k == j) ? (1.0 / sqrt(2.0)) : 1;
					embedded_v[id] = embedded_hyperplanes[i][j] * embedded_hyperplanes[i][k] * ratio;
					id++;
				}
			}

			embedded_hyperplanes[i] = embedded_v;
		}

		flann::Matrix<double> dataset = flann::Matrix<double>(new double[embedded_dim * embedded_hyperplanes.size()], embedded_hyperplanes.size(), embedded_dim);

		for(std::size_t i = 0; i < embedded_hyperplanes.size(); ++i)
		{
			for(std::size_t j = 0; j < embedded_dim; ++j)
				dataset[i][j] = embedded_hyperplanes[i][j];
		}

		index = new flann::Index<FLANN_WRAPPER::DistanceRN>(dataset, flann::KDTreeIndexParams());
		index->buildIndex();


	}

	double MulticonlitronDistanceToDecisionBoundary_EmbedKNN::distance(const DataVector& v) const
	{
		DataVector extended_v(original_dim + 2);
		for(std::size_t i = 0; i < original_dim; ++i) extended_v[i] = v[i];
		extended_v[original_dim] = 1;
		extended_v[original_dim + 1] = 0;

		flann::Matrix<double> queryset = flann::Matrix<double>(new double[embedded_dim], 1, embedded_dim);

		std::size_t id = 0;

		for(std::size_t i = 0; i < original_dim + 2; ++i)
		{
			for(std::size_t j = i; j < original_dim + 2; ++j)
			{
				double ratio = (i == j) ? (1.0 / sqrt(2.0)) : 1;
				queryset[0][id] = -ratio * extended_v[i] * extended_v[j];
				id++;
			}
		}

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		std::size_t min_id = indices[0][0];

		return hyperplanes[min_id]->distance(v);
	}

	double MulticonlitronDistanceToDecisionBoundary_EmbedKNN::distance(const DataVector& v, DataVector& closest_v) const
	{
		DataVector extended_v(original_dim + 2);
		for(std::size_t i = 0; i < original_dim; ++i) extended_v[i] = v[i];
		extended_v[original_dim] = 1;
		extended_v[original_dim + 1] = 0;

		flann::Matrix<double> queryset = flann::Matrix<double>(new double[embedded_dim], 1, embedded_dim);

		std::size_t id = 0;

		for(std::size_t i = 0; i < original_dim + 2; ++i)
		{
			for(std::size_t j = i; j < original_dim + 2; ++j)
			{
				double ratio = (i == j) ? (1.0 / sqrt(2.0)) : 1;
				queryset[0][id] = -ratio * extended_v[i] * extended_v[j];
				id++;
			}
		}

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		std::size_t min_id = indices[0][0];

		closest_v = hyperplanes[min_id]->project(v);

		return hyperplanes[min_id]->distance(v);
	}




}