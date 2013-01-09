#include "decision_boundary_sampler.h"


namespace APDL
{
	std::vector<DataVector> sampleSelectionInitialize(const std::vector<DataVector>& samples, std::size_t n)
	{
		// assume n > samples.size()
		std::size_t dim = samples[0].dim();
		std::list<std::size_t> non_selected;
		std::list<std::size_t> selected;
		selected.push_back(0);
		for(std::size_t i = 1; i < samples.size(); ++i) non_selected.push_back(i);

		for(std::size_t i = 1; i < n; ++i)
		{
			std::list<std::size_t>::iterator max_it;
			double max_dist = -1;

			for(std::list<std::size_t>::iterator j_it = non_selected.begin(); j_it != non_selected.end(); ++j_it)
			{
				std::size_t id_non_selected = *j_it;

				double min_dist_to_selected = (std::numeric_limits<double>::max)();

				for(std::list<std::size_t>::iterator k_it = selected.begin(); k_it != selected.end(); ++k_it)
				{
					std::size_t id_selected = *k_it;

					double dist = 0;
					for(std::size_t d = 0; d < dim; ++d)
					{
						double v = samples[id_selected][d] - samples[id_non_selected][d];
						dist += v * v;
					}

					if(dist < min_dist_to_selected) min_dist_to_selected = dist;
				}

				if(min_dist_to_selected > max_dist) { max_dist = min_dist_to_selected; max_it = j_it; }
			}

			selected.push_back(*max_it);
			non_selected.erase(max_it);
		}

		std::vector<DataVector> res;
		for(std::list<std::size_t>::iterator it = selected.begin(); it != selected.end(); ++it)
		{
			res.push_back(samples[*it]);
		}

		return res;
	}


	std::vector<DataVector> sampleSelectionKMeans(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter)
	{
		if(samples.size() <= n) return samples;
		std::size_t dim = samples[0].dim();
		std::vector<DataVector> centers = sampleSelectionInitialize(samples, n);

		std::size_t n_clusters = centers.size();

		// clusters[i] contains the id of points in the i-th cluster
		std::vector<std::vector<std::size_t> > clusters(n_clusters);
		for(std::size_t iter = 0; iter < max_iter; ++iter)
		{
			for(std::size_t i = 0; i < samples.size(); ++i)
			{
				std::size_t min_cluster_id = -1;
				double min_cluster_dist = (std::numeric_limits<double>::max)();

				for(std::size_t j = 0; j < centers.size(); ++j)
				{
					double dist = 0;
					for(std::size_t d = 0; d < dim; ++d)
					{
						double v = centers[j][d] - samples[i][d];
						dist += v * v;
					}

					if(dist < min_cluster_dist) { min_cluster_dist = dist; min_cluster_id = j; }
				}

				clusters[min_cluster_id].push_back(i);
			}

			// for each cluster, compute the mean as the new center
			for(std::size_t i = 0; i < clusters.size(); ++i)
			{
				if(clusters[i].size() == 0) continue;

				DataVector mean(dim);
				mean.setZero();

				for(std::size_t j = 0; j < clusters[i].size(); ++j)
					mean += samples[clusters[i][j]];
				mean *= (1.0 / clusters[i].size());

				centers[i] = mean;
			}
		}

		return centers;
	}

	std::vector<DataVector> sampleSelectionKCentroids(std::vector<DataVector>& samples, std::size_t n, std::size_t max_iter)
	{
		if(samples.size() <= n) return samples;
		std::size_t dim = samples[0].dim();
		std::vector<DataVector> centers = sampleSelectionInitialize(samples, n);

		std::size_t n_clusters = centers.size();

		// clusters[i] contains the id of points in the i-th cluster
		std::vector<std::vector<std::size_t> > clusters(n_clusters);

		double old_cluster_cost;
		double new_cluster_cost;

		for(std::size_t iter = 0; iter < max_iter; ++iter)
		{
			if(iter != 0) old_cluster_cost = new_cluster_cost;

			new_cluster_cost = 0;

			for(std::size_t i = 0; i < samples.size(); ++i)
			{
				std::size_t min_cluster_id = -1;
				double min_cluster_dist = (std::numeric_limits<double>::max)();

				for(std::size_t j = 0; j < centers.size(); ++j)
				{
					double dist = 0;
					for(std::size_t d = 0; d < dim; ++d)
					{
						double v = centers[j][d] - samples[i][d];
						dist += v * v;
					}

					if(dist < min_cluster_dist) { min_cluster_dist = dist; min_cluster_id = j; }
				}

				clusters[min_cluster_id].push_back(i);

				new_cluster_cost += min_cluster_dist; 
			}

			if(iter != 0)
			{
				if(new_cluster_cost >= old_cluster_cost) break;
			}


			// for each cluster, compute the mean and then use the point closest to the mean as the new center
			for(std::size_t i = 0; i < clusters.size(); ++i)
			{
				if(clusters[i].size() == 0) continue;

				DataVector mean(dim);
				mean.setZero();

				for(std::size_t j = 0; j < clusters[i].size(); ++j)
					mean += samples[clusters[i][j]];
				mean *= (1.0 / clusters[i].size());

				double min_dist = (std::numeric_limits<double>::max)();
				std::size_t min_id = -1;
				for(std::size_t j = 0; j < clusters[i].size(); ++j)
				{
					std::size_t id = clusters[i][j];
					double dist = 0;
					for(std::size_t d = 0; d < dim; ++d)
					{
						double v = samples[id][d] - mean[d];
						dist += v * v;
					}

					if(dist < min_dist) { min_dist = dist; min_id = id; }
				}

				centers[i] = samples[min_id];
			}
		}

		return centers;
	}

	void sample_decision_boundary_interpolation2(const SVMLearner& learner,
												std::vector<DataVector>& samples,
												std::size_t search_num)
	{
		const struct svm_model* model = learner.model;
		std::size_t dim = learner.feature_dim;
		double inv_w_norm = sqrt(learner.hyperw_normsqr);

		double* upper = new double[dim];
		double* lower = new double[dim];
		for(std::size_t i = 0; i < dim; ++i)
		{
			upper[i] = learner.scaler->v_min[i];
			lower[i] = learner.scaler->v_max[i];
		}

		std::size_t start[2];
		std::size_t number[2];
		start[0] = 0;
		number[0] = model->nSV[0];
		start[1] = start[0] + number[0];
		number[1] = model->nSV[1];

		if(number[1] < number[0])
		{
			std::size_t tmp = number[1];
			number[1] = number[0];
			number[0] = tmp;

			tmp = start[1];
			start[1] = start[0];
			start[0] = tmp;
		}

		svm_node* midpoint = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i) midpoint[i].index = i + 1;
		midpoint[dim].index = -1;

		DataVector p(dim);

		flann::Index<FLANN_WRAPPER::DistanceRN>* index = NULL;
		if(number[1] < number[0]) 
			index = learner.constructIndexOfSupportVectorsClass0();
		else 
			index = learner.constructIndexOfSupportVectorsClass1();


		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = start[0]; i < start[0] + number[0]; ++i)
		{
			svm_node* x = model->SV[i];
			for(std::size_t j = 0; j < dim; ++j)
				queryset[0][j] = x[j].value;

			std::vector<std::vector<int> > indices;
			std::vector<std::vector<double> > dists;

			index->knnSearch(queryset, indices, dists, search_num, flann::SearchParams());

			for(std::size_t j = 0; j < indices[0].size(); ++j)
			{
				svm_node* y = model->SV[start[1] + indices[0][j]];

				feature_space_midpoint_with_gradient(model, x, y, midpoint, upper, lower);

				//double f_value = svm_predict_values_twoclass(model, midpoint);
				//if(abs(f_value) > 1) continue;

				for(std::size_t k = 0; k < dim; ++k)
					p[k] = midpoint[k].value;

				samples.push_back(p);
			}
		}

		delete index;
		delete [] midpoint;
		delete [] upper;
		delete [] lower;
	}

	void sample_decision_boundary_interpolation(const SVMLearner& learner,
												std::vector<DataVector>& samples,
												std::size_t search_num)
	{
		const struct svm_model* model = learner.model;
		std::size_t dim = learner.feature_dim;

		std::size_t start[2];
		std::size_t number[2];
		start[0] = 0;
		number[0] = model->nSV[0];
		start[1] = start[0] + number[0];
		number[1] = model->nSV[1];

		if(number[1] < number[0])
		{
			std::size_t tmp = number[1];
			number[1] = number[0];
			number[0] = tmp;

			tmp = start[1];
			start[1] = start[0];
			start[0] = tmp;
		}

		svm_node* midpoint = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i) midpoint[i].index = i + 1;
		midpoint[dim].index = -1;

		DataVector p(dim);

		flann::Index<FLANN_WRAPPER::DistanceRN>* index = NULL;
		if(number[1] < number[0]) 
			index = learner.constructIndexOfSupportVectorsClass0();
		else 
			index = learner.constructIndexOfSupportVectorsClass1();

		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = start[0]; i < start[0] + number[0]; ++i)
		{
			svm_node* x = model->SV[i];
			for(std::size_t j = 0; j < dim; ++j)
				queryset[0][j] = x[j].value;

			std::vector<std::vector<int> > indices;
			std::vector<std::vector<double> > dists;

			index->knnSearch(queryset, indices, dists, search_num, flann::SearchParams());

			for(std::size_t j = 0; j < indices[0].size(); ++j)
			{
				svm_node* y = model->SV[start[1] + indices[0][j]];

				// feature_space_midpoint(model, x, y, midpoint);

				for(std::size_t k = 0; k < dim; ++k)
					midpoint[k].value = 0.5 * (x[k].value + y[k].value);

				// double f_value = svm_predict_values_twoclass(model, midpoint);
				// if(abs(f_value) > 1) continue;

				for(std::size_t k = 0; k < dim; ++k)
					p[k] = midpoint[k].value;

				samples.push_back(p);
			}
		}

		delete index;

		delete [] midpoint;
	}


	void sample_decision_boundary_interpolation(const MulticonlitronLearner& learner, 
		                                        std::vector<DataVector>& samples,
												std::size_t search_num)
	{
		std::size_t dim = learner.feature_dim;
		std::size_t data_num = learner.model.numOfHyperPlanes();
		flann::Matrix<double> dataset = flann::Matrix<double>(new double[dim * data_num], data_num, dim);
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim * data_num], data_num, dim);
		std::size_t id = 0;
	
		for(std::size_t i = 0; i < learner.model.size(); ++i)
		{
			for(std::size_t j = 0; j < learner.model[i].size(); ++j)
			{
				for(std::size_t d = 0; d < dim; ++d)
				{
					dataset[id][d] = learner.model[i][j].supp1[d];
					queryset[id][d] = learner.model[i][j].supp2[d];
				}
				id++;
			}
		}	



		flann::Index<FLANN_WRAPPER::DistanceRN>* index = new flann::Index<FLANN_WRAPPER::DistanceRN>(dataset, flann::KDTreeIndexParams());
		index->buildIndex();

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;

		index->knnSearch(queryset, indices, dists, search_num, flann::SearchParams());


		for(std::size_t i = 0; i < indices.size(); ++i)
		{
			for(std::size_t j = 0; j < indices[i].size(); ++j)
			{
				DataVector p(dim);
				for(std::size_t d = 0; d < dim; ++d)
					p[d] = 0.5 * (queryset[i][d] + dataset[indices[i][j]][d]);

				samples.push_back(p);
			}
		}

		delete index;
	}


	std::vector<DataVector> filter(const SVMEvaluator& eval, std::vector<DataVector>& samples, double threshold)
	{
		std::vector<DataVector> res;
		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			if(abs(eval.evaluate(samples[i])) < threshold)
				res.push_back(samples[i]);
		}

		return res;
	}

	std::vector<DataVector> filter(const MulticonlitronEvaluator& eval, std::vector<DataVector>& samples, double threshold)
	{
		std::vector<DataVector> res;
		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			if(abs(eval.evaluate(samples[i])) < threshold)
				res.push_back(samples[i]);
		}

		return res;
	}
}