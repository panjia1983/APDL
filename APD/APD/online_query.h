#ifndef ONLINE_QUERY_H
#define ONLINE_QUERY_H

#include "contact_space_learning.h"

namespace APDL
{
	struct QueryResult
	{
		DataVector v; // can be used for collision
		bool col;
		double PD;

		QueryResult()
		{
			col = false;
			PD = 0;
		}
	};

	// only use support vector
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query(
		const SVMLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];


		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = learner.model->SV[indices[0][0]][i].value;

		res.v = v;
		res.col = true;
		res.PD = sqrt(dists[0][0]);

		return res;
	}

	// use extended model
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query(
		const SVMLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const std::vector<DataVector>& index_samples,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];


		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = index_samples[indices[0][0]][i];

		res.v = v;
		res.col = true;
		res.PD = sqrt(dists[0][0]);

		return res;
	}

	// only use support vectors
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query2(
		const SVMLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];


		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		svm_node* query_node = new svm_node[dim + 1];
		svm_node* closest_node = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i)
		{
			query_node[i].index = i + 1;
			query_node[i].value = query[i];
			closest_node[i].index = i + 1;
		}
		query_node[dim].index = -1;
		closest_node[dim].index = -1;

		double* initial_guess = new double[dim];
		for(std::size_t i = 0; i < dim; ++i)
			initial_guess[i] = learner.model->SV[indices[0][0]][i].value;

		double d = dist_to_decision_boundary_with_gradient(learner.model, query_node, initial_guess, NULL, NULL, closest_node);

		DataVector v(contactspace.data_dim());
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = closest_node[i].value;

		ContactSpace::DistanceType metric;
		res.v = v;
		res.col = true;
		DataVector query_tmp(query);
		res.PD = sqrt(metric(&v[0], &query_tmp[0], dim));

		delete [] query_node;
		delete [] closest_node;
		delete [] initial_guess;

		return res;
	}

	// use extended model
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query2(
		const SVMLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const std::vector<DataVector>& index_samples,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];


		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		svm_node* query_node = new svm_node[dim + 1];
		svm_node* closest_node = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i)
		{
			query_node[i].index = i + 1;
			query_node[i].value = query[i];
			closest_node[i].index = i + 1;
		}
		query_node[dim].index = -1;
		closest_node[dim].index = -1;

		double* initial_guess = new double[dim];
		for(std::size_t i = 0; i < dim; ++i)
			initial_guess[i] = index_samples[indices[0][0]][i];

		double d = dist_to_decision_boundary_with_gradient(learner.model, query_node, initial_guess, NULL, NULL, closest_node);

		DataVector v(contactspace.data_dim());
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = closest_node[i].value;

		ContactSpace::DistanceType metric;
		res.v = v;
		res.col = true;
		DataVector query_tmp(query);
		res.PD = sqrt(metric(&v[0], &query_tmp[0], dim));

		delete [] query_node;
		delete [] closest_node;
		delete [] initial_guess;

		return res;
	}

	// use support vectors
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query(
		const MulticonlitronLearner& learner, const ContactSpace& contactspace,
		IndexType<typename ContactSpace::DistanceType>* index,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		std::size_t id = indices[0][0];
		if(id % 2) id = (id - 1) / 2;
		else id = id / 2;
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = 0.5 * (learner.hyperplanes[id]->supp1[i] + learner.hyperplanes[id]->supp2[i]);

		res.v = v;
		res.col = true;
		res.PD = sqrt(dists[0][0]);
	}

	// use extended model
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query(
		const MulticonlitronLearner& learner, const ContactSpace& contactspace,
		IndexType<typename ContactSpace::DistanceType>* index,
		const std::vector<DataVector>& index_samples,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		std::size_t id = indices[0][0];
		for(std::size_t i = 0; i < dim; ++i)
			v[i] = index_samples[id][i];

		res.v = v;
		res.col = true;
		res.PD = sqrt(dists[0][0]);
	}

	// use support vectors
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query2(
		const MulticonlitronLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		std::size_t id = indices[0][0];
		if(id % 2)
			id = (id - 1) / 2;
		else
			id = id / 2;

		double* initial_guess = new double[dim];
		double* x = new double[dim];
		double* closest_x = new double[dim];
		for(std::size_t i = 0; i < dim; ++i)
		{
			initial_guess[i] = 0.5 * (learner.hyperplanes[id]->supp2[i] + learner.hyperplanes[id]->supp1[i]);
			x[i] = query[i];
		}

		conlitron_dist_to_decision_boundary((void*)&learner.model, x, initial_guess, NULL, NULL, closest_x, conlitron_distanceF, dim);

		for(std::size_t i = 0; i < dim; ++i)
			v[i] = closest_x[i];

		ContactSpace::DistanceType metric;
		res.v = v;
		res.col = true;
		DataVector query_tmp(query);
		res.PD = sqrt(metric(&v[0], &query_tmp[0], dim));

		delete [] x;
		delete [] closest_x;
		delete [] initial_guess;

		return res;
	}

	// use extended model
	template<typename ContactSpace, template <typename> class IndexType>
	QueryResult PD_query2(
		const MulticonlitronLearner& learner, const ContactSpace& contactspace, 
		IndexType<typename ContactSpace::DistanceType>* index,
		const std::vector<DataVector>& index_samples,
		const DataVector& query)
	{
		QueryResult res;
		if(learner.predict(query).label < 0)
			return res;

		std::size_t dim = learner.feature_dim;
		flann::Matrix<double> queryset = flann::Matrix<double>(new double[dim], 1, dim);
		for(std::size_t i = 0; i < dim; ++i)
			queryset[0][i] = query[i];

		std::vector<std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		index->knnSearch(queryset, indices, dists, 1, flann::SearchParams());

		DataVector v(contactspace.data_dim());
		std::size_t id = indices[0][0];

		double* initial_guess = new double[dim];
		double* x = new double[dim];
		double* closest_x = new double[dim];
		for(std::size_t i = 0; i < dim; ++i)
		{
			initial_guess[i] = index_samples[id][i];
			x[i] = query[i];
		}

		conlitron_dist_to_decision_boundary((void*)&learner.model, x, initial_guess, NULL, NULL, closest_x, conlitron_distanceF, dim);

		for(std::size_t i = 0; i < dim; ++i)
			v[i] = closest_x[i];

		ContactSpace::DistanceType metric;
		res.v = v;
		res.col = true;
		DataVector query_tmp(query);
		res.PD = sqrt(metric(&v[0], &query_tmp[0], dim));

		delete [] x;
		delete [] closest_x;
		delete [] initial_guess;

		return res;
	}
}

#endif