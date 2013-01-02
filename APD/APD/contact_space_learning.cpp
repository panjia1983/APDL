#include "contact_space_learning.h"
#include "distance_proxy.h"



namespace APDL
{
	PredictResult MulticonlitronLearner::predict(const DataVector& query) const
	{
		PredictResult predict_result;

		const DataVector& v = (scaler && use_scaler) ? scaler->scale(query) : query;

		double predict_label = model.evaluate(v);
		if(predict_label > 0) predict_label = 1;
		else predict_label = 0;

		return PredictResult(predict_label, 1);
	}

	std::vector<PredictResult> MulticonlitronLearner::predict(const std::vector<DataVector>& queries) const
	{
		std::vector<PredictResult> predict_results;

		for(std::size_t i = 0; i < queries.size(); ++i)
		{
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(queries[i]) : queries[i];
			double predict_label = model.evaluate(v);
			if(predict_label > 0) predict_label = 1;
			else predict_label = 0;

			predict_results.push_back(PredictResult(predict_label, 1));
		}

		return predict_results;
	}

	std::vector<PredictResult> MulticonlitronLearner::predict(const std::vector<ContactSpaceSampleData>& queries) const
	{
		std::vector<PredictResult> predict_results;

		for(std::size_t i = 0; i < queries.size(); ++i)
		{
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(queries[i].v) : queries[i].v;
			double predict_label = model.evaluate(v);
			if(predict_label > 0) predict_label = 1;
			else predict_label = 0;
			
			predict_results.push_back(PredictResult(predict_label, 1));
		}

		return predict_results;
	}


	void MulticonlitronLearner::learn(const std::vector<ContactSpaceSampleData>& data, std::size_t active_dim)
	{
		if(active_dim  != feature_dim)
		{
			std::cout << "The dimension of the learner is not correct!" << std::endl;
			return;
		}

		std::vector<DataVector> X, Y;
		for(std::size_t i = 0; i < data.size(); ++i)
		{
			const DataVector& data_v = (scaler && use_scaler) ? scaler->scale(data[i].v) : data[i].v;
			DataVector v(distancer.dim());
			for(std::size_t j = 0; j < distancer.dim(); ++j)
				v[j] = data_v[j];

			if(data[i].col) 
				X.push_back(v);
			else
				Y.push_back(v);
		}


		if(use_approximate_dist)
		{
			if(knn_index) delete knn_index;
			generateIndex<FLANN_WRAPPER::DistanceRN, DataVector>(Y, active_dim, knn_index, flann::KDTreeIndexParams());
		}

		model = SMA(X, Y);

		hyperplanes.clear();

		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				hyperplanes.push_back(&(model[i][j]));
			}
		}
	}

	void MulticonlitronLearner::incremental_learn(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
	{
		if(data_.size() == 0) return;
		std::vector<ContactSpaceSampleData> data = data_;


		DataVector v(data[0].v.dim());
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				const DataVector& scaled_supp1 = (scaler && use_scaler) ? scaler->unscale(model[i][j].supp1) : scaler->unscale(model[i][j].supp1);
				for(std::size_t k = 0; k < active_dim; ++k)
					v[k] = scaled_supp1[k];
				data.push_back(ContactSpaceSampleData(v, true));

				const DataVector& scaled_supp2 = (scaler && use_scaler) ? scaler->unscale(model[i][j].supp2) : scaler->unscale(model[i][j].supp1);
				for(std::size_t k = 0; k < active_dim; ++k)
					v[k] = scaled_supp2[k];
				data.push_back(ContactSpaceSampleData(v, false));
			}
		}

		learn(data, active_dim);
	}


	HyperPlane MulticonlitronLearner::CDMA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const
	{
		if(X.size() == 1 && Y.size() == 1)
		{

			DataVector x_star(X[0]);
			DataVector y_star(Y[0]);

			do
			{
				DataVector x(x_star), y(y_star);

				DataVector best_x(X[0]);
				double best_d = distancer.sqrDistance(best_x, y_star);

				double d12 = distancer.sqrDistance(X[0], x);
				if(d12 > 0)
				{
					double lambda = dot_prod(X[0] - x, y_star - x) / d12;
					if(lambda > 0 && lambda < 1)
					{
						DataVector z = x + lambda * (X[0] - x);
						double d = distancer.sqrDistance(z, y_star);
						if(d < best_d)
						{
							best_x = z;
							best_d = d;
						}
					}
				}

				x_star = best_x;

				DataVector best_y(Y[0]);
				best_d = distancer.sqrDistance(best_y, x_star);

				d12 = distancer.sqrDistance(Y[0], y);
				if(d12 > 0)
				{
					double lambda = dot_prod(Y[0] - y, x_star - y) / d12;
					if(lambda > 0 && lambda < 1)
					{
						DataVector z = y + lambda * (Y[0] - y);
						double d = distancer.sqrDistance(z, x_star);
						if(d < best_d)
						{
							best_y = z;
							best_d = d;
						}
					}
				}

				if(abs(distancer.sqrDistance(x, y) - distancer.sqrDistance(x_star, y_star)) < epsilon) break;
			}
			while(1);

			DataVector w = x_star - y_star;
			double b = (dot_prod(y_star, y_star) - dot_prod(x_star, x_star)) * 0.5;

			// check whether the data is separable (should be removable for our case)
			bool separable = true;
			
			if(dot_prod(X[0], w) + b <= -func_epsilon)
			{
				separable = false;
			}

			if(dot_prod(Y[0], w) + b >= func_epsilon)
			{
				separable = false;
			}

			if(!separable)
			{
				w.setZero();
				b = 0;
			}

			return HyperPlane(w, b, x_star, y_star);
		}
		else
		{
			DataVector x_star(X[rng.uniformInt(0, X.size() - 1)]);
			DataVector y_star(Y[rng.uniformInt(0, Y.size() - 1)]);

			do
			{
				DataVector x(x_star), y(y_star);

				double best_d = 0;

				DataVector best_x(X[0]);
				best_d = distancer.sqrDistance(best_x, y_star);
				for(std::size_t i = 1; i < X.size(); ++i)
				{
					double d = distancer.sqrDistance(X[i], y_star);
					if(d < best_d)
					{
						best_x = X[i];
						best_d = d;
					}
				}

				for(std::size_t i = 0; i < X.size(); ++i)
				{
					double d12 = distancer.sqrDistance(X[i], x);
					if(d12 > 0)
					{
						double lambda = dot_prod(X[i] - x, y_star - x) / d12;
						if(lambda > 0 && lambda < 1)
						{
							DataVector z = x + lambda * (X[i] - x);
							double d = distancer.sqrDistance(z, y_star);
							if(d < best_d)
							{
								best_x = z;
								best_d = d;
							}
						}
					}
				}

				x_star = best_x;

				DataVector best_y(Y[0]);
				best_d = distancer.sqrDistance(best_y, x_star);
				for(std::size_t i = 1; i < Y.size(); ++i)
				{
					double d = distancer.sqrDistance(Y[i], x_star);
					if(d < best_d)
					{
						best_y = Y[i];
						best_d = d;
					}
				}

				for(std::size_t i = 0; i < Y.size(); ++i)
				{
					double d12 = distancer.sqrDistance(Y[i], y);
					if(d12 > 0)
					{
						double lambda = dot_prod(Y[i] - y, x_star - y) / d12;
						if(lambda > 0 && lambda < 1)
						{
							DataVector z = y + lambda * (Y[i] - y);
							double d = distancer.sqrDistance(z, x_star);
							if(d < best_d)
							{
								best_y = z;
								best_d = d;
							}
						}
					}
				}

				if(abs(distancer.sqrDistance(x, y) - distancer.sqrDistance(x_star, y_star)) < epsilon) break;
			}
			while(1);

			DataVector w = x_star - y_star;
			double b = (dot_prod(y_star, y_star) - dot_prod(x_star, x_star)) * 0.5;

			// check whether the data is separable (should be removable for our case)
			bool separable = true;
			for(std::size_t i = 0; i < X.size(); ++i)
			{
				if(dot_prod(X[i], w) + b <= -func_epsilon)
				{
					separable = false;
					break;
				}
			}

			for(std::size_t i = 0; i < Y.size(); ++i)
			{
				if(dot_prod(Y[i], w) + b >= func_epsilon)
				{
					separable = false;
					break;
				}
			}

			if(!separable)
			{
				w.setZero();
				b = 0;
			}

			return HyperPlane(w, b, x_star, y_star);
		}

	}

	MulticonlitronLearner::Conlitron MulticonlitronLearner::SCA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const 
	{
		Conlitron g;
		g.resize(Y.size());
		for(std::size_t i = 0; i < Y.size(); ++i)
			g[i] = CDMA(X, std::vector<DataVector>(1, Y[i]));

		Conlitron f;

		std::vector<std::size_t> Y_old, Y_new;
		Y_old.resize(Y.size());
		for(std::size_t i = 0; i < Y.size(); ++i)
			Y_old[i] = i;

		// cache for g[j].evaluate(Y[j]);
		std::vector<double> evaluate_cache(Y.size());
		for(std::size_t i = 0; i < Y.size(); ++i)
			evaluate_cache[i] = g[i].evaluate(Y[i]);

		do
		{
			std::size_t best_id = -1;
			double best_v = -(std::numeric_limits<double>::max)();
			for(std::size_t i = 0; i < Y_old.size(); ++i)
			{
				std::size_t j = Y_old[i];

				double v = evaluate_cache[j];

				if(v > best_v) 
				{
					best_v = v;
					best_id = j;
				}
			}

			f.push_back(g[best_id]);

			for(std::size_t i = 0; i < Y_old.size(); ++i)
			{
				std::size_t j = Y_old[i];
				if(g[best_id].evaluate(Y[j]) > best_v + func_epsilon2)
					Y_new.push_back(j);
			}

			if(Y_new.size() == 0) break;

			Y_new.swap(Y_old);
			Y_new.clear();
		}
		while(1);

		bool separable = true;

		for(std::size_t i = 0; i < X.size(); ++i)
		{
			for(std::size_t j = 0; j < f.size(); ++j)
			{
				if(f[j].evaluate(X[i]) <= -func_epsilon)
				{
					separable = false;
					goto after_separable_test;
				}
			}
		}

		for(std::size_t i = 0; i < Y.size(); ++i)
		{
			bool all_fail = true;
			for(std::size_t j = 0; j < f.size(); ++j)
			{
				if(f[j].evaluate(Y[i]) <= func_epsilon)
				{
					all_fail = false;
					break;
				}
			}

			if(all_fail)
			{
				separable = false;
				goto after_separable_test;
			}
		}

after_separable_test:

		if(!separable) f.clear();

		return f;

	}

	double MulticonlitronLearner::sqrDistance(const DataVector& query, const std::vector<DataVector>& data) const
	{
		double min_d = (std::numeric_limits<double>::max)();
		for(std::size_t i = 0; i < data.size(); ++i)
		{
			double d = distancer.sqrDistance(query, data[i]);
			if(d < min_d) min_d = d;
		}

		return min_d;
	}

	

	MulticonlitronLearner::MultiConlitron MulticonlitronLearner::SMA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const
	{

		MultiConlitron res;

		std::vector<std::size_t> I_old(X.size()), I_new;
		for(std::size_t i = 0; i < I_old.size(); ++i)
			I_old[i] = i;

		do
		{
			double min_dist = (std::numeric_limits<double>::max)();
			std::size_t min_id = 0;
			for(std::size_t i = 0; i < I_old.size(); ++i)
			{
				std::size_t j = I_old[i];

				double dist;

				if(use_approximate_dist)
				{
					std::vector<int> indices;
					std::vector<double> dists;
					knnSearch<FLANN_WRAPPER::DistanceRN, DataVector>(X[j], feature_dim, knn_index, indices, dists, 1, flann::SearchParams());
					dist = dists[0];
				}
				else
					dist = sqrDistance(X[j], Y);

				if(dist < min_dist)
				{
					min_dist = dist;
					min_id = j;
				}
			}

			std::vector<DataVector> X_;
			X_.push_back(X[min_id]);

			Conlitron CLP = SCA(X_, Y);

			if(CLP.size() == 0)
			{
				res.clear();
				break;
			}

			std::vector<std::pair<double, bool> > CLP_evaluate_cache(CLP.size(), std::make_pair(0, false));

			for(std::size_t i = 0; i < I_old.size(); ++i)
			{
				std::size_t j = I_old[i];
				bool exist = false;
				for(std::size_t k = 0; k < CLP.size(); ++k)
				{
					if(!CLP_evaluate_cache[k].second) CLP_evaluate_cache[k].first = CLP[k].evaluate(X[min_id]);
					if(CLP[k].evaluate(X[j]) < CLP_evaluate_cache[k].first - func_epsilon2)
					{
						exist = true;
						break;
					}
				}

				if(exist)
					I_new.push_back(j);
			}


			if(I_new.size() == 0) break;

			res.push_back(CLP);

			I_new.swap(I_old);
			I_new.clear();

		}while(1);

		return res;
	}


	void MulticonlitronLearner::collectSupportVectors(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp1[k];
				samples.push_back(ContactSpaceSampleData(v, false));
				id++;
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp2[k];
				samples.push_back(ContactSpaceSampleData(v, true));
				id++;
			}
		}	
	}

	void MulticonlitronLearner::collectSupportVectorsClass0(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp1[k];
				samples.push_back(ContactSpaceSampleData(v, false));
				id++;
			}
		}
	}

	void MulticonlitronLearner::collectSupportVectorsClass1(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp2[k];
				samples.push_back(ContactSpaceSampleData(v, true));
				id++;
			}
		}
	}

	void MulticonlitronLearner::collectSupportVectors(std::vector<DataVector>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp1[k];
				samples.push_back(v);
				id++;
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp2[k];
				samples.push_back(v);
				id++;
			}
		}
	}

	void MulticonlitronLearner::collectSupportVectorsClass0(std::vector<DataVector>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp1[k];
				samples.push_back(v);
				id++;
			}
		}
	}

	void MulticonlitronLearner::collectSupportVectorsClass1(std::vector<DataVector>& samples) const
	{
		std::size_t id = 0;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < model.size(); ++i)
		{
			for(std::size_t j = 0; j < model[i].size(); ++j)
			{
				for(std::size_t k = 0; k < feature_dim; ++k)
					v[k] = model[i][j].supp2[k];
				samples.push_back(v);
				id++;
			}
		}
	}

	flann::Index<FLANN_WRAPPER::DistanceRN>* MulticonlitronLearner::constructIndexOfSupportVectors() const
	{
		std::vector<DataVector> samples;
		collectSupportVectors(samples);
		return constructIndexForLearning(samples);
	}

	flann::Index<FLANN_WRAPPER::DistanceRN>* MulticonlitronLearner::constructIndexOfSupportVectorsClass0() const
	{
		std::vector<DataVector> samples;
		collectSupportVectorsClass0(samples);
		return constructIndexForLearning(samples);
	}

	flann::Index<FLANN_WRAPPER::DistanceRN>* MulticonlitronLearner::constructIndexOfSupportVectorsClass1() const
	{
		std::vector<DataVector> samples;
		collectSupportVectorsClass1(samples);
		return constructIndexForLearning(samples);
	}


	PredictResult SVMLearner::predict(const DataVector& query) const
	{
		int svm_type = svm_get_svm_type(model);
		int nr_class = svm_get_nr_class(model);
		double* prob_estimates = NULL;

		svm_node* x = (svm_node *)malloc((feature_dim + 1)*sizeof(svm_node));

		if(param.probability)
		{
			prob_estimates = (double *) malloc(nr_class*sizeof(double));
		}

		const DataVector& v = (scaler && use_scaler) ? scaler->scale(query) : query;

		for(std::size_t j = 0; j < feature_dim; ++j)
		{
			x[j].index = j + 1;
			x[j].value = v[j];
		}
		x[feature_dim].index = -1;

		PredictResult predict_result;
		double predict_label;
		if(param.probability)
		{
			predict_label = svm_predict_probability(model, x, prob_estimates);
			if(predict_label > 0) predict_label = 1;
			else predict_label = 0;
			predict_result = PredictResult(predict_label, *prob_estimates);
		}
		else
		{	
			predict_label = svm_predict(model, x);
			if(predict_label > 0) predict_label = 1;
			else predict_label = 0;
			predict_result = PredictResult(predict_label, 1);
		}

		if(param.probability) free(prob_estimates);
		free(x);

		return predict_result;
	}

	std::vector<PredictResult> SVMLearner::predict(const std::vector<DataVector>& queries) const
	{
		std::vector<PredictResult> predict_results;

		int svm_type = svm_get_svm_type(model);
		int nr_class = svm_get_nr_class(model);
		double* prob_estimates = NULL;

		svm_node* x = (svm_node *)malloc((feature_dim + 1)*sizeof(svm_node));

		if(param.probability)
		{
			prob_estimates = (double *) malloc(nr_class*sizeof(double));
		}

		for(std::size_t i = 0; i < queries.size(); ++i)
		{
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(queries[i]) : queries[i];
			for(std::size_t j = 0; j < feature_dim; ++j)
			{
				x[j].index = j + 1;
				x[j].value = v[j];
			}
			x[feature_dim].index = -1;

			double predict_label;
			if(param.probability)
			{
				predict_label = svm_predict_probability(model, x, prob_estimates);
				if(predict_label > 0) predict_label = 1;
				else predict_label = 0;
				predict_results.push_back(PredictResult(predict_label, *prob_estimates));
			}
			else
			{	
				predict_label = svm_predict(model, x);
				if(predict_label > 0) predict_label = 1;
				else predict_label = 0;
				predict_results.push_back(PredictResult(predict_label, 1));
			}
		}

		if(param.probability) free(prob_estimates);
		free(x);

		return predict_results;
	}


	std::vector<PredictResult> SVMLearner::predict(const std::vector<ContactSpaceSampleData>& queries) const
	{
		std::vector<PredictResult> predict_results;

		int svm_type = svm_get_svm_type(model);
		int nr_class = svm_get_nr_class(model);
		double* prob_estimates = NULL;

		svm_node* x = (svm_node *)malloc((feature_dim + 1)*sizeof(svm_node));

		if(param.probability)
		{
			prob_estimates = (double *) malloc(nr_class*sizeof(double));
		}

		for(std::size_t i = 0; i < queries.size(); ++i)
		{
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(queries[i].v) : queries[i].v;
			for(std::size_t j = 0; j < feature_dim; ++j)
			{
				x[j].index = j + 1;
				x[j].value = v[j];
			}
			x[feature_dim].index = -1;

			double predict_label;
			if(param.probability)
			{
				predict_label = svm_predict_probability(model, x, prob_estimates);
				if(predict_label > 0) predict_label = 1;
				else predict_label = 0;
				predict_results.push_back(PredictResult(predict_label, *prob_estimates));
			}
			else
			{	
				predict_label = svm_predict(model, x);
				if(predict_label > 0) predict_label = 1;
				else predict_label = 0;
				predict_results.push_back(PredictResult(predict_label, 1));
			}
		}

		if(param.probability) free(prob_estimates);
		free(x);

		return predict_results;
	}

	void SVMLearner::learn(const std::vector<ContactSpaceSampleData>& data, const std::vector<double>& weights, std::size_t active_dim)
	{
		if(data.size() == 0) return;
		
		if(model) svm_free_and_destroy_model(&model);
		
		feature_dim = active_dim;
		
		if(param.gamma == 0) param.gamma = 1.0 / feature_dim;
		
		// set svm_problem
		problem.l = data.size();
		if(problem.y) delete [] problem.y;
		problem.y = new double[problem.l];
		if(problem.x) delete [] problem.x;
		problem.x = new svm_node* [problem.l];
		if(problem.W) delete [] problem.W;
		problem.W = new double[problem.l];
		if(x_space) delete [] x_space;
		x_space = new svm_node[(feature_dim + 1)* problem.l];
		
		for(std::size_t i = 0; i < data.size(); ++i)
		{
			svm_node* cur_x_space = x_space + (feature_dim + 1) * i;
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(data[i].v) : data[i].v;
			for(std::size_t j = 0; j < feature_dim; ++j)
			{
				cur_x_space[j].index = j + 1;
				cur_x_space[j].value = v[j];
			}
			cur_x_space[feature_dim].index = -1;
			
			problem.x[i] = cur_x_space;
			problem.y[i] = (data[i].col ? 1 : -1);
			problem.W[i] = weights[i];
		}
		
		// build model & classify
		model = svm_train(&problem, &param);
		hyperw_normsqr = svm_hyper_w_normsqr_twoclass(model);
	}

	void SVMLearner::learn(const std::vector<ContactSpaceSampleData>& data, std::size_t active_dim)
	{
		if(data.size() == 0) return;
		
		if(model) svm_free_and_destroy_model(&model);
		
		feature_dim = active_dim;
		
		if(param.gamma == 0) param.gamma = 1.0 / feature_dim;
		
		// set svm_problem
		problem.l = data.size();
		if(problem.y) delete [] problem.y;
		problem.y = new double[problem.l];
		if(problem.x) delete [] problem.x;
		problem.x = new svm_node* [problem.l];
		if(problem.W) delete [] problem.W;
		problem.W = new double[problem.l];
		if(x_space) delete [] x_space;
		x_space = new svm_node[(feature_dim + 1)* problem.l];
		
		for(std::size_t i = 0; i < data.size(); ++i)
		{
			svm_node* cur_x_space = x_space + (feature_dim + 1) * i;
			const DataVector& v = (scaler && use_scaler) ? scaler->scale(data[i].v) : data[i].v;
			for(std::size_t j = 0; j < feature_dim; ++j)
			{
				cur_x_space[j].index = j + 1;
				cur_x_space[j].value = v[j];
			}
			cur_x_space[feature_dim].index = -1;
			
			problem.x[i] = cur_x_space;
			problem.y[i] = (data[i].col ? 1 : -1);
			problem.W[i] = 1;
		}
		
		// build model & classify
		model = svm_train(&problem, &param);
		hyperw_normsqr = svm_hyper_w_normsqr_twoclass(model);
	}

	void SVMLearner::incremental_learn(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
	{
		if(!model || data_.size() == 0) return;

		std::size_t dim = data_[0].v.dim();

		std::vector<ContactSpaceSampleData> data = data_;
	
		for(std::size_t i = 0; i < model->l; ++i)
		{
			DataVector v_(active_dim);
			for(std::size_t j = 0; j < active_dim; ++j)
				v_[j] = model->SV[i][j].value;

			const DataVector& scaled_supp = (scaler && use_scaler) ? scaler->unscale(v_) : v_;

			DataVector v(dim);
			for(std::size_t j = 0; j < active_dim; ++j)
				v[j] = scaled_supp[j];

			int id = model->sv_indices[i] - 1;			
			data.push_back(ContactSpaceSampleData(v, (problem.y[id] > 0))); // lucky, scaled_supp works (in fact it should be extended).
		}

		learn(data, active_dim);
	}

	HyperPlane SVMLearner::getLinearModel() const
	{
		double b = -*(model->rho);
		DataVector w(feature_dim);

		std::size_t numSV = model->l;

		for(std::size_t i = 0; i < feature_dim; ++i)
		{
			w[i] = 0;
			for(std::size_t j = 0; j < numSV; ++j)
				w[i] += model->sv_coef[0][j] * model->SV[j][i].value;
		}

		if(model->label[0] == -1)
		{
			for(std::size_t i = 0; i < feature_dim; ++i)
				w[i] = -w[i];
			b = -b;
		}

		return HyperPlane(w, b);
	}


	void SVMLearner::collectSupportVectors(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			int id = model->sv_indices[i] - 1;
			bool col = (problem.y[id] > 0);
			samples.push_back(ContactSpaceSampleData(v, col));
		}
	}

	void SVMLearner::collectSupportVectorsClass0(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			int id = model->sv_indices[i] - 1;
			bool col = (problem.y[id] > 0);
			if(!col)
				samples.push_back(ContactSpaceSampleData(v, col));
		}
	}

	void SVMLearner::collectSupportVectorsClass1(std::vector<ContactSpaceSampleData>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			int id = model->sv_indices[i] - 1;
			bool col = (problem.y[id] > 0);
			if(col)
				samples.push_back(ContactSpaceSampleData(v, col));
		}
	}

	void SVMLearner::collectSupportVectors(std::vector<DataVector>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			samples.push_back(v);
		}
	}

	void SVMLearner::collectSupportVectorsClass0(std::vector<DataVector>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			int id = model->sv_indices[i] - 1;
			bool col = (problem.y[id] > 0);
			if(!col)
				samples.push_back(v);
		}
	}

	void SVMLearner::collectSupportVectorsClass1(std::vector<DataVector>& samples) const
	{
		std::size_t data_num = model->l;
		DataVector v(feature_dim);
		for(std::size_t i = 0; i < data_num; ++i)
		{
			for(std::size_t j = 0; j < feature_dim; ++j)
				v[j] = model->SV[i][j].value;
			int id = model->sv_indices[i] - 1;
			bool col = (problem.y[id] > 0);
			if(col)
				samples.push_back(v);
		}
	}

	flann::Index<FLANN_WRAPPER::DistanceRN>* SVMLearner::constructIndexOfSupportVectors() const
	{
		std::vector<DataVector> samples;
		collectSupportVectors(samples);
		return constructIndexForLearning(samples);
	}

	flann::Index<FLANN_WRAPPER::DistanceRN>* SVMLearner::constructIndexOfSupportVectorsClass0() const
	{
		std::vector<DataVector> samples;
		collectSupportVectorsClass0(samples);
		return constructIndexForLearning(samples);
	}
	
	flann::Index<FLANN_WRAPPER::DistanceRN>* SVMLearner::constructIndexOfSupportVectorsClass1() const
	{
		std::vector<DataVector> samples;
		collectSupportVectorsClass1(samples);
		return constructIndexForLearning(samples);
	}
}




