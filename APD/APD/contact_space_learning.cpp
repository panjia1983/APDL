#include "contact_space_learning.h"
#include "distance_proxy.h"

namespace APDL
{
	std::vector<PredictResult> MulticonlitronLearner::predict(const std::vector<ContactSpaceSampleData>& queries) const
	{
		std::vector<PredictResult> predict_results;

		for(std::size_t i = 0; i < queries.size(); ++i)
		{
			const DataVector& v = queries[i].v;

			double predict_label = model.evaluate2(v);
			if(predict_label > 0) predict_label = 0;
			else predict_label = 1;
			
			predict_results.push_back(PredictResult(predict_label, 1));
		}

		return predict_results;
	}

	HyperPlane MulticonlitronLearner::CDMA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y, double epsilon) const
	{
		RNG rng;
		DataVector x_star = X[rng.uniformInt(0, X.size() - 1)];
		DataVector y_star = Y[rng.uniformInt(0, Y.size() - 1)];


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

			if(std::abs(distancer.sqrDistance(x, y) - distancer.sqrDistance(x_star, y_star)) < epsilon) break;
		}
		while(1);

		DataVector w = x_star - y_star;
		double b = (dot_prod(y_star, y_star) - dot_prod(x_star, x_star)) * 0.5;

		// check whether the data is separable (should be removable for our case)
		bool separable = true;
		for(std::size_t i = 0; i < X.size(); ++i)
		{
			if(dot_prod(X[i], w) + b < 0)
			{
				separable = false;
				break;
			}
		}

		for(std::size_t i = 0; i < Y.size(); ++i)
		{
			if(dot_prod(Y[i], w) + b > 0)
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

	MulticonlitronLearner::Conlitron MulticonlitronLearner::SCA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const 
	{
		Conlitron g;
		for(std::size_t i = 0; i < Y.size(); ++i)
		{
			std::vector<DataVector> y;
			y.push_back(Y[i]);
			g.push_back(CDMA(X, y));
		}

		Conlitron f;

		std::vector<std::size_t> Y_old, Y_new;
		Y_old.resize(Y.size());
		for(std::size_t i = 0; i < Y.size(); ++i)
			Y_old[i] = i;

		do
		{
			std::size_t best_id = -1;
			double best_v = -(std::numeric_limits<double>::max)();
			for(std::size_t i = 0; i < Y_old.size(); ++i)
			{
				std::size_t j = Y_old[i];
				double v = g[j].evaluate(Y[j]);
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
				if(g[best_id].evaluate(Y[j]) > best_v)
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
				if(f[j].evaluate(X[i]) < 0)
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
				if(f[j].evaluate(Y[i]) <= 0)
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
				double dist = sqrDistance(X[j], Y);
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

			for(std::size_t i = 0; i < I_old.size(); ++i)
			{
				std::size_t j = I_old[i];
				bool exist = false;
				for(std::size_t k = 0; k < CLP.size(); ++k)
				{
					if(CLP[k].evaluate(X[j]) < CLP[k].evaluate(X[min_id]))
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

	//std::vector<PredictResult> SVMLearner::predict2(const std::vector<ContactSpaceSampleData>& queries, std::size_t active_dim) const
	//{
	//	std::vector<PredictResult> predict_results;
	//	int svm_type = svm_get_svm_type(model);
	//	int nr_class = svm_get_nr_class(model);

	//	assert((svm_type == C_SVC) || (svm_type == NU_SVC));
	//	assert(nr_class == 2);

	//	svm_node* x = (svm_node *)malloc((active_dim + 1)*sizeof(svm_node));

	//	flann::Index<FLANN_WRAPPER::DistanceRN>* index;
	//	std::vector<std::vector<double> > boundary;
	//	if(scaler)
	//	{
	//		double* lower = new double[active_dim];
	//		double* upper = new double[active_dim];
	//		for(int j = 0; j < active_dim; ++j)
	//		{
	//			lower[j] = use_scaler ? 0 : scaler->v_min[j];
	//			upper[j] = use_scaler ? 1 : scaler->v_max[j];
	//		}
	//		sample_decision_boundary_brute_force(model, upper, lower, active_dim, 0.005, 100, boundary, index);
	//	}

	//	std::ofstream file("collect.txt");

	//	double inv_hyperw = 1 / sqrt(hyperw_normsqr);
	//	// for(std::size_t i = 0; i < queries.size(); ++i)
	//	for(std::size_t i = 0; i < queries.size(); ++i)
	//	{
	//		const DataVector& v = (scaler && use_scaler) ? scaler->scale(queries[i].v) : queries[i].v;
	//		for(std::size_t j = 0; j < active_dim; ++j)
	//		{
	//			x[j].index = j + 1;
	//			x[j].value = v[j];
	//		}
	//		x[active_dim].index = -1;

	//		double predict_label = svm_predict_values_twoclass(model, x);

	//		double s = (scaler && use_scaler) ? scaler->getScale() : 1;
	//		double f_dist = predict_label * inv_hyperw;
	//		double i_dist = sqrt(-1/param.gamma * log(1 - f_dist * f_dist * 0.5)) * s;


	//		double real_i_dist1 = dist_to_decision_boundary(model, x) * s;
	//		double real_i_dist2 = dist_to_decision_boundary_with_gradient(model, x) * s;
	//		double real_i_dist3 = dist_to_decision_boundary_constrain_free(model, sqrt(hyperw_normsqr), x) * s;

	//		bool is_sv = false;
	//		for(int j = 0; j < model->l; ++j)
	//		{
	//			if(i == model->sv_indices[j])
	//			{
	//				is_sv = true;
	//				break;
	//			}
	//		}

	//		double brute_force_dist = dist_to_decision_boundary_brute_force(boundary, index, x) * s;
	//		
	//		file << is_sv << " " << brute_force_dist << " 1) " << f_dist << " 2) " << i_dist << " 3) " << real_i_dist1 << " 4) " << real_i_dist2 << " 5) " << real_i_dist3 << std::endl;

	//		if(predict_label > 0) predict_label = 1;
	//		else predict_label = 0;

	//		predict_results.push_back(PredictResult(predict_label, 1));
	//	}
	//	std::cout << std::endl;
	//	file.close();


	//	// interpolate and check
	//	std::vector<std::vector<double> > samples;
	//	sample_decision_boundary_interpolation(model, active_dim, sqrt(hyperw_normsqr), samples);

	//	svm_node* sample = new svm_node[active_dim + 1];
	//	for(std::size_t i = 0; i < active_dim; ++i) sample[i].index = i + 1;
	//	sample[active_dim].index = -1;

	//	std::cout << samples.size() << std::endl;
	//	for(std::size_t i = 0; i < samples.size(); ++i)
	//	{
	//		for(std::size_t j = 0; j < active_dim; ++j) sample[j].value = samples[i][j];
	//		
	//		double val = dist_to_decision_boundary_brute_force(boundary, index, sample);
	//	}

	//	delete [] sample;

	//	free(x);
	//	delete index;


	//	return predict_results;
	//}


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

}