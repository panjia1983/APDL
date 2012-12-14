#ifndef CONTACT_SPACE_LEARNING
#define CONTACT_SPACE_LEARNING


#include "contact_space.h"
#include "distance_proxy.h"
#include <limits>
#include "learning/svm.h"
#include "learning/incsvm.h"



namespace APDL
{
					
	struct HyperPlane
	{
		DataVector w;
		double b;

		DataVector supp1, supp2;

		HyperPlane(const DataVector& w_, double b_) : w(w_), b(b_) 
		{
			double inv_w_len = 1.0 / sqrt(dot_prod(w, w));
			w = inv_w_len * w;
			b = inv_w_len * b;
		}

		HyperPlane(const DataVector& w_, double b_, const DataVector& supp1_, const DataVector& supp2_) : w(w_), b(b_), supp1(supp1_), supp2(supp2_)
		{
			double inv_w_len = 1.0 / sqrt(dot_prod(w, w));
			w = inv_w_len * w;
			b = inv_w_len * b;
		}

		inline double evaluate(const DataVector& v) const
		{
			return dot_prod(w, v) + b;
		}

		double distance(const DataVector& v) const
		{
			return abs(evaluate(v));
		}

		DataVector project(const DataVector& v) const
		{
			return v - evaluate(v) * w;
		}

		std::size_t dim() const { return w.dim(); }


	};

	struct WeightedHyperPlane
	{
		HyperPlane hp;
		double a;

		WeightedHyperPlane(const HyperPlane& hp_, double a_) : hp(hp_), a(a_) {}

		double evaluate(const DataVector& v) const
		{
			if(hp.evaluate(v) > 0)
				return a;
			else return -a;
		}
	};

	struct MultipleWeightedHyperPlane
	{
		std::vector<WeightedHyperPlane> whps;

		double evaluate(const DataVector& v) const
		{
			double value = 0;
			for(std::size_t i = 0; i < whps.size(); ++i)
				value += whps[i].evaluate(v);
			return value;
		}
	};

	struct PredictResult
	{
		int label;
		double prob;

		PredictResult(int label_, double prob_)
		{
			label = label_;
			prob = prob_;
		}
	};

	class SVMLearner
	{
	public:


		SVMLearner()
		{
			param.svm_type = C_SVC;
			param.kernel_type = RBF;
			param.degree = 3;
			param.gamma = 0;	// 1/num_features
			param.coef0 = 0;
			param.nu = 0.5;
			param.cache_size = 100; // can change
			param.C = 1;
			param.eps = 1e-3;
			param.p = 0.1;
			param.shrinking = 1;    // use shrinking
			param.probability = 0;
			param.nr_weight = 0;
			param.weight_label = NULL;
			param.weight = NULL;
			
			param.nr_weight = 2;
			param.weight_label = (int *)realloc(param.weight_label, sizeof(int) * param.nr_weight);
			param.weight = (double *)realloc(param.weight, sizeof(double) * param.nr_weight);
			param.weight_label[0] = -1;
			param.weight_label[1] = 1;
			param.weight[0] = 1;
			param.weight[1] = 1;
			
			model = NULL;
			x_space = NULL;
			problem.x = NULL;
			problem.y = NULL;
			problem.W = NULL;

			scaler = NULL;
			use_scaler = false;
		}

		void setCSVM() { param.svm_type = C_SVC; }
		void setNuSVM() { param.svm_type = NU_SVC; }
		void setC(double C) { param.C = C; }
		void setNu(double nu) { param.nu = nu; }
		void setLinearClassifier() { param.kernel_type = LINEAR; }
		void setNonLinearClassifier() { param.kernel_type = RBF; }
		void setProbability(bool use_probability) { param.probability = use_probability; }
		void setScaler(const Scaler& scaler_)
		{
			if(scaler) delete scaler;
			scaler = new Scaler(scaler_.v_min, scaler_.v_max, scaler_.active_dim);
		}

		void setNegativeWeight(double c)
		{
			param.weight[0] = c;
		}

		void setPositiveWeight(double c)
		{
			param.weight[1] = c;
		}

		void setEPS(double e)
		{
			param.eps = e;
		}

		void setUseScaler(bool use_)
		{
			use_scaler = use_;
		}

		void setGamma(double gamma)
		{
			param.gamma = gamma;
		}
		
		~SVMLearner()
		{
			svm_destroy_param(&param); 
			svm_free_and_destroy_model(&model);
			delete [] x_space;
			delete [] problem.x;
			delete [] problem.y;
			delete [] problem.W;
			delete scaler;
		}
		
		void learn(const std::vector<ContactSpaceSampleData>& data, std::size_t active_dim);

		void learn(const std::vector<ContactSpaceSampleData>& data, const std::vector<double>& weights, std::size_t active_dim);

		// std::vector<PredictResult> predict2(const std::vector<ContactSpaceSampleData>& queries, std::size_t active_dim) const;

		std::vector<PredictResult> predict(const std::vector<ContactSpaceSampleData>& queries) const;
		
		void save(const std::string& file_name) const
		{
			if(model)
				svm_save_model(file_name.c_str(), model);
		}

		HyperPlane getLinearModel() const;


		void saveVisualizeData(const std::string& file_name, const Scaler& scaler, std::size_t grid_n) const
		{
			std::size_t dim = feature_dim;
			std::ofstream os(file_name.c_str());

			DataVector upper = scaler.v_max;
			DataVector lower = scaler.v_min;

			DataVector scale = (upper - lower) * (1.0 / grid_n);
			upper = upper + scale * 10;
			lower = lower - scale * 10;

			grid_n += 20;

			std::size_t n = 1;
			for(std::size_t i = 0; i < dim; ++i)
				n *= grid_n;

			std::size_t* ids = new std::size_t[dim];

			svm_node* v = new svm_node[dim + 1];
			for(std::size_t i = 0; i < dim; ++i) v[i].index = i + 1;
			v[dim].index = -1;

			for(std::size_t i = 0; i < n; ++i)
			{
				std::size_t i_ = i;
				for(std::size_t j = 0; j < dim; ++j)
				{
					ids[j] = i_ % grid_n;
					i_ = i_ / grid_n;
				}

				for(std::size_t j = 0; j < dim; ++j)
					v[j].value = (upper[j] - lower[j]) * ids[j] / (double)grid_n + lower[j];

				double eval = svm_predict_values_twoclass(model, v);

				//for(std::size_t j = 0; j < dim; ++j)
				//	os << v[j].value << " ";
				// os << eval << endl;

				os << eval << " ";
				if((i % grid_n) == (grid_n - 1))
					os << std::endl;
			}

			os.close();

			delete [] v;
		}

		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectors() const;
		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectorsClass0() const;
		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectorsClass1() const;

		svm_parameter param;
		
		svm_problem problem;
		svm_node* x_space;
		
		svm_model* model;
		double hyperw_normsqr;

		std::size_t feature_dim;

		Scaler* scaler;

		bool use_scaler;
	};

	

	class MulticonlitronLearner
	{
	public:

		struct Conlitron : public std::vector<HyperPlane>
		{
			double evaluate2(const DataVector& v) const
			{
				double min_v = (std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < size(); ++i)
				{
					double eval = this->operator[](i).evaluate(v);
					if(eval < min_v) min_v = eval;
				}

				return min_v;
			}

			int evaluate(const DataVector& v) const
			{
				for(std::size_t i = 0; i < size(); ++i)
				{
					if(this->operator[](i).evaluate(v) < 0) return -1;
				}

				return 1;
			}

			std::size_t numOfHyperPlanes() const { return size(); }
		};

		struct MultiConlitron : public std::vector<Conlitron>
		{
			double evaluate2(const DataVector& v) const
			{
				double max_v = -(std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < size(); ++i)
				{
					double eval = this->operator[](i).evaluate2(v);
					if(eval > max_v) max_v = eval;
				}

				return max_v;
			}

			int evaluate(const DataVector& v) const
			{
				for(std::size_t i = 0; i < size(); ++i)
				{
					if(this->operator[](i).evaluate2(v) > 0) return 1;
				}
				return -1;
			}

			std::size_t numOfConlitrons() const { return size(); }
			std::size_t numOfHyperPlanes() const
			{
				std::size_t n = 0;
				for(std::size_t i = 0; i < size(); ++i)
					n += this->operator[](i).numOfHyperPlanes();

				return n;
			}
		};

		MulticonlitronLearner(const DataVector& weight, double epsilon_) : distancer(weight), epsilon(epsilon_)
		{
			feature_dim = weight.dim();
			scaler = NULL;

		}

		HyperPlane CDMA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y, double epsilon = 0.01) const;

		Conlitron SCA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const;

		double sqrDistance(const DataVector& query, const std::vector<DataVector>& data) const;

		MultiConlitron SMA(const std::vector<DataVector>& X, const std::vector<DataVector>& Y) const;

		std::vector<PredictResult> predict(const std::vector<ContactSpaceSampleData>& queries) const;

		void learn(const std::vector<ContactSpaceSampleData>& data, std::size_t active_dim)
		{
			if(active_dim  != feature_dim)
			{
				std::cout << "The dimension of the learner is not correct!" << std::endl;
				return;
			}

			std::vector<DataVector> X, Y;
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				DataVector v(distancer.dim());
				for(std::size_t j = 0; j < distancer.dim(); ++j)
					v[j] = data[i].v[j];

				if(data[i].col == false) 
					X.push_back(v);
				else
					Y.push_back(v);
			}

			model = SMA(X, Y);
		}

		void setScaler(const Scaler& scaler_)
		{
			if(scaler) delete scaler;
			scaler = new Scaler(scaler_.v_min, scaler_.v_max, scaler_.active_dim);
		}

		void saveVisualizeData(const std::string& file_name, const Scaler& scaler, std::size_t grid_n) const
		{
			std::size_t dim = model[0][0].w.dim();
			std::ofstream os(file_name.c_str());

			DataVector upper = scaler.v_max;
			DataVector lower = scaler.v_min;

			DataVector scale = (upper - lower) * (1.0 / grid_n);
			upper = upper + scale * 10;
			lower = lower - scale * 10;

			grid_n += 20;

			std::size_t n = 1;
			for(std::size_t i = 0; i < dim; ++i)
				n *= grid_n;

			std::size_t* ids = new std::size_t[dim];
			DataVector v(dim);
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
					v[j] = (upper[j] - lower[j]) * ids[j] / (double)grid_n + lower[j];
				}

				double eval = model.evaluate2(v);

				for(std::size_t j = 0; j < dim; ++j)
					os << v[j] << " ";
				os << eval << std::endl;	
			}

			os.close();
		}

		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectors() const;
		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectorsClass0() const;
		flann::Index<FLANN_WRAPPER::DistanceRN>* constructIndexOfSupportVectorsClass1() const;

		DistanceProxyRN distancer;

		MultiConlitron model;

		std::size_t feature_dim;

		double epsilon;

		Scaler* scaler;
	};

	class AdaBoostLearner
	{
	public:
		typedef FLANN_WRAPPER::DistanceRN ClassifierDistancer;
		AdaBoostLearner()
		{
			knn_index = NULL;
		}

		~AdaBoostLearner()
		{
			delete knn_index;
		}

		std::vector<PredictResult> predict(const std::vector<ContactSpaceSampleData>& queries) const
		{
			std::vector<PredictResult> results;
			for(std::size_t i = 0; i < queries.size(); ++i)
			{
				DataVector v(feature_dim);
				for(std::size_t j = 0; j < feature_dim; ++j)
					v[j] = queries[i].v[j];

				double value = model.evaluate(v);

				int pred = (value > 0)? 1 : 0;

				results.push_back(PredictResult(pred, 1));
			}
			return results;
		}

		std::vector<ContactSpaceSampleData> balanceSample(const std::vector<ContactSpaceSampleData>& data_) const
		{
			std::vector<ContactSpaceSampleData> sampled_data;
			int n_pos = 0;
			int n_neg = 0;
			for(std::size_t i = 0; i < data_.size(); ++i)
			{
				if(data_[i].col == 0) n_neg++;
				else n_pos++;
			}

			if(n_pos > n_neg)
			{
				double ratio = n_neg / (double)n_pos;
				for(std::size_t i = 0; i < data_.size(); ++i)
				{
					if(!data_[i].col) sampled_data.push_back(data_[i]);
					else
					{
						double t = rng.uniformReal(0, 1);
						if(t < ratio)
							sampled_data.push_back(data_[i]);
					}
				}
			}
			else
			{
				double ratio = n_pos / (double)n_neg;
				for(std::size_t i = 0; i < data_.size(); ++i)
				{
					if(data_[i].col) sampled_data.push_back(data_[i]);
					else
					{
						double t = rng.uniformReal(0, 1);
						if(t < ratio)
							sampled_data.push_back(data_[i]);
					}
				}
			}

			return sampled_data;
		}

		void learn(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
		{
			// data = balanceSample(data_);
			data = data_;

			if(knn_index) delete knn_index;
			generateIndex<ContactSpaceR2::DistanceType>(data, 
				active_dim, 
				knn_index,
				flann::KDTreeIndexParams());

			positive_indices.clear();
			negative_indices.clear();
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				if(data[i].col) positive_indices.push_back(i);
				else negative_indices.push_back(i);
			}

			weak_classifier_id = 0;
			
			learnC1(data, active_dim);

			std::vector<PredictResult> results = predict(data);

			std::size_t error_n = 0;
			for(std::size_t i = 0; i < results.size(); ++i)
			{
				error_n += (results[i].label != data[i].col);
				std::cout << "(" << results[i].label << "," << data[i].col << ")";
			}
			std::cout << std::endl;
			std::cout << error_n / (double)data.size() << std::endl;
			
		}

		void learnC0(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
		{
			data = data_;
			feature_dim = active_dim;
			double initial_weight = 1;
			weights.resize(data.size(), initial_weight);
			weights_sum = data.size() * initial_weight;

			std::vector<WeightedHyperPlane> whps;
			std::vector<bool> preds(data.size());

			for(std::size_t iter = 0; iter < 100; ++iter)
			{
				HyperPlane hp = weakClassifier3();

				for(std::size_t i = 0; i < data.size(); ++i)
				{
					DataVector v(feature_dim);
					for(std::size_t j = 0; j < feature_dim; ++j)
						v[j] = data[i].v[j];
					if(hp.evaluate(v) > 0)
						preds[i] = true;
					else preds[i] = false;
				}

				double error = 0;
				for(std::size_t i = 0; i < data.size(); ++i)
				{
					if(preds[i] != data[i].col)
						error += weights[i] * 1;
				}
				error /= weights_sum;

				if(error > 0.5)
				{
					error = 1 - error;
					hp.w = -hp.w;
					hp.b = -hp.b;
				}

				std::cout << error << " ";

				if((0.5 - error < 0.001) || error < 0.001) 
				{
					if(whps.empty())
						whps.push_back(WeightedHyperPlane(hp, 1));
					break;
				}


				double alpha = 0.5 * log((1-error)/error);

				std::cout << "alpha " << alpha << std::endl;

				whps.push_back(WeightedHyperPlane(hp, alpha));

				double weights_max = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
				{
					weights[i] *= exp(-alpha * (2 * preds[i] - 1) * (2 * data[i].col - 1));
					if(weights[i] > weights_max) weights_max = weights[i];
				}

				weights_max /= initial_weight;

				for(std::size_t i = 0; i < weights.size(); ++i)
					weights[i] /= weights_max;

				weights_sum = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
					weights_sum += weights[i];
			}

			model.whps = whps;
		}

		std::pair<double, double> getCost(const std::vector<ContactSpaceSampleData>& data) const
		{
			std::size_t n_pos, n_neg = 0;
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				if(data[i].col == 0) n_neg++;
			}

			n_pos = data.size() - n_neg;

			double C_pos, C_neg;
			if(n_pos > n_neg)
			{
				C_neg = 1;
				C_pos = n_neg / (double)n_pos * 0.8;
			}
			else
			{
				C_pos = 1;
				C_neg = n_pos / (double)n_neg * 0.8;
			}

			return std::make_pair(C_pos, C_neg);
		}

		void learnC1(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
		{
			data = data_;
			feature_dim = active_dim;
			double initial_weight = 1.0;
			weights.resize(data.size(), initial_weight);
			weights_sum = data.size() * initial_weight;

			std::pair<double, double> costs = getCost(data);
			double C_pos = costs.first;
			double C_neg = costs.second;

			std::cout << C_pos << " " << C_neg << std::endl;

			std::vector<WeightedHyperPlane> whps;
			std::vector<bool> preds(data.size());

			for(std::size_t iter = 0; iter < 100; ++iter)
			{
				HyperPlane hp = weakClassifier3();

				for(std::size_t i = 0; i < data.size(); ++i)
				{
					DataVector v(feature_dim);
					for(std::size_t j = 0; j < feature_dim; ++j)
						v[j] = data[i].v[j];
					if(hp.evaluate(v) > 0)
						preds[i] = true;
					else preds[i] = false;
				}

				double delta = 0;
				for(std::size_t i = 0; i < data.size(); ++i)
				{
					if(preds[i] != data[i].col)
					{
						if(data[i].col == 0)
							delta -= weights[i] * C_neg;
						else
							delta -= weights[i] * C_pos;
					}
					else
					{
						if(data[i].col == 0)
							delta += weights[i] * C_neg;
						else
							delta += weights[i] * C_pos;
					}
				}

				delta /= weights_sum;

				if(delta < 0)
				{
					delta = -delta;
					hp.w = -hp.w;
					hp.b = -hp.b;
				}

				std::cout << delta << std::endl;

				if(delta < 0.001 || delta > 1 - 0.001) break;


				double alpha = 0.5 * log((1+delta)/(1-delta));
				std::cout << "delta " << delta << " alpha " << alpha << std::endl;

				whps.push_back(WeightedHyperPlane(hp, alpha));

				double weights_max = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
				{
					double C = (data[i].col == 0) ? C_neg : C_pos;
					weights[i] *= exp(-alpha * C * (2 * preds[i] - 1) * (2 * data[i].col - 1));
					if(weights[i] > weights_max) weights_max = weights[i];
				}

				for(std::size_t i = 0; i < weights.size(); ++i)
					weights[i] /= weights_max;

				weights_sum = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
					weights_sum += weights[i];
			}

			model.whps = whps;
		}


		void learnC2(const std::vector<ContactSpaceSampleData>& data_, std::size_t active_dim)
		{
			data = data_;
			feature_dim = active_dim;
			double initial_weight = 1.0;
			weights.resize(data.size(), initial_weight);
			weights_sum = data.size() * initial_weight;

			std::pair<double, double> costs = getCost(data);
			double C_pos = costs.first;
			double C_neg = costs.second;

			std::cout << C_pos << " " << C_neg << std::endl;

			std::vector<WeightedHyperPlane> whps;
			std::vector<bool> preds(data.size());

			for(std::size_t iter = 0; iter < 100; ++iter)
			{
				HyperPlane hp = weakClassifier();

				for(std::size_t i = 0; i < data.size(); ++i)
				{
					DataVector v(feature_dim);
					for(std::size_t j = 0; j < feature_dim; ++j)
						v[j] = data[i].v[j];
					if(hp.evaluate(v) > 0)
						preds[i] = true;
					else preds[i] = false;
				}

				double delta1 = 0, delta2 = 0;
				for(std::size_t i = 0; i < data.size(); ++i)
				{
					if(preds[i] != data[i].col)
					{
						if(data[i].col == 0)
							delta1 += weights[i] * C_neg;
						else
							delta1 += weights[i] * C_pos;
					}
					else
					{
						if(data[i].col == 0)
							delta2 += weights[i] * C_neg;
						else
							delta2 += weights[i] * C_pos;
					}
				}

				delta2 /= weights_sum;
				delta1 /= weights_sum;

				if(delta2 < delta1)
				{
					std::swap(delta1, delta2);
					hp.w = -hp.w;
					hp.b = -hp.b;
				}

				std::cout << delta2 << " " << delta1 << std::endl;

				if(delta2 - delta1 < 0.001 || delta1 < 0.001) break;


				double alpha = 0.5 * log(delta2/delta1);

				whps.push_back(WeightedHyperPlane(hp, alpha));

				double weights_max = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
				{
					double C = (data[i].col == 0) ? C_neg : C_pos;
					weights[i] *= (C * exp(-alpha * (2 * preds[i] - 1) * (2 * data[i].col - 1)));
					if(weights[i] > weights_max) weights_max = weights[i];
				}

				for(std::size_t i = 0; i < weights.size(); ++i)
					weights[i] /= weights_max;

				weights_sum = 0;
				for(std::size_t i = 0; i < weights.size(); ++i)
					weights_sum += weights[i];
			}

			model.whps = whps;
		}

		HyperPlane weakClassifier()
		{
			SVMLearner learner;
			learner.setLinearClassifier();

			learner.learn(data, weights, feature_dim);

			return learner.getLinearModel();
		}

		HyperPlane weakClassifier2()
		{	
			do
			{
				if(weak_classifier_id >= negative_indices.size()) weak_classifier_id = 0;

				const ContactSpaceSampleData& query = data[negative_indices[weak_classifier_id]];

				std::vector<int> knn_indices;
				std::vector<double> knn_dists;

				knnSearch<ClassifierDistancer>(query, feature_dim, knn_index, knn_indices, knn_dists, 20, flann::SearchParams());


				SVMLearner learner;
				learner.setLinearClassifier();

				std::vector<ContactSpaceSampleData> sub_data;
				std::vector<double> sub_weights(knn_indices.size());
				bool all_one_type = true;
				for(std::size_t i = 0; i < knn_indices.size(); ++i)
				{
					sub_data.push_back(data[knn_indices[i]]);
					sub_weights[i] = weights[knn_indices[i]];
					if(data[knn_indices[i]].col != data[knn_indices[0]].col) all_one_type = false;
				}

				if(all_one_type)
				{
					weak_classifier_id++;
					continue;
				}

				learner.learn(sub_data, sub_weights, feature_dim);

				weak_classifier_id++;

				return learner.getLinearModel();
			}while(1);
		}

		struct DimensionSorter
		{
			DimensionSorter(const std::vector<ContactSpaceSampleData>& data_, std::size_t dim_id_) : data(data_), dim_id(dim_id_)
			{

			}

			bool operator () (std::size_t i, std::size_t j) const
			{
				return data[i].v[dim_id] < data[j].v[dim_id];
			}

			const std::vector<ContactSpaceSampleData>& data;
			std::size_t dim_id;
		};

		HyperPlane weakClassifier3()
		{
			if(weak_classifier_id >= feature_dim) weak_classifier_id = 0;

			std::vector<std::size_t> v(data.size());
			for(std::size_t i = 0; i < data.size(); ++i)
				v[i] = i;

			DimensionSorter comp(data, weak_classifier_id);
			std::sort(v.begin(), v.end(), comp);

			double error = 0;
			for(std::size_t i = 0; i < data.size(); ++i)
			{
				if(!data[i].col) error += weights[i] * 1;
			}	

			double min_error = error;
			std::cout << "error0 " << min_error << std::endl;
			std::size_t min_id = -1;
			for(std::size_t i = 0; i < v.size(); ++i)
			{
				std::size_t id = v[i];
				double weight = weights[id];
				if(data[id].col) error += weight;
				else error -= weight;
				if(error < min_error) { min_error = error; min_id = i; }
			}

			DataVector w(feature_dim);
			for(std::size_t i = 0; i < feature_dim; ++i)
				w[i] = 0;
			w[weak_classifier_id] = 1;
			double b;
			if(min_id == (std::size_t)-1)
				b = -data[v[0]].v[weak_classifier_id];
			else if(min_id == data.size() - 1)
				b = -data[v[data.size() - 1]].v[weak_classifier_id];
			else
				b = -(data[v[min_id]].v[weak_classifier_id] + data[v[min_id+1]].v[weak_classifier_id]) * 0.5;
		
			for(std::size_t i = 0; i < w.dim(); ++i)
				std::cout << w[i] << " ";
			std::cout << b << ": " << min_error << std::endl;
			weak_classifier_id++;
			return HyperPlane(w, b);
		}
		

		std::vector<double> weights;
		double weights_sum; 
		std::vector<ContactSpaceSampleData> data;

		std::vector<std::size_t> positive_indices;
		std::vector<std::size_t> negative_indices;

		flann::Index<ClassifierDistancer>* knn_index;

		std::size_t feature_dim;

		MultipleWeightedHyperPlane model;

		mutable RNG rng;

		std::size_t weak_classifier_id;
	};
	
}



#endif