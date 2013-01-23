#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>


namespace APDL
{
	void test_svm_learner_3d_cupspoon()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_cupspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}


	void test_svm_learner_3d_cupspoon_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_cupspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat_inc.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(2, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}










	void test_svm_learner_3d_starspoon()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/spoon.off");
		readObjFile(Q, "../data/star_large.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_starspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}

	void test_svm_learner_3d_starspoon_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/spoon.off");
		readObjFile(Q, "../data/star_large.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_starspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(2, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}



	void test_svm_learner_3d_torus()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Ringz/ringz.obj");
		readObjFile(Q, "../data/models/Ringz/ringz2.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_ringz.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		//// original learner, scaled
		//{
		//	std::vector<ContactSpaceSampleData> contactspace_samples;
		//	int n_iter = 10;
		//	
		//	for(int i = 0; i < n_iter; ++i)
		//	{
		//		tools::Profiler::Begin("original learner");
		//		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
		//		for(std::size_t j = 0; j < samples.size(); ++j)
		//			contactspace_samples.push_back(samples[j]);

		//		SVMLearner learner;
		//		learner.setDim(contactspace.active_data_dim());
		//		learner.setC(20);
		//		learner.setProbability(true);
		//		learner.setScaler(contactspace.getScaler());
		//		learner.setUseScaler(true);
		//		learner.setGamma(50); 

		//		learner.learn(contactspace_samples, contactspace.active_data_dim());

		//		tools::Profiler::End("original learner");

		//		active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
		//	}
		//	active_learning_stat_file << std::endl;	
		//}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(10000, 5000, 5000, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}


	void test_svm_learner_3d_torus_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Bullet/ringz.obj");
		readObjFile(Q, "../data/models/Bullet/ringz.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_ringz.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat_inc.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(2, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}



	void test_svm_learner_3d_teeth()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_teeth.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(20);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;	
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}


	void test_svm_learner_3d_teeth_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_teeth.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		//// original learner, scaled
		//{
		//	std::vector<ContactSpaceSampleData> contactspace_samples;
		//	int n_iter = 10;
		//	
		//	for(int i = 0; i < n_iter; ++i)
		//	{
		//		tools::Profiler::Begin("original learner");
		//		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(1000);
		//		for(std::size_t j = 0; j < samples.size(); ++j)
		//			contactspace_samples.push_back(samples[j]);

		//		SVMLearner learner;
		//		learner.setDim(contactspace.active_data_dim());
		//		learner.setC(20);
		//		learner.setProbability(true);
		//		learner.setScaler(contactspace.getScaler());
		//		learner.setUseScaler(true);
		//		learner.setGamma(50); 

		//		learner.learn(contactspace_samples, contactspace.active_data_dim());

		//		tools::Profiler::End("original learner");

		//		active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
		//	}
		//	active_learning_stat_file << std::endl;	
		//}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}
}

void main()
{
	APDL::tools::Profiler::Start();

	// APDL::test_svm_learner_3d_cupspoon();
	// APDL::test_svm_learner_3d_cupspoon_inc();

	// APDL::test_svm_learner_3d_starspoon();
	// APDL::test_svm_learner_3d_starspoon_inc();

	APDL::test_svm_learner_3d_torus();
	// APDL::test_svm_learner_3d_torus_inc();

	// APDL::test_svm_learner_3d_teeth_inc();


	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}