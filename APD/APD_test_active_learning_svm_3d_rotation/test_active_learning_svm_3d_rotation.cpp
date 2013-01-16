#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>


namespace APDL
{
	void test_svm_learner_3d_rotation_cupspoon()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Quat contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_cupspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(10000);
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

				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = "original_model_" + convert.str() + ".txt";
				learner.save(model_file_name);

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
				active_learning_stat_file.flush();
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
			param.max_depth = 3;
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
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
			param.max_depth = 3;
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


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d__rotation_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}


	void test_svm_learner_3d_rotation_cupspoon_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_cupspoon.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}










	void test_svm_learner_3d_rotation_starspoon()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/spoon.off");
		readObjFile(Q, "../data/star_large.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}

	void test_svm_learner_3d_rotation_starspoon_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/spoon.off");
		readObjFile(Q, "../data/star_large.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}



	void test_svm_learner_3d_rotation_torus()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Bullet/ringz.obj");
		readObjFile(Q, "../data/models/Bullet/ringz.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		AABB3D aabb = computeAABB(P);

		// ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		double scale_x = aabb.b_max[0] - aabb.b_min[0];
		double scale_y = aabb.b_max[1] - aabb.b_min[1];
		double scale_z = aabb.b_max[2] - aabb.b_min[2];
		ContactSpaceSE3Euler2 contactspace(P, Q, scale_x, scale_y, scale_z);
		std::ofstream scaler_file("scaler_3d_rotation_ringz.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 10;
			
			for(int i = 0; i < n_iter; ++i)
			{
				std::cout << i << std::endl;
				tools::Profiler::Begin("original learner");
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(10000);
				for(std::size_t j = 0; j < samples.size(); ++j)
					contactspace_samples.push_back(samples[j]);

				SVMLearner learner;
				learner.setDim(contactspace.active_data_dim());
				learner.setC(50);
				learner.setProbability(true);
				learner.setScaler(contactspace.getScaler());
				learner.setUseScaler(true);
				learner.setGamma(50); 

				learner.learn(contactspace_samples, contactspace.active_data_dim());

				tools::Profiler::End("original learner");

				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = "original_model_" + convert.str() + ".txt";
				learner.save(model_file_name);

				// active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";

				std::vector<ContactSpaceSampleData> test_samples = contactspace.uniform_sample(20000);
				active_learning_stat_file << contactspace_samples.size() << " ";
				double r1 = empiricalErrorRatio(contactspace_samples, learner);
				// double r2 = errorRatioOnGrid(contactspace, learner, 5);
				double r2 = empiricalErrorRatio(test_samples, learner);
				active_learning_stat_file << r1 << " " << r2 << std::endl;
				active_learning_stat_file.flush();



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
			learner.setGamma(20); 

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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(20); 

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


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}


	void test_svm_learner_3d_rotation_torus_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Bullet/ringz.obj");
		readObjFile(Q, "../data/models/Bullet/ringz.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_ringz.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}







	void test_svm_learner_3d_rotation_teethsmall()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new_small.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new_small.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = "original_model_" + convert.str() + ".txt";
				learner.save(model_file_name);

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}

	void test_svm_learner_3d_rotation_teethsmall_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new_small.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new_small.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}





	void test_svm_learner_3d_rotation_teeth()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				std::ostringstream  convert;
				convert << i;
				std::string model_file_name = "original_model_" + convert.str() + ".txt";
				learner.save(model_file_name);

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
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
			param.max_depth = 3;
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
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(50, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}

	void test_svm_learner_3d_rotation_teeth_inc()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_starspoon.txt");
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << " ";
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
			aparam.num_grid = 5;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_3d_rotation_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}
}

void main()
{
	APDL::tools::Profiler::Start();

	APDL::test_svm_learner_3d_rotation_cupspoon();
	// APDL::test_svm_learner_3d_rotation_cupspoon_inc();

	// APDL::test_svm_learner_3d_rotation_starspoon();
	// APDL::test_svm_learner_3d_rotation_starspoon_inc();

	// APDL::test_svm_learner_3d_rotation_torus();
	// APDL::test_svm_learner_3d_rotation_torus_inc();

	// APDL::test_svm_learner_3d_rotation_teethsmall();


	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}