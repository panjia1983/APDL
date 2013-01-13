#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>


std::ofstream test_file; 

namespace APDL
{
	void test_svm_active_learning()
	{
		std::ifstream room_file("../data/rooms_star.dat");

		if(!room_file.is_open())
		{
			std::cerr << "Failed to open the input file." << std::endl;
			return;
		}

		Minkowski_Cspace_2D::Polygon_2 P, Q;

		room_file >> P >> Q;
		room_file.close();

		Polygon p1 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(P);
		Polygon p2 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(Q);

		ContactSpaceSE2 contactspace(p1, p2, 2);
		std::ofstream scaler_file("scaler_2d_rotation.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat.txt");

		//// original learner, not scaled
		//{
		//	std::vector<ContactSpaceSampleData> contactspace_samples;
		//	int n_iter = 30;
		//	
		//	for(int i = 0; i < n_iter; ++i)
		//	{
		//		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
		//		for(std::size_t j = 0; j < samples.size(); ++j)
		//			contactspace_samples.push_back(samples[j]);

		//		SVMLearner learner;
		//		learner.setDim(contactspace.active_data_dim());
		//		learner.setC(10);
		//		learner.setProbability(true);
		//		learner.setScaler(contactspace.getScaler());

		//		learner.learn(contactspace_samples, contactspace.active_data_dim());

		//		active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
		//	}
		//	active_learning_stat_file << std::endl;			
		//}

		//// active learner, not scaled
		//{
		//	SVMLearner learner;
		//	learner.setDim(contactspace.active_data_dim());
		//	learner.setC(10);
		//	learner.setProbability(true);
		//	learner.setScaler(contactspace.getScaler());

		//	SpatialTreeEParam param;
		//	param.max_depth = 8;
		//	param.initial_depth = 4;
		//	param.stop_abs_diff = 0.2;
		//	param.stop_related_diff = 0.1;
		//	param.epsilon = 0;
		//	param.result_eps = 0;


		//	SVMEvaluator evaluator(learner);
		//	FilterParam fparam;
		//	fparam.enable_filter = false;
		//	DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

		//	ActiveLearningParam aparam(100, 50, 50, 29);
		//	aparam.debug = true;
		//	aparam.debug_os = &active_learning_stat_file;
		//	aparam.model_name = "model_2d_active";
		//	aparam.num_grid = 20;
		//	active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		//}

		//// active learner2, not scaled
		//{
		//	SVMLearner learner;
		//	learner.setDim(contactspace.active_data_dim());
		//	learner.setC(10);
		//	learner.setProbability(true);
		//	learner.setScaler(contactspace.getScaler());

		//	SpatialTreeEParam param;
		//	param.max_depth = 8;
		//	param.initial_depth = 4;
		//	param.stop_abs_diff = 0.2;
		//	param.stop_related_diff = 0.1;
		//	param.epsilon = 0;
		//	param.result_eps = 0;

		//	SVMEvaluator evaluator(learner);
		//	FilterParam fparam;
		//	fparam.enable_filter = false;
		//	DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

		//	ActiveLearningParam aparam(100, 50, 50, 58);
		//	aparam.debug = true;
		//	aparam.debug_os = &active_learning_stat_file;
		//	aparam.model_name = "model_2d_activeb";
		//	aparam.num_grid = 20;
		//	active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		//}


		// original learner, scaled
		{

			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 19);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 38);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}





	void test_svm_active_learning_inc()
	{
		std::ifstream room_file("../data/rooms_star.dat");

		if(!room_file.is_open())
		{
			std::cerr << "Failed to open the input file." << std::endl;
			return;
		}

		Minkowski_Cspace_2D::Polygon_2 P, Q;

		room_file >> P >> Q;
		room_file.close();

		Polygon p1 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(P);
		Polygon p2 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(Q);

		ContactSpaceSE2 contactspace(p1, p2, 2);
		std::ofstream scaler_file("scaler_2d.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("stat_inc.txt");

		//// original learner, not scaled
		//{
		//	std::vector<ContactSpaceSampleData> contactspace_samples;
		//	int n_iter = 10;

		//	for(int i = 0; i < n_iter; ++i)
		//	{
		//		std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
		//		for(std::size_t j = 0; j < samples.size(); ++j)
		//			contactspace_samples.push_back(samples[j]);

		//		SVMLearner learner;
		//		learner.setDim(contactspace.active_data_dim());
		//		learner.setC(10);
		//		learner.setProbability(true);
		//		learner.setScaler(contactspace.getScaler());

		//		learner.learn(contactspace_samples, contactspace.active_data_dim());

		//		active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << " ";
		//	}
		//	active_learning_stat_file << std::endl;
		//}


		//// active learner, not scaled
		//{
		//	SVMLearner learner;
		//	learner.setDim(contactspace.active_data_dim());
		//	learner.setC(10);
		//	learner.setProbability(true);
		//	learner.setScaler(contactspace.getScaler());

		//	SpatialTreeEParam param;
		//	param.max_depth = 8;
		//	param.initial_depth = 4;
		//	param.stop_abs_diff = 0.2;
		//	param.stop_related_diff = 0.1;
		//	param.epsilon = 0;
		//	param.result_eps = 0;

		//	SVMEvaluator evaluator(learner);
		//	FilterParam fparam;
		//	DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

		//	ActiveLearningParam aparam(100, 50, 50, 9);
		//	aparam.debug = true;
		//	aparam.debug_os = &active_learning_stat_file;
		//	active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		//}

		// original learner, scaled
		{
			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << " ";
			}
			active_learning_stat_file << std::endl;
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 8;
			param.initial_depth = 4;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0.002;
			param.result_eps = 1e-3;


			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 18);
			aparam.debug = true;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_inc";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}





	void test_svm_active_learning_poly_spider()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");


		ContactSpaceSE2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::ofstream scaler_file("scaler_2d_spider.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("spider_stat.txt");

		// original learner, scaled
		{

			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 19);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 38);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}

	void test_svm_active_learning_poly_spider_inc()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");


		ContactSpaceSE2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::ofstream scaler_file("scaler_2d_spider.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("spider_stat_inc.txt");

		// original learner, scaled
		{

			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 19);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}


	void test_svm_active_learning_poly_monkey()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");


		ContactSpaceSE2 contactspace(polys1, polys2, 2);
		std::ofstream scaler_file("scaler_2d_monkey.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("monkey_stat.txt");

		// original learner, scaled
		{

			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 19);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_active";
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 38);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_activeb";
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}

	void test_svm_active_learning_poly_monkey_inc()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");


		ContactSpaceSE2 contactspace(polys1, polys2, 2);
		std::ofstream scaler_file("scaler_2d_monkey.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::ofstream active_learning_stat_file("monkey_stat.txt");

		// original learner, scaled
		{

			std::vector<ContactSpaceSampleData> contactspace_samples;
			int n_iter = 20;

			for(int i = 0; i < n_iter; ++i)
			{
				std::vector<ContactSpaceSampleData> samples = contactspace.uniform_sample(100);
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

				active_learning_stat_file << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << " ";
			}
			active_learning_stat_file << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(50);

			SpatialTreeEParam param;
			//param.max_depth = 8;
			//param.initial_depth = 4;
			//param.stop_abs_diff = 0.2;
			//param.stop_related_diff = 0.1;
			//param.epsilon = 0;
			//param.result_eps = 0;

			param.max_depth = 6;
			param.initial_depth = 3;
			param.stop_abs_diff = 1e-4;
			param.stop_related_diff = 0.02;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.enable_filter = false;
			// DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

			ActiveLearningParam aparam(100, 50, 50, 19);
			aparam.debug = true;
			aparam.num_grid = 20;
			aparam.debug_os = &active_learning_stat_file;
			aparam.model_name = "model_2d_scaled_active";
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}
}

void main()
{
	APDL::tools::Profiler::Start();

	test_file.open("active_test.txt");
	
	// APDL::test_svm_active_learning();
	// APDL::test_svm_active_learning_inc();

	// APDL::test_svm_active_learning_poly_spider();
	// APDL::test_svm_active_learning_poly_spider_inc();

	// APDL::test_svm_active_learning_poly_monkey();
	APDL::test_svm_active_learning_poly_monkey_inc();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}