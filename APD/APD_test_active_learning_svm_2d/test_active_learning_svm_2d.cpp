#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>

void* user_conlitron_model;
double* user_conlitron_data;

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

		ContactSpaceR2 contactspace(p1, p2, 2);
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

		std::ofstream out("space_test_2d.txt");
		asciiWriter(out, contactspace_samples);

		// original learner, not scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;


			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}

		// active learner, not scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			SpatialTreeEParam param;
			param.max_depth = 8;
			param.initial_depth = 4;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;


			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// original learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(20);

			 

			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;


			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}


		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(20);

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
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
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

		ContactSpaceR2 contactspace(p1, p2, 2);
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

		std::ofstream out("space_test_2d.txt");
		asciiWriter(out, contactspace_samples);

		// original learner, not scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;


			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}


		// active learner, not scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			SpatialTreeEParam param;
			param.max_depth = 8;
			param.initial_depth = 4;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// original learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(5);
			
			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;


			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setDim(contactspace.active_data_dim());
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);
			learner.setGamma(5);

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
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}





	void test_svm_active_learning_poly_spider()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/tooth.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");


		ContactSpaceR2 contactspace(polys1, polys2, 2);




	}
}

void main()
{
	APDL::test_svm_active_learning_poly_spider();
	return;
	APDL::test_svm_active_learning();
	APDL::test_svm_active_learning_inc();
}