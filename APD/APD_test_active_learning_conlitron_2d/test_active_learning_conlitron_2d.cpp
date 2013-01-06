#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>



namespace APDL
{
	void test_conlitron_active_learning()
	{
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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setScaler(contactspace.getScaler());
			learner.use_approximate_dist = true;
			learner.learn(contactspace_samples, 2);

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}


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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setDim(contactspace.active_data_dim());
			learner.setScaler(contactspace.getScaler());
			learner.use_approximate_dist = true;
			learner.learn(contactspace_samples, 2);

			SpatialTreeEParam param;
			param.max_depth = 10;
			param.initial_depth = 3;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;	
			param.result_eps = 0;

			MulticonlitronLearner evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<MulticonlitronLearner, MulticonlitronEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setDim(contactspace.active_data_dim());
			learner.setScaler(contactspace.getScaler());
			learner.use_approximate_dist = true;
			learner.learn(contactspace_samples, 2);

			SpatialTreeEParam param;
			param.max_depth = 10;
			param.initial_depth = 3;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;	
			param.result_eps = 0;

			MulticonlitronLearner evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<MulticonlitronLearner, MulticonlitronEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 18);
			aparam.debug = true;
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}



	void test_conlitron_active_learning_inc()
	{
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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setScaler(contactspace.getScaler());
			learner.learn(contactspace_samples, 2);

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}


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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setDim(contactspace.active_data_dim());
			learner.setScaler(contactspace.getScaler());
			learner.learn(contactspace_samples, 2);

			SpatialTreeEParam param;
			param.max_depth = 10;
			param.initial_depth = 3;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;	
			param.result_eps = 0;

			MulticonlitronLearner evaluator(learner);
			FilterParam fparam;
			DecisionBoundaryHierarchialTreeESampler<MulticonlitronLearner, MulticonlitronEvaluator> decision_boundary_sampler(param, fparam, learner);

			ActiveLearningParam aparam(100, 50, 50, 9);
			aparam.debug = true;
			active_learning_incremental(contactspace, learner, decision_boundary_sampler, aparam);
		}
	}
}

void main()
{
	APDL::test_conlitron_active_learning();
	// APDL::test_conlitron_active_learning_inc();
}