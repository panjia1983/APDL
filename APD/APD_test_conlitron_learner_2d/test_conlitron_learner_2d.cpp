#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>

#include <APD/profile.h>


namespace APDL
{
	void test_conlitron_learner_2d()
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
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace_samples);

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			///////////////////////////////////////
			learner.use_approximate_dist = false;
			///////////////////////////////////////

			tools::Profiler::Begin("learn without approximate knn");
			learner.learn(contactspace_samples, 2);
			tools::Profiler::End("learn without approximate knn");

			std::cout << learner.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);
		}
	}

	void test_conlitron_learner_2d_approximate_knn()
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
			///////////////////////////////////////
			learner.use_approximate_dist = true;
			///////////////////////////////////////

			tools::Profiler::Begin("learn with approximate knn");
			learner.learn(contactspace_samples, 2);
			tools::Profiler::End("learn with approximate knn");

			learner.save("model.txt");

			std::cout << learner.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;

			MulticonlitronLearner learner_load(w, 0.01);
			learner_load.load("model.txt", "", false, 2);
			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner_load) << " " << errorRatioOnGrid(contactspace, learner_load, 100) << std::endl;



			MulticonlitronLearner learner2(w, 0.01, 0.01);
			///////////////////////////////////////
			learner.use_approximate_dist = true;
			///////////////////////////////////////

			tools::Profiler::Begin("learn with approximate knn2");
			learner2.learn(contactspace_samples, 2);
			tools::Profiler::End("learn with approximate knn2");

			std::cout << learner2.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner2) << " " << errorRatioOnGrid(contactspace, learner2, 100) << std::endl;



			MulticonlitronLearner learner3(w, 0.05, 0.01);
			///////////////////////////////////////
			learner3.use_approximate_dist = true;
			///////////////////////////////////////

			tools::Profiler::Begin("learn with approximate knn3");
			learner3.learn(contactspace_samples, 2);
			tools::Profiler::End("learn with approximate knn3");

			std::cout << learner3.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner3) << " " << errorRatioOnGrid(contactspace, learner3, 100) << std::endl;

		}
	}

}

void main()
{
	APDL::tools::Profiler::Start();

	// APDL::test_conlitron_learner_2d();
	APDL::test_conlitron_learner_2d_approximate_knn();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}