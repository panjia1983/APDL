#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>

#include <APD/profile.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{
	void test_conlitron_learner_2d_rotation()
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
			
			ContactSpaceSE2 contactspace(p1, p2, 2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace_samples);

			DataVector w(3);
			w[0] = 1; w[1] = 1; w[2] = 1;
			MulticonlitronLearner learner(w, 0.01);
			///////////////////////////////////////
			learner.use_approximate_dist = false;
			///////////////////////////////////////

			tools::Profiler::Begin("learn without approximate knn");
			learner.learn(contactspace_samples, 3);
			tools::Profiler::End("learn without approximate knn");

			std::cout << learner.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 50) << std::endl;
		}
	}

	void test_conlitron_learner_2d_rotation_approximate_knn()
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

			ContactSpaceSE2 contactspace(p1, p2, 2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace_samples);

			DataVector w(3);
			w[0] = 1; w[1] = 1; w[2] = 1;
			MulticonlitronLearner learner(w, 0.01);
			///////////////////////////////////////
			learner.use_approximate_dist = true;
			///////////////////////////////////////

			tools::Profiler::Begin("learn with approximate knn");
			learner.learn(contactspace_samples, 3);
			tools::Profiler::End("learn with approximate knn");

			std::cout << learner.model.numOfHyperPlanes() << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 50) << std::endl;


		}
	}
}

void main()
{
	APDL::tools::Profiler::Start();
	
	APDL::test_conlitron_learner_2d_rotation();
	APDL::test_conlitron_learner_2d_rotation_approximate_knn();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}