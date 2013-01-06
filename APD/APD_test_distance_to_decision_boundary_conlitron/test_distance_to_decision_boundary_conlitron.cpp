#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_distance.h>
#include <APD/active_learning.h>


namespace APDL
{	
	void test_distance_to_decision_boundary_conlitron()
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
			learner.learn(contactspace_samples, 2);

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);

			MulticonlitronDistanceToDecisionBoundary_BruteForce distancer1(learner);
			MulticonlitronDistanceToDecisionBoundary_KNN distancer2(learner);
			MulticonlitronDistanceToDecisionBoundary_EmbedKNN distancer3(learner);

			for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			{
				const DataVector& v = contactspace_samples[i].v;
				double d1 = distancer1.distance(v);
				double d2 = distancer2.distance(v);
				double d3 = distancer3.distance(v);
				if(d1 > 5 || d2 > 5 || d3 > 5) 
					std::cout << d1 << " " << d2 << " " << d3 << std::endl;
			}

		}
	}
}

void main()
{
	APDL::test_distance_to_decision_boundary_conlitron();
}