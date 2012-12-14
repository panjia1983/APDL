#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_distance.h>

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
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.learn(contactspace.data, 2);

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			}


			int error_num = 0;
			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);

			MulticonlitronDistanceToDecisionBoundary_BruteForce distancer1(learner);
			MulticonlitronDistanceToDecisionBoundary_KNN distancer2(learner);
			MulticonlitronDistanceToDecisionBoundary_EmbedKNN distancer3(learner);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				const DataVector& v = contactspace.data[i].v;
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