#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_distance.h>

namespace APDL
{	
	void test_distance_to_decision_boundary_svm()
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
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);

			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model.txt");

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			std::cout << std::endl;

			std::cout << learner.hyperw_normsqr << std::endl;


			SVMDistanceToDecisionBoundary_Bruteforce distancer1(learner);
			SVMDistanceToDecisionBoundary_OptimizationGradient distancer2(learner);
		    SVMDistanceToDecisionBoundary_Projection distancer3(learner);
			SVMDistanceToDecisionBoundary_Optimization distancer4(learner);
			SVMDistanceToDecisionBoundary_RoughLowerBound distancer5(learner);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				DataVector v = (learner.scaler && learner.use_scaler) ? contactspace.getScaler().scale(contactspace.data[i].v) : contactspace.data[i].v;
				double d1 = distancer1.distance(v);
				double d2 = distancer2.distance(v);
				double d3 = distancer3.distance(v);
				double d4 = distancer4.distance(v);
				double d5 = distancer5.distance(v);
				std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << " " << d5 << std::endl;
			}

			return;
		}

	}

}

void main()
{
	APDL::test_distance_to_decision_boundary_svm();
}