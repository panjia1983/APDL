#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>

namespace APDL
{
	void test_boosting_learner_2d()
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
			
			ContactSpaceR2 contactspace(p1, p2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);
			
			AdaBoostLearner learner;
			learner.learn(contactspace.data, contactspace.active_data_dim());

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			//for(std::size_t i = 0; i < results.size(); ++i)
			//{
			//	std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			//}
			//std::cout << std::endl;
		}
	}

}

void main()
{
	APDL::test_boosting_learner_2d();
}