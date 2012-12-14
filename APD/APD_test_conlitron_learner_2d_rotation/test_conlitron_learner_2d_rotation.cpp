#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>

namespace APDL
{
	void test_conlitron_learner_2d_rotation()
	{
		// scaled
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
			for(int i = 0; i < 10000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace);

			DataVector w(3);
			w[0] = 1; w[1] = 1; w[2] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.learn(contactspace.data, 3);

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

			learner.saveVisualizeData("conlitron_2d_rotation_vis.txt", contactspace.getScaler(), 100);
		}
	}
}

void main()
{
	APDL::test_conlitron_learner_2d_rotation();
}