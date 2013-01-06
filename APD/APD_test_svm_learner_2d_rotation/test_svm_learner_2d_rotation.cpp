#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>


namespace APDL
{
	void test_svm_learner_2d_rotation()
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
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			/// we need to change gamma!!! for scaled version
			learner.setGamma(0.5 * contactspace.getScaler().getScale() * 10);


			std::ofstream scaler_file("scaler_2d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d_rotation.txt");

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;


		}


		return;

		// no scaled
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
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			//learner.setUseScaler(true);

			std::ofstream scaler_file("scaler_2d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d_rotation.txt");

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
		}
	}
}

void main()
{
	APDL::test_svm_learner_2d_rotation();
}