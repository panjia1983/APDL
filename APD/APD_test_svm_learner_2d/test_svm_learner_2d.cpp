#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{
	void test_svm_learner_2d()
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
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			// learner.setUseScaler(true);


			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
				std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			std::cout << std::endl;
		}

		return;

		// linear model
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
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);
				
			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setLinearClassifier();
			learner.setC(10);
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			HyperPlane hp = learner.getLinearModel();
			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			{
				DataVector v(contactspace.active_data_dim());
				for(std::size_t j = 0; j < contactspace.active_data_dim(); ++j)
					v[j] = contactspace_samples[i].v[j];
				double pred = hp.evaluate(v);
				if(pred > 0) 
					std::cout << 1 << " ";
				else 
					std::cout << 0 << " ";
			}
			std::cout << std::endl;

			for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
				std::cout << results[i].label << " ";
			std::cout << std::endl;
		}

	}

}

void main()
{
	APDL::test_svm_learner_2d();
}