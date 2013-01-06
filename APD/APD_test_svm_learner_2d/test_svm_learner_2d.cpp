#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>


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
			learner.setUseScaler(true);


			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;


			SVMLearner learner_load;
			learner_load.load("model_2d.txt", "scaler_2d.txt", true, 2);
			std::cout << contactspace_samples.size() << ": " << std::endl;
			std::cout << empiricalErrorRatio(contactspace_samples, learner_load) << std::endl;
			std::cout << errorRatioOnGrid(contactspace, learner_load, 100) << std::endl;


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



	void test_svm_learner_poly_spiders()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");


		ContactSpaceR2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_2d_spiders.txt");
		asciiWriter(out, contactspace_samples);

		SVMLearner learner;
		learner.setC(10);
		learner.setProbability(true);
		learner.setScaler(contactspace.getScaler());
		learner.setGamma(20);

		std::ofstream scaler_file("scaler_2d_spiders.txt");
		scaler_file << contactspace.getScaler() << std::endl;


		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model_2d_spiders.txt");

		std::vector<PredictResult> results = learner.predict(contactspace_samples);

		//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
		//std::cout << std::endl;


		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
	}

	void test_svm_learner_poly_spidertooth()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/tooth.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys");


		ContactSpaceR2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_2d_spidertooth.txt");
		asciiWriter(out, contactspace_samples);

		SVMLearner learner;
		learner.setC(10);
		learner.setProbability(true);
		learner.setScaler(contactspace.getScaler());
		learner.setGamma(20);

		std::ofstream scaler_file("scaler_2d_spidertooth.txt");
		scaler_file << contactspace.getScaler() << std::endl;


		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model_2d_spidertooth.txt");

		std::vector<PredictResult> results = learner.predict(contactspace_samples);

		//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
		//std::cout << std::endl;


		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
	}


	void test_svm_learner_poly_monkeys()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");


		ContactSpaceR2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_2d_monkeys.txt");
		asciiWriter(out, contactspace_samples);

		SVMLearner learner;
		learner.setC(10);
		learner.setProbability(true);
		learner.setScaler(contactspace.getScaler());
		learner.setGamma(20);

		std::ofstream scaler_file("scaler_2d_monkeys.txt");
		scaler_file << contactspace.getScaler() << std::endl;


		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model_2d_monkeys.txt");

		std::vector<PredictResult> results = learner.predict(contactspace_samples);

		//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
		//std::cout << std::endl;


		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
	}

	void test_svm_learner_poly_monkeytooth()
	{
		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/tooth.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_monkey550.polys");


		ContactSpaceR2 contactspace(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_2d_monkeytooth.txt");
		asciiWriter(out, contactspace_samples);

		SVMLearner learner;
		learner.setC(10);
		learner.setProbability(true);
		learner.setScaler(contactspace.getScaler());
		learner.setGamma(20);

		std::ofstream scaler_file("scaler_2d_monkeytooth.txt");
		scaler_file << contactspace.getScaler() << std::endl;


		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model_2d_monkeytooth.txt");

		std::vector<PredictResult> results = learner.predict(contactspace_samples);

		//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
		//std::cout << std::endl;


		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;
	}
}

void main()
{
	APDL::test_svm_learner_2d();
	// APDL::test_svm_learner_poly_spiders();
	// APDL::test_svm_learner_poly_spidertooth();
	// APDL::test_svm_learner_poly_monkeys();
	// APDL::test_svm_learner_poly_monkeytooth();
}