#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{
	void test_online_query_svm()
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

			Minkowski_Cspace_2D::Polygon_with_holes_2 cspace_R2 = Minkowski_Cspace_2D::Minkowski_Cobstacle_R2(P, Q);

			Polygon p1 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(P);
			Polygon p2 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(Q);

			ContactSpaceR2 contactspace(p1, p2, 2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace_samples);

			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;


			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;


			flann::Index<ContactSpaceR2::DistanceType>* query_index = learner.constructIndexForQuery<ContactSpaceR2, flann::KDTreeIndexParams>();

			std::vector<ContactSpaceSampleData> query_samples = contactspace.uniform_sample(100);

			for(std::size_t i = 0; i < query_samples.size(); ++i)
			{
				QueryResult apprx_PD = PD_query(learner, contactspace, query_index, query_samples[i].v);
				QueryResult apprx_PD2 = PD_query2(learner, contactspace, query_index, query_samples[i].v);
				std::pair<DataVector, double> exact_PD = Minkowski_Cspace_2D::Exact_PD_R2(query_samples[i].v, cspace_R2);

				std::cout << apprx_PD.v[0] << " " << apprx_PD.v[1] << " " << apprx_PD.PD << std::endl;
				std::cout << apprx_PD2.v[0] << " " << apprx_PD2.v[1] << " " << apprx_PD2.PD << std::endl;
				std::cout << exact_PD.first[0] << " " << exact_PD.first[1] << " " << exact_PD.second << std::endl;
				std::cout << std::endl;
			}

			delete query_index;
		}
	}

	void test_online_query_conlitron()
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

			Minkowski_Cspace_2D::Polygon_with_holes_2 cspace_R2 = Minkowski_Cspace_2D::Minkowski_Cobstacle_R2(P, Q);

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

			std::vector<PredictResult> results = learner.predict(contactspace_samples);

			std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;


			flann::Index<ContactSpaceR2::DistanceType>* query_index = learner.constructIndexForQuery<ContactSpaceR2, flann::KDTreeIndexParams>();

			std::vector<ContactSpaceSampleData> query_samples = contactspace.uniform_sample(100);

			for(std::size_t i = 0; i < query_samples.size(); ++i)
			{
				QueryResult apprx_PD = PD_query(learner, contactspace, query_index, query_samples[i].v);
				QueryResult apprx_PD2 = PD_query2(learner, contactspace, query_index, query_samples[i].v);
				std::pair<DataVector, double> exact_PD = Minkowski_Cspace_2D::Exact_PD_R2(query_samples[i].v, cspace_R2);

				std::cout << "query: " << query_samples[i].v[0] << " " << query_samples[i].v[1] << std::endl;
				std::cout << apprx_PD.v[0] << " " << apprx_PD.v[1] << " " << apprx_PD.PD << std::endl;
				std::cout << apprx_PD2.v[0] << " " << apprx_PD2.v[1] << " " << apprx_PD2.PD << std::endl;
				std::cout << exact_PD.first[0] << " " << exact_PD.first[1] << " " << exact_PD.second << std::endl;
				std::cout << std::endl;
			}

			delete query_index;
		}
	}
}

void main()
{
	//APDL::test_online_query_svm();
	APDL::test_online_query_conlitron();
}