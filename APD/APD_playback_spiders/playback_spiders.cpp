#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/math_utility.h>

#include <APD/profile.h>

#include <boost/timer.hpp>

namespace APDL
{
	extern double distance_weight[7];

	void playback_exact_CSpace()
	{
		std::ofstream timing_file("timing_exact_SE2.txt");
		std::ofstream timing_construct_file("timing_construct_exact_SE2.txt");

		std::string base_name = "../data/models/Box2D/spider_cinfigs/dump_transform";

		std::vector<Polygon> poly = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);

		for(std::size_t i = 0; i < 7; ++i)
			distance_weight[i] = 1;

		{
			double Ix, Iy;
			inertia(poly, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = rotation_weight;
		}

		boost::timer t_construct;
		int n_angle = 30;
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(poly, poly, n_angle);

		timing_construct_file << t_construct.elapsed();
		timing_construct_file.flush();

		Collider2D collider(poly, poly);


		for(std::size_t i = 1001; i < 9999; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t j = 0; j < 4 - len; ++j)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::string PD_file_name= std::string("PD_exact_SE2") + ret + ".txt";


			std::ofstream PD_file(PD_file_name.c_str());

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile(filename);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;

				if(frame[j][0].first == "Wall" || frame[j][1].first == "Wall") continue;

				boost::timer t;
				DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
				if(collider.isCollide(q))
				{
					std::pair<DataVector, double> PD_result = Minkowski_Cspace_2D::Exact_PD_SE2(q, cspace);
					PD_file << PD_result.second << " ";
				}
				else
				{
					PD_file << 0 << " ";
				}

				timing_per_frame += t.elapsed();
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}


	void playback_local_PD()
	{
		std::ofstream timing_file("timing_local_PD.txt");

		std::string base_name = "../data/models/Box2D/spider_cinfigs/dump_transform";

		std::vector<Polygon> poly = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);

		for(std::size_t i = 1001; i < 9999; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t j = 0; j < 4 - len; ++j)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::string PD_file_name= std::string("PD_local") + ret + ".txt";


			std::ofstream PD_file(PD_file_name.c_str());

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile(filename);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;
				if(frame[j][0].first == "Wall" || frame[j][1].first == "Wall") continue;

				DataVector q = relative2D(frame[j][0].second, frame[j][1].second);

				boost::timer t;
				double pd = Collider2D::PDt(poly, poly, q);
				timing_per_frame += t.elapsed();
				PD_file << pd << " ";
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}


	void playback1()
	{
		std::ofstream timing_file("timing_APD.txt");
		std::ofstream timing_construct_file("timing_construct_APD.txt");

		std::string base_name = "../data/models/Box2D/spider_cinfigs/dump_transform";

		std::vector<Polygon> poly = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);

		for(std::size_t i = 0; i < 7; ++i)
			distance_weight[i] = 1;

		{
			double Ix, Iy;
			inertia(poly, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = rotation_weight;
		}


		ContactSpaceSE2 contactspace(poly, poly, 0.2 * (getCircle(poly).second + getCircle(poly).second));

		SVMLearner learner;

		learner.setC(50);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50);
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(100000);
		std::ofstream scaler_file("scaler.txt");
		scaler_file << contactspace.getScaler() << std::endl;
		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");
		std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << std::endl;

		Collider2D collider(poly, poly);

		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();


		for(std::size_t i = 1001; i < 9999; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t j = 0; j < 4 - len; ++j)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::string PD_file_name= std::string("PD_APD") + ret + ".txt";


			std::ofstream PD_file(PD_file_name.c_str());

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile(filename);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;

				if(frame[j][0].first == "Wall" || frame[j][1].first == "Wall") continue;

				boost::timer t;
				DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
				if(collider.isCollide(q))
				{
					QueryResult PD_result = PD_query(learner, contactspace, query_index, q);
					PD_file << PD_result.PD << " ";
				}
				else
				{
					PD_file << 0 << " ";
				}

				timing_per_frame += t.elapsed();
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}

	void playback()
	{
		std::ofstream timing_file("timing_APD.txt");
		std::ofstream timing_construct_file("timing_construct_APD.txt");

		std::string base_name = "../data/models/Box2D/spider_cinfigs/dump_transform";

		std::vector<Polygon> poly = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);

		for(std::size_t i = 0; i < 7; ++i)
			distance_weight[i] = 1;

		{
			double Ix, Iy;
			inertia(poly, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = rotation_weight;
		}


		ContactSpaceSE2 contactspace(poly, poly, 0.2 * (getCircle(poly).second + getCircle(poly).second));

		SVMLearner learner;

		learner.setC(50);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50);
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(100000);
		std::ofstream scaler_file("scaler.txt");
		scaler_file << contactspace.getScaler() << std::endl;
		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");
		std::cout << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << std::endl;

		Collider2D collider(poly, poly);

		// flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();


		int knn_k = 50;

		std::vector<ContactSpaceSampleData> support_samples;
		learner.collectSupportVectors(support_samples);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(contactspace, learner, support_samples, 0.01, knn_k);

		for(std::size_t i = 1001; i < 9999; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t j = 0; j < 4 - len; ++j)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::string PD_file_name= std::string("PD_APD") + ret + ".txt";


			std::ofstream PD_file(PD_file_name.c_str());

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile(filename);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;

				if(frame[j][0].first == "Wall" || frame[j][1].first == "Wall") continue;

				boost::timer t;
				DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
				if(collider.isCollide(q))
				{
					// QueryResult PD_result = PD_query(learner, contactspace, query_index, q);
					QueryResult PD_result = PD_query(learner, contactspace, extended_model.index, extended_model.samples, q);
					PD_file << PD_result.PD << " ";
				}
				else
				{
					PD_file << 0 << " ";
				}

				timing_per_frame += t.elapsed();
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}
}

void main()
{
	APDL::playback();

	// APDL::playback_local_PD();
	// APDL::playback_exact_CSpace();
}