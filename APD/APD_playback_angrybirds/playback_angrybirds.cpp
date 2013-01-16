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

		std::string base_name = "../data/models/Box2D/angrybird_configs/dump_transform";

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/GreenPig32.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/BigRedBird17.polys");
		std::vector<Polygon> polys3 = readPolyFile("../data/models/Box2D/WhiteBird30.polys");

		double weights[3][7];

		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 7; ++j)
				weights[i][j] = 1;
		}

		{
			double Ix, Iy;
			inertia(polys1, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[0][0] = 1; weights[0][1] = 1; weights[0][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys2, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[1][0] = 1; weights[1][1] = 1; weights[1][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys3, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[2][0] = 1; weights[2][1] = 1; weights[2][2] = rotation_weight;
		}


		boost::timer t_construct;
		int n_angle = 30;
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_11 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys1, polys1, n_angle);
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_12 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys1, polys2, n_angle);
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_13 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys1, polys3, n_angle);
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_22 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys2, polys2, n_angle);
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_23 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys2, polys3, n_angle);
		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >  cspace_SE2_33 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(polys3, polys3, n_angle);

		timing_construct_file << t_construct.elapsed();
		timing_construct_file.flush();

		Collider2D collider11(polys1, polys1);
		Collider2D collider12(polys1, polys2);
		Collider2D collider13(polys1, polys3);
		Collider2D collider22(polys2, polys2);
		Collider2D collider23(polys2, polys3);
		Collider2D collider33(polys3, polys3);



		std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> >* cspaces[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				cspaces[i][j] = NULL;

		cspaces[0][0] = &cspace_SE2_11;
		cspaces[1][1] = &cspace_SE2_22;
		cspaces[2][2] = &cspace_SE2_33;
		cspaces[0][1] = &cspace_SE2_12;
		cspaces[0][2] = &cspace_SE2_13;
		cspaces[1][2] = &cspace_SE2_23;

		Collider2D* colliders[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				colliders[i][j] = NULL;

		colliders[0][0] = &collider11;
		colliders[0][1] = &collider12;
		colliders[0][2] = &collider13;
		colliders[1][1] = &collider22;
		colliders[1][2] = &collider23;
		colliders[2][2] = &collider33;


		for(std::size_t i = 0; i < 999; ++i)
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
				std::vector<Polygon>* obj1 = NULL;
				std::vector<Polygon>* obj2 = NULL;

				int id1, id2;

				if(frame[j][0].first == "GreenPig") { obj1 = &polys1; id1 = 0; }
				else if(frame[j][0].first == "RedBird") { obj1 = &polys2; id1 = 1; }
				else if(frame[j][0].first == "WhiteBird") { obj1 = &polys3; id1 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(frame[j][1].first == "GreenPig") { obj2 = &polys1; id2 = 0; }
				else if(frame[j][1].first == "RedBird") { obj2 = &polys2; id2 = 1; }
				else if(frame[j][1].first == "WhiteBird") { obj2 = &polys3; id2 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(id1 > id2)
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id1][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][1].second, frame[j][0].second);
					if(colliders[id2][id1]->isCollide(q))
					{
						std::pair<DataVector, double> PD_result = Minkowski_Cspace_2D::Exact_PD_SE2(q, *cspaces[id2][id1]);
						PD_file << PD_result.second << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
					
				}
				else
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id2][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
					if(colliders[id1][id2]->isCollide(q))
					{
						std::pair<DataVector, double> PD_result = Minkowski_Cspace_2D::Exact_PD_SE2(q, *cspaces[id1][id2]);
						PD_file << PD_result.second << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
				}
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}

	void playback_local_PD()
	{
		std::ofstream timing_file("timing_local_PD.txt");

		std::string base_name = "../data/models/Box2D/angrybird_configs/dump_transform";

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/GreenPig32.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/BigRedBird17.polys");
		std::vector<Polygon> polys3 = readPolyFile("../data/models/Box2D/WhiteBird30.polys");


		for(std::size_t i = 0; i < 999; ++i)
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
				std::vector<Polygon>* obj1 = NULL;
				std::vector<Polygon>* obj2 = NULL;

				int id1, id2;

				if(frame[j][0].first == "GreenPig") { obj1 = &polys1; id1 = 0; }
				else if(frame[j][0].first == "RedBird") { obj1 = &polys2; id1 = 1; }
				else if(frame[j][0].first == "WhiteBird") { obj1 = &polys3; id1 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(frame[j][1].first == "GreenPig") { obj2 = &polys1; id2 = 0; }
				else if(frame[j][1].first == "RedBird") { obj2 = &polys2; id2 = 1; }
				else if(frame[j][1].first == "WhiteBird") { obj2 = &polys3; id2 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				DataVector q = relative2D(frame[j][0].second, frame[j][1].second);

				boost::timer t;
				double pd = Collider2D::PDt(*obj1, *obj2, q);
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
		std::ofstream timing_file("timing_exact_APD.txt");
		std::ofstream timing_construct_file("timing_construct_APD.txt");

		std::string base_name = "../data/models/Box2D/angrybird_configs/dump_transform";

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/GreenPig32.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/BigRedBird17.polys");
		std::vector<Polygon> polys3 = readPolyFile("../data/models/Box2D/WhiteBird30.polys");

		double weights[3][7];

		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 7; ++j)
				weights[i][j] = 1;
		}

		{
			double Ix, Iy;
			inertia(polys1, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[0][0] = 1; weights[0][1] = 1; weights[0][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys2, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[1][0] = 1; weights[1][1] = 1; weights[1][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys3, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[2][0] = 1; weights[2][1] = 1; weights[2][2] = rotation_weight;
		}


		ContactSpaceSE2 contactspace11(polys1, polys1, 0.2 * (getCircle(polys1).second + getCircle(polys1).second));
		ContactSpaceSE2 contactspace12(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace13(polys1, polys3, 0.2 * (getCircle(polys1).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace22(polys2, polys2, 0.2 * (getCircle(polys2).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace23(polys2, polys3, 0.2 * (getCircle(polys2).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace33(polys3, polys3, 0.2 * (getCircle(polys3).second + getCircle(polys3).second));

		ContactSpaceSE2* contactspaces[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				contactspaces[i][j] = NULL;
		contactspaces[0][0] = &contactspace11;
		contactspaces[0][1] = &contactspace12;
		contactspaces[0][2] = &contactspace13;
		contactspaces[1][1] = &contactspace22;
		contactspaces[1][2] = &contactspace23;
		contactspaces[2][2] = &contactspace33;


		SVMLearner learner11, learner12, learner13, learner22, learner23, learner33;
		SVMLearner* learners[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				learners[i][j] = NULL;

		learners[0][0] = &learner11;
		learners[0][1] = &learner12;
		learners[0][2] = &learner13;
		learners[1][1] = &learner22;
		learners[1][2] = &learner23;
		learners[2][2] = &learner33;

		std::ofstream active_learning_stat_file("stat.txt");

		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				if(learners[i][j])
				{
					learners[i][j]->setC(50);
					learners[i][j]->setScaler(contactspaces[i][j]->getScaler());
					learners[i][j]->setUseScaler(true);
					learners[i][j]->setGamma(50);
					std::vector<ContactSpaceSampleData> contactspace_samples = contactspaces[i][j]->uniform_sample(100000);
					std::ostringstream  convert;
					convert << i << j;
					std::string scaler_file_name = "scaler_" + convert.str() + ".txt";
					std::string model_file_name = "model_" + convert.str() + ".txt";
					std::ofstream scaler_file(scaler_file_name.c_str());
					scaler_file << contactspaces[i][j]->getScaler() << std::endl;
					learners[i][j]->learn(contactspace_samples, contactspaces[i][j]->active_data_dim());
					learners[i][j]->save(model_file_name.c_str());
					std::cout << empiricalErrorRatio(contactspace_samples, *learners[i][j]) << " " << errorRatioOnGrid(*contactspaces[i][j], *learners[i][j], 20) << std::endl;

					//std::ostringstream  convert;
					//convert << i << j;
					//std::string scaler_file_name = "scaler_" + convert.str() + ".txt";
					//std::string model_file_name = "model_" + convert.str() + "_active";
					//learners[i][j]->setDim(contactspaces[i][j]->active_data_dim());
					//learners[i][j]->setC(50);
					//learners[i][j]->setScaler(contactspaces[i][j]->getScaler());

					//learners[i][j]->setUseScaler(true);
					//learners[i][j]->setGamma(50);

					//std::ofstream scaler_file(scaler_file_name.c_str());
					//scaler_file << contactspaces[i][j]->getScaler() << std::endl;

					//SpatialTreeEParam param;

					//param.max_depth = 6;
					//param.initial_depth = 3;
					//param.stop_abs_diff = 1e-4;
					//param.stop_related_diff = 0.02;
					//param.epsilon = 0;
					//param.result_eps = 0;

					//SVMEvaluator evaluator(*learners[i][j]);
					//FilterParam fparam;
					//// fparam.enable_filter = false;
					//DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, *learners[i][j]);
					//// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

					//ActiveLearningParam aparam(1000, 500, 500, 9);
					//aparam.debug = true;
					//aparam.num_grid = 20;
					//aparam.debug_os = &active_learning_stat_file;
					//aparam.model_name = model_file_name;
					//active_learning(*contactspaces[i][j], *learners[i][j], decision_boundary_sampler, aparam);
				}
			}
		}


		Collider2D* colliders[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				colliders[i][j] = NULL;

		Collider2D collider11(polys1, polys1);
		Collider2D collider12(polys1, polys2);
		Collider2D collider13(polys1, polys3);
		Collider2D collider22(polys2, polys2);
		Collider2D collider23(polys2, polys3);
		Collider2D collider33(polys3, polys3);

		colliders[0][0] = &collider11;
		colliders[0][1] = &collider12;
		colliders[0][2] = &collider13;
		colliders[1][1] = &collider22;
		colliders[1][2] = &collider23;
		colliders[2][2] = &collider33;

		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index11 = learner11.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index12 = learner12.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index13 = learner13.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index22 = learner22.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index23 = learner23.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index33 = learner33.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();


		flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_indices[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				query_indices[i][j] = NULL;
		query_indices[0][0] = query_index11;
		query_indices[0][1] = query_index12;
		query_indices[0][2] = query_index13;
		query_indices[1][1] = query_index22;
		query_indices[1][2] = query_index23;
		query_indices[2][2] = query_index33;


		for(std::size_t i = 0; i < 999; ++i)
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
				std::vector<Polygon>* obj1 = NULL;
				std::vector<Polygon>* obj2 = NULL;

				int id1, id2;

				if(frame[j][0].first == "GreenPig") { obj1 = &polys1; id1 = 0; }
				else if(frame[j][0].first == "RedBird") { obj1 = &polys2; id1 = 1; }
				else if(frame[j][0].first == "WhiteBird") { obj1 = &polys3; id1 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(frame[j][1].first == "GreenPig") { obj2 = &polys1; id2 = 0; }
				else if(frame[j][1].first == "RedBird") { obj2 = &polys2; id2 = 1; }
				else if(frame[j][1].first == "WhiteBird") { obj2 = &polys3; id2 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(id1 > id2)
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id1][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][1].second, frame[j][0].second);
					if(colliders[id2][id1]->isCollide(q))
					{
						QueryResult PD_result = PD_query(*learners[id2][id1], *contactspaces[id2][id1], query_indices[id2][id1], q);
						PD_file << PD_result.PD << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
					
				}
				else
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id2][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
					if(colliders[id1][id2]->isCollide(q))
					{
						QueryResult PD_result = PD_query(*learners[id1][id2], *contactspaces[id1][id2], query_indices[id1][id2], q);
						PD_file << PD_result.PD << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
				}
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}
	}


	void playback()
	{
		std::ofstream timing_file("timing_exact_APD.txt");
		std::ofstream timing_construct_file("timing_construct_APD.txt");

		std::string base_name = "../data/models/Box2D/angrybird_configs/dump_transform";

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/GreenPig32.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/BigRedBird17.polys");
		std::vector<Polygon> polys3 = readPolyFile("../data/models/Box2D/WhiteBird30.polys");

		double weights[3][7];

		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 7; ++j)
				weights[i][j] = 1;
		}

		{
			double Ix, Iy;
			inertia(polys1, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[0][0] = 1; weights[0][1] = 1; weights[0][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys2, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[1][0] = 1; weights[1][1] = 1; weights[1][2] = rotation_weight;
		}

		{
			double Ix, Iy;
			inertia(polys3, Ix, Iy);
			double rotation_weight = sqrt(Ix * Ix + Iy * Iy);
			weights[2][0] = 1; weights[2][1] = 1; weights[2][2] = rotation_weight;
		}


		ContactSpaceSE2 contactspace11(polys1, polys1, 0.2 * (getCircle(polys1).second + getCircle(polys1).second));
		ContactSpaceSE2 contactspace12(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace13(polys1, polys3, 0.2 * (getCircle(polys1).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace22(polys2, polys2, 0.2 * (getCircle(polys2).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace23(polys2, polys3, 0.2 * (getCircle(polys2).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace33(polys3, polys3, 0.2 * (getCircle(polys3).second + getCircle(polys3).second));

		ContactSpaceSE2* contactspaces[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				contactspaces[i][j] = NULL;
		contactspaces[0][0] = &contactspace11;
		contactspaces[0][1] = &contactspace12;
		contactspaces[0][2] = &contactspace13;
		contactspaces[1][1] = &contactspace22;
		contactspaces[1][2] = &contactspace23;
		contactspaces[2][2] = &contactspace33;


		SVMLearner learner11, learner12, learner13, learner22, learner23, learner33;
		SVMLearner* learners[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				learners[i][j] = NULL;

		learners[0][0] = &learner11;
		learners[0][1] = &learner12;
		learners[0][2] = &learner13;
		learners[1][1] = &learner22;
		learners[1][2] = &learner23;
		learners[2][2] = &learner33;

		std::ofstream active_learning_stat_file("stat.txt");

		int n_samples = 100000;
		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				if(learners[i][j])
				{
					learners[i][j]->setC(50);
					learners[i][j]->setScaler(contactspaces[i][j]->getScaler());
					learners[i][j]->setUseScaler(true);
					learners[i][j]->setGamma(50);
					std::vector<ContactSpaceSampleData> contactspace_samples = contactspaces[i][j]->uniform_sample(n_samples);
					std::ostringstream  convert;
					convert << i << j;
					std::string scaler_file_name = "scaler_" + convert.str() + ".txt";
					std::string model_file_name = "model_" + convert.str() + ".txt";
					std::ofstream scaler_file(scaler_file_name.c_str());
					scaler_file << contactspaces[i][j]->getScaler() << std::endl;
					learners[i][j]->learn(contactspace_samples, contactspaces[i][j]->active_data_dim());
					learners[i][j]->save(model_file_name.c_str());
					std::cout << empiricalErrorRatio(contactspace_samples, *learners[i][j]) << " " << errorRatioOnGrid(*contactspaces[i][j], *learners[i][j], 20) << std::endl;

					//std::ostringstream  convert;
					//convert << i << j;
					//std::string scaler_file_name = "scaler_" + convert.str() + ".txt";
					//std::string model_file_name = "model_" + convert.str() + "_active";
					//learners[i][j]->setDim(contactspaces[i][j]->active_data_dim());
					//learners[i][j]->setC(50);
					//learners[i][j]->setScaler(contactspaces[i][j]->getScaler());

					//learners[i][j]->setUseScaler(true);
					//learners[i][j]->setGamma(50);

					//std::ofstream scaler_file(scaler_file_name.c_str());
					//scaler_file << contactspaces[i][j]->getScaler() << std::endl;

					//SpatialTreeEParam param;

					//param.max_depth = 6;
					//param.initial_depth = 3;
					//param.stop_abs_diff = 1e-4;
					//param.stop_related_diff = 0.02;
					//param.epsilon = 0;
					//param.result_eps = 0;

					//SVMEvaluator evaluator(*learners[i][j]);
					//FilterParam fparam;
					//// fparam.enable_filter = false;
					//DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, *learners[i][j]);
					//// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);

					//ActiveLearningParam aparam(1000, 500, 500, 9);
					//aparam.debug = true;
					//aparam.num_grid = 20;
					//aparam.debug_os = &active_learning_stat_file;
					//aparam.model_name = model_file_name;
					//active_learning(*contactspaces[i][j], *learners[i][j], decision_boundary_sampler, aparam);
				}
			}
		}


		Collider2D* colliders[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				colliders[i][j] = NULL;

		Collider2D collider11(polys1, polys1);
		Collider2D collider12(polys1, polys2);
		Collider2D collider13(polys1, polys3);
		Collider2D collider22(polys2, polys2);
		Collider2D collider23(polys2, polys3);
		Collider2D collider33(polys3, polys3);

		colliders[0][0] = &collider11;
		colliders[0][1] = &collider12;
		colliders[0][2] = &collider13;
		colliders[1][1] = &collider22;
		colliders[1][2] = &collider23;
		colliders[2][2] = &collider33;

		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index11 = learner11.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index12 = learner12.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index13 = learner13.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index22 = learner22.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index23 = learner23.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();
		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_index33 = learner33.constructIndexOfSupportVectorsForQuery<ContactSpaceSE2, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();


		//flann::HierarchicalClusteringIndex<ContactSpaceSE2::DistanceType>* query_indices[3][3];
		//for(int i = 0; i < 3; ++i)
		//	for(int j = 0; j < 3; ++j)
		//		query_indices[i][j] = NULL;
		//query_indices[0][0] = query_index11;
		//query_indices[0][1] = query_index12;
		//query_indices[0][2] = query_index13;
		//query_indices[1][1] = query_index22;
		//query_indices[1][2] = query_index23;
		//query_indices[2][2] = query_index33;




		std::vector<ContactSpaceSampleData> support_samples11;
		std::vector<ContactSpaceSampleData> support_samples12;
		std::vector<ContactSpaceSampleData> support_samples13;
		std::vector<ContactSpaceSampleData> support_samples22;
		std::vector<ContactSpaceSampleData> support_samples23;
		std::vector<ContactSpaceSampleData> support_samples33;

		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex>* extended_models[3][3];
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 3; ++j)
				extended_models[i][j] = NULL;

		int knn_k = 50;

		learners[0][0]->collectSupportVectors(support_samples11);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model11 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[0][0], *learners[0][0], support_samples11, 0.01, knn_k);
		learners[0][1]->collectSupportVectors(support_samples12);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model12 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[0][1], *learners[0][1], support_samples12, 0.01, knn_k);
		learners[0][2]->collectSupportVectors(support_samples13);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model13 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[0][2], *learners[0][2], support_samples13, 0.01, knn_k);
		learners[1][1]->collectSupportVectors(support_samples22);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model22 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[1][1], *learners[1][1], support_samples22, 0.01, knn_k);
		learners[1][2]->collectSupportVectors(support_samples23);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model23 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[1][2], *learners[1][2], support_samples23, 0.01, knn_k);
		learners[2][2]->collectSupportVectors(support_samples33);
		ExtendedModel<ContactSpaceSE2, flann::HierarchicalClusteringIndex> extended_model33 = constructExtendedModelForModelDecisionBoundary<ContactSpaceSE2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(*contactspaces[2][2], *learners[2][2], support_samples33, 0.01, knn_k);

		extended_models[0][0] = &extended_model11;
		extended_models[0][1] = &extended_model12;
		extended_models[0][2] = &extended_model13;
		extended_models[1][1] = &extended_model22;
		extended_models[1][2] = &extended_model23;
		extended_models[2][2] = &extended_model33;


		for(std::size_t i = 0; i < 999; ++i)
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
				std::vector<Polygon>* obj1 = NULL;
				std::vector<Polygon>* obj2 = NULL;

				int id1, id2;

				if(frame[j][0].first == "GreenPig") { obj1 = &polys1; id1 = 0; }
				else if(frame[j][0].first == "RedBird") { obj1 = &polys2; id1 = 1; }
				else if(frame[j][0].first == "WhiteBird") { obj1 = &polys3; id1 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(frame[j][1].first == "GreenPig") { obj2 = &polys1; id2 = 0; }
				else if(frame[j][1].first == "RedBird") { obj2 = &polys2; id2 = 1; }
				else if(frame[j][1].first == "WhiteBird") { obj2 = &polys3; id2 = 2; }
				else 
				{
					std::cout << "error" << std::endl;
					return;
				}

				if(id1 > id2)
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id1][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][1].second, frame[j][0].second);
					if(colliders[id2][id1]->isCollide(q))
					{
						// QueryResult PD_result = PD_query(*learners[id2][id1], *contactspaces[id2][id1], query_indices[id2][id1], q);
						QueryResult PD_result = PD_query(*learners[id2][id1], *contactspaces[id2][id1], extended_models[id2][id1]->index, extended_models[id2][id1]->samples, q);
						PD_file << PD_result.PD << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
					
				}
				else
				{
					for(int k = 0; k < 7; ++k)
						distance_weight[k] = weights[id2][k];

					boost::timer t;
					DataVector q = relative2D(frame[j][0].second, frame[j][1].second);
					if(colliders[id1][id2]->isCollide(q))
					{
						// QueryResult PD_result = PD_query(*learners[id1][id2], *contactspaces[id1][id2], query_indices[id1][id2], q);
						QueryResult PD_result = PD_query(*learners[id1][id2], *contactspaces[id1][id2], extended_models[id1][id2]->index, extended_models[id1][id2]->samples, q);
						PD_file << PD_result.PD << " ";
					}
					else
					{
						PD_file << 0 << " ";
					}

					timing_per_frame += t.elapsed();
				}
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
	// APDL::playback_exact_CSpace();
	// APDL::playback_local_PD();
}