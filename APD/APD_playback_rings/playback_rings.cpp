#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/math_utility.h>

#include <APD/profile.h>

#include <boost/timer.hpp>
#include <APD/stopwatch.h>

Stopwatch aTimer;

namespace APDL
{
	extern double distance_weight[7];

	void playback_exact_CSpace()
	{
		bool use_euler = true;
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/ringz/ringz.obj");
		readObjFile(Q, "../data/models/ringz/ringz.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		{
			double Ix, Iy, Iz;
			inertia_weight(Q, Ix, Iy, Iz);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = 1;
			std::cout << Ix << " " << Iy << " " << Iz << std::endl;
			distance_weight[3] = Ix; distance_weight[4] = Iy; distance_weight[5] = Iz;
		}

		std::vector<std::pair<C2A_Model*, Quaternion> > CSpace;

		std::ofstream timing_construct_file("timing_construct_exact_SE3.txt");
		boost::timer t_construct;
		int n_max = 200;

		std::ifstream rotation_stream("../data/ringz_exact_cspace/ringz_rotations.txt");

		int i = 0;
		while(!rotation_stream.eof())
		{
			std::string line;
			std::getline(rotation_stream, line, '\n');
			if(line.size() == 0) continue;
			std::istringstream reader(line);

			double a, b, c, d;
			reader >> a >> b >> c >> d;

			Quaternion q(a, b, c, d);

			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::string obj_file_name = std::string("../data/ringz_exact_cspace/ringz_") + ret + ".obj";

			C2A_Model* model;
			readObjFile(model, obj_file_name);

			CSpace.push_back(std::make_pair(model, q));

			i++;
			if(i >= n_max) break; 
		}

		timing_construct_file << t_construct.elapsed();
		timing_construct_file.flush();

		std::ofstream timing_file("timing_exact_SE3.txt");
		

		std::string base_name = "../data/models/Ringz/ConcaveTotusRainFallDemo/dump_transform";



		//for(std::size_t i = 750; i <= 950; ++i)
		for(std::size_t i = 771; i <= 950; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t j = 0; j < 4 - len; ++j)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::string PD_file_name= std::string("PD_exact_SE3") + ret + ".txt";


			std::ofstream PD_file(PD_file_name.c_str());

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile3D(filename, use_euler);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;

				if(frame[j][0].first != "8" || frame[j][1].first != "8") continue;

				boost::timer t;
				DataVector q = relative3D(frame[j][0].second, frame[j][1].second);
				if(collider.isCollide(q))
				{
					std::pair<DataVector, double> PD_result = Minkowski_Cspace_3D::Exact_PD_SE3(q, CSpace);
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

		std::string base_name = "../data/models/Ringz/ConcaveTotusRainFallDemo/dump_transform";

		bool use_euler = true;

		std::vector<C2A_Model*> P;
		std::vector<C2A_Model*> Q;

		readObjFiles(P, "../data/models/ringz/ringz_convex.obj");
		readObjFiles(Q, "../data/models/ringz/ringz_convex.obj");

		for(std::size_t i = 750; i <= 950; ++i)
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

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile3D(filename, use_euler);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;
				if(frame[j][0].first != "8" || frame[j][1].first != "8") continue;

				DataVector q = relative3D(frame[j][0].second, frame[j][1].second);

				boost::timer t;
				double pd = Collider3D::PDt(P, Q, q);
				timing_per_frame += t.elapsed();
				PD_file << pd << " ";
			}

			timing_file << timing_per_frame << " ";
			timing_file.flush();

			PD_file.flush();
		}

	}

	void playback()
	{
		std::vector<std::vector<DataVector> > frames;

		// std::string base_name = "../data/models/Ringz/ConcaveTotusRainFallDemo/dump_transform";
		std::string base_name = "../data/models/Ringz/dump_transform";

		bool use_euler = true;

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/ringz/ringz.obj");
		readObjFile(Q, "../data/models/ringz/ringz.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		{
			double Ix, Iy, Iz;
			inertia_weight(Q, Ix, Iy, Iz);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = 1;
			std::cout << Ix << " " << Iy << " " << Iz << std::endl;
			distance_weight[3] = Ix; distance_weight[4] = Iy; distance_weight[5] = Iz;
		}


		AABB3D aabb = computeAABB(P);

		// ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		double scale_x = aabb.b_max[0] - aabb.b_min[0];
		double scale_y = aabb.b_max[1] - aabb.b_min[1];
		double scale_z = aabb.b_max[2] - aabb.b_min[2];
		ContactSpaceSE3Euler2 contactspace(P, Q, scale_x, scale_y, scale_z);

		std::ofstream scaler_file("scaler_3d_rotation_ringz.txt");
		scaler_file << contactspace.getScaler() << std::endl;
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(100000);

		SVMLearner learner;
		learner.setDim(contactspace.active_data_dim());
		learner.setC(20);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50); 

		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");

		//std::vector<ContactSpaceSampleData> test_samples = contactspace.uniform_sample(20000);
		//std::cout << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << empiricalErrorRatio(test_samples, learner) << std::endl;

		// flann::HierarchicalClusteringIndex<ContactSpaceSE3Euler::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();

		std::vector<ContactSpaceSampleData> support_samples;
		learner.collectSupportVectors(support_samples);
		ExtendedModel<ContactSpaceSE3Euler2, flann::HierarchicalClusteringIndex> extended_model = 
			constructExtendedModelForModelDecisionBoundary<ContactSpaceSE3Euler2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(contactspace, learner, support_samples, 0.01, 50);


		std::ofstream timing_file("timing_APD.txt");

		// for(std::size_t i = 750; i <= 950; ++i)
		int i = 1050;
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

			std::vector<std::vector<std::pair<std::string, DataVector> > > frame = readDumpFile3D(filename, use_euler);

			double timing_per_frame = 0;

			for(std::size_t j = 0; j < frame.size(); ++j)
			{
				std::cout << i << " " << j << std::endl;

				if(frame[j][0].first != "8" || frame[j][1].first != "8") continue;

				boost::timer t;
				DataVector q = relative3D(frame[j][0].second, frame[j][1].second);
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

	void generateSamples()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Ringz/ringz.obj");
		readObjFile(Q, "../data/models/Ringz/ringz2.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_3d.txt");
		asciiWriter(out, contactspace_samples);
	}

	void playback_exact_CSpace_R3()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Ringz/ringz.obj");
		readObjFile(Q, "../data/models/Ringz/ringz2.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		C2A_Model* CSpace;
		readObjFile(CSpace, "../data/models/Ringz/ringzpair.obj");

		std::vector<ContactSpaceSampleData> contactspace_samples;
		std::ifstream in("space_test_3d.txt");
		asciiReader(in, contactspace_samples);

		std::ofstream timing_file("timing_exact_R3.txt");
		std::ofstream PD_file("PD_exact_R3.txt");

		for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		{
			std::cout << i << std::endl;
			DataVector q_col(6);
			DataVector q(3);
			for(std::size_t j = 0; j < 3; ++j)
				q[j] = contactspace_samples[i].v[j];

			for(std::size_t j = 0; j < 6; ++j)
				q_col[j] = contactspace_samples[i].v[j];


			boost::timer t;
			aTimer.Reset();
			aTimer.Start();
			std::pair<DataVector, double> pd_result;
			if(!collider.isCollide(q_col)) 
			{
				pd_result.second = 0;
			}
			else
			{
				pd_result = Minkowski_Cspace_3D::Exact_PD_R3(q, CSpace);
			}
			aTimer.Stop();
			PD_file << pd_result.second << " ";	
			//timing_file << t.elapsed() << " ";
			timing_file << aTimer.GetTime() * 1000 << " ";

			
			timing_file.flush();
			PD_file.flush();
		}
	}

	void playback_local_PD_R3()
	{
		std::ofstream timing_file("timing_local_PD_R3.txt");
		std::ofstream PD_file("PD_local_R3.txt");

		std::vector<C2A_Model*> P;
		std::vector<C2A_Model*> Q;

		readObjFiles(P, "../data/models/Ringz/ringz_convex.obj");
		readObjFiles(Q, "../data/models/Ringz/ringz2_convex.obj");


		std::vector<ContactSpaceSampleData> contactspace_samples;
		std::ifstream in("space_test_3d.txt");
		asciiReader(in, contactspace_samples);

		for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		{
			std::cout << i << std::endl;
			DataVector q_col(6);
			DataVector q(3);
			for(std::size_t j = 0; j < 3; ++j)
				q[j] = contactspace_samples[i].v[j];

			for(std::size_t j = 0; j < 6; ++j)
				q_col[j] = contactspace_samples[i].v[j];

			boost::timer t;
			aTimer.Reset();
			aTimer.Start();
			double pd = Collider3D::PDt(P, Q, q_col);
			aTimer.Stop();
			PD_file << pd << " ";	
			// timing_file << t.elapsed() << " ";
			timing_file << aTimer.GetTime() * 1000 << " ";

			
			timing_file.flush();
			PD_file.flush();
		}
	}

	void playback_R3()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Ringz/ringz.obj");
		readObjFile(Q, "../data/models/Ringz/ringz2.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		std::vector<ContactSpaceSampleData> contactspace_samples;
		std::ifstream in("space_test_3d.txt");
		asciiReader(in, contactspace_samples);

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_ringz.txt");
		scaler_file << contactspace.getScaler() << std::endl;
		std::vector<ContactSpaceSampleData> train_samples = contactspace.uniform_sample(100000);

		SVMLearner learner;
		learner.setDim(contactspace.active_data_dim());
		learner.setC(20);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50); 

		learner.learn(train_samples, contactspace.active_data_dim());
		learner.save("model_R3.txt");

		// flann::HierarchicalClusteringIndex<ContactSpaceSE3Euler::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();

		std::vector<ContactSpaceSampleData> support_samples;
		learner.collectSupportVectors(support_samples);
		ExtendedModel<ContactSpaceR3, flann::Index> extended_model = 
			constructExtendedModelForModelDecisionBoundary<ContactSpaceR3, SVMLearner, flann::Index, flann::KDTreeIndexParams>(contactspace, learner, support_samples, 0.01, 50);


		std::ofstream timing_file("timing_APD_R3.txt");
		std::ofstream PD_file("PD_APD_R3.txt");

		for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
		{
			std::cout << i << std::endl;
			DataVector q_col(6);
			DataVector q(3);
			for(std::size_t j = 0; j < 3; ++j)
				q[j] = contactspace_samples[i].v[j];

			for(std::size_t j = 0; j < 6; ++j)
				q_col[j] = contactspace_samples[i].v[j];

			boost::timer t;
			aTimer.Reset();
			aTimer.Start();
			if(!collider.isCollide(q_col)) 
			{
				PD_file << 0 << " ";
			}
			else
			{
				QueryResult pd_result = PD_query(learner, contactspace, extended_model.index, extended_model.samples, q);
				PD_file << pd_result.PD << " ";
			}
			aTimer.Stop();
			// timing_file << t.elapsed() << " ";
			timing_file << aTimer.GetTime() * 1000 << " ";

			
			timing_file.flush();
			PD_file.flush();
		}
	}
}

void main()
{
	// APDL::playback();
	// APDL::playback_local_PD();
	// APDL::playback_exact_CSpace();

	// APDL::generateSamples();
	// APDL::playback_exact_CSpace_R3();
	// APDL::playback_local_PD_R3();
	APDL::playback_R3();
}