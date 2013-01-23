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
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

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

		int n_max = 200;

		std::ifstream rotation_stream("../data/cupspoon_exact_cspace/cupspoon_rotations.txt");

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
			std::string obj_file_name = std::string("../data/cupspoon_exact_cspace/cupspoon_") + ret + ".obj";

			C2A_Model* model;
			readObjFile(model, obj_file_name);

			CSpace.push_back(std::make_pair(model, q));

			i++;
			if(i >= n_max) break; 
		}

		std::ofstream timing_file("timing_exact_SE2.txt");
		std::ofstream PD_file("PD_exact_SE2.txt");

		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon_NoJump.ani";

		frames = readAnimationFile(frame_file_name, use_euler);

		std::vector<DataVector> result_qs;

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			//double R[3][3];
			//double T[3];

			//if(q.dim() == 6)
			//{
			//	ConfigEuler2RotTran(R, T, q);
			//}
			//else if(q.dim() == 7)
			//{
			//	ConfigQuat2RotTrans(R, T, q);
			//}

			//for(int j = 0; j < 3; ++j)
			//	for(int k = 0; k < 3; ++k)
			//		std::cout << R[j][k] << " ";
			//for(int j = 0; j < 3; ++j)
			//	std::cout << T[j] << " ";
			//std::cout << std::endl;

			boost::timer t;
			std::pair<DataVector, double> pd_result;
			if(!collider.isCollide(q)) 
			{
				pd_result.first = frames[i][1];
				
				//DataVector q2 = unrelative3D(frames[i][0], q);
				//for(int j = 0; j < q2.dim(); ++j)
				//	std::cout << q2[j] << " ";
				//std::cout << std::endl;
				//for(int j = 0; j < q2.dim(); ++j)
				//	std::cout << frames[i][1][j] << " ";
				//std::cout << std::endl;

				pd_result.second = 0;
				result_qs.push_back(frames[i][1]);
			}
			else
			{
				pd_result = Minkowski_Cspace_3D::Exact_PD_SE3(q, CSpace);
				DataVector q2 = unrelative3D(frames[i][0], pd_result.first);
				result_qs.push_back(q2);
			}
			PD_file << pd_result.second << " ";	
			timing_file << t.elapsed() << " ";

			
			timing_file.flush();
			PD_file.flush();
		}

		std::ofstream new_animation_stream("new_animation.ani");
		new_animation_stream << frames.size() << "f" << std::endl;
		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			new_animation_stream << i << "f" << std::endl;

			double R[3][3];
			double T[3];

			if(frames[i][0].dim() == 6)
			{
				ConfigEuler2RotTran(R, T, frames[i][0]);
			}
			else if(frames[i][0].dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, frames[i][0]);
			}

			for(int j = 0; j < 3; ++j)
				for(int k = 0; k < 3; ++k)
				{
					// new_animation_stream << R[j][k] << " ";
					new_animation_stream << R[k][j] << " ";
				}
			for(int j = 0; j < 3; ++j)
				new_animation_stream << T[j] << " ";
			new_animation_stream << std::endl;

			if(result_qs[i].dim() == 6)
			{
				ConfigEuler2RotTran(R, T, result_qs[i]);
			}
			else if(result_qs[i].dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, result_qs[i]);
			}

			for(int j = 0; j < 3; ++j)
				for(int k = 0; k < 3; ++k)
				{
					// new_animation_stream << R[j][k] << " ";
					new_animation_stream << R[k][j] << " ";
				}
			for(int j = 0; j < 3; ++j)
				new_animation_stream << T[j] << " ";
			new_animation_stream << std::endl;
		}
	}

	void playback_local_PD()
	{
		std::ofstream timing_file("timing_local_PD.txt");
		std::ofstream PD_file("PD_local.txt");

		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon_NoJump.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		std::vector<C2A_Model*> P;
		std::vector<C2A_Model*> Q;

		readObjFiles(P, "../data/models/CupSpoon/cup_convex.obj");
		readObjFiles(Q, "../data/models/CupSpoon/spoon_convex.obj");

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			boost::timer t;
			double pd = Collider3D::PDt(P, Q, q);
			PD_file << pd << " ";	
			timing_file << t.elapsed() << " ";

			timing_file.flush();
			PD_file.flush();
		}

	}

	void playback()
	{
		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon_NoJump.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

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


		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_cupspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(500000);

		SVMLearner learner;
		learner.setDim(contactspace.active_data_dim());
		learner.setC(20);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50); 

		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");

		// flann::HierarchicalClusteringIndex<ContactSpaceSE3Euler::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();

		std::vector<ContactSpaceSampleData> support_samples;
		learner.collectSupportVectors(support_samples);
		ExtendedModel<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex> extended_model = 
			constructExtendedModelForModelDecisionBoundary<ContactSpaceSE3Euler, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(contactspace, learner, support_samples, 0.01, 50);


		std::ofstream timing_file("timing_APD.txt");
		std::ofstream PD_file("PD_APD.txt");

		std::vector<DataVector> result_qs;

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			boost::timer t;
			if(!collider.isCollide(q)) 
			{
				result_qs.push_back(frames[i][1]);
				PD_file << 0 << " ";
			}
			else
			{
				// QueryResult pd_result = PD_query(learner, contactspace, query_index, q);
				QueryResult pd_result = PD_query(learner, contactspace, extended_model.index, extended_model.samples, q);
				DataVector q2 = unrelative3D(frames[i][0], pd_result.v);
				result_qs.push_back(q2);
				PD_file << pd_result.PD << " ";
			}

			timing_file << t.elapsed() << " ";

			
			timing_file.flush();
			PD_file.flush();
		}

		std::ofstream new_animation_stream("new_animation_APD.ani");
		new_animation_stream << frames.size() << "f" << std::endl;
		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			new_animation_stream << i << "f" << std::endl;

			double R[3][3];
			double T[3];

			if(frames[i][0].dim() == 6)
			{
				ConfigEuler2RotTran(R, T, frames[i][0]);
			}
			else if(frames[i][0].dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, frames[i][0]);
			}

			for(int j = 0; j < 3; ++j)
				for(int k = 0; k < 3; ++k)
				{
					// new_animation_stream << R[j][k] << " ";
					new_animation_stream << R[k][j] << " ";
				}
			for(int j = 0; j < 3; ++j)
				new_animation_stream << T[j] << " ";
			new_animation_stream << std::endl;

			if(result_qs[i].dim() == 6)
			{
				ConfigEuler2RotTran(R, T, result_qs[i]);
			}
			else if(result_qs[i].dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, result_qs[i]);
			}

			for(int j = 0; j < 3; ++j)
				for(int k = 0; k < 3; ++k)
				{
					// new_animation_stream << R[j][k] << " ";
					new_animation_stream << R[k][j] << " ";
				}
			for(int j = 0; j < 3; ++j)
				new_animation_stream << T[j] << " ";
			new_animation_stream << std::endl;
		}

	}






	void generateSamples()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

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
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		C2A_Model* CSpace;
		readObjFile(CSpace, "../data/cupspoon.obj");

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
			//aTimer.Reset();
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

		readObjFiles(P, "../data/models/CupSpoon/cup_convex.obj");
		readObjFiles(Q, "../data/models/CupSpoon/spoon_convex.obj");


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
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		std::vector<ContactSpaceSampleData> contactspace_samples;
		std::ifstream in("space_test_3d.txt");
		asciiReader(in, contactspace_samples);

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_cupspoon.txt");
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
	APDL::playback_exact_CSpace_R3();
	// APDL::playback_local_PD_R3();
	// APDL::playback_R3();

}