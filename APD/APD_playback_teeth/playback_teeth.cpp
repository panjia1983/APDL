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
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

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

		int n_max = 100;

		std::ifstream rotation_stream("../data/teeth_exact_cspace/teeth_rotations.txt");

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
			std::string obj_file_name = std::string("../data/teeth_exact_cspace/teeth_") + ret + ".obj";

			C2A_Model* model;
			readObjFile(model, obj_file_name);

			CSpace.push_back(std::make_pair(model, q));

			i++;
			if(i >= n_max) break; 
		}

		std::ofstream timing_file("timing_exact_SE3.txt");
		std::ofstream PD_file("PD_exact_SE3.txt");

		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/Teeth/teeth0116.ani";

		frames = readAnimationFile(frame_file_name, use_euler);

		std::vector<DataVector> result_qs;

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			boost::timer t;
			std::pair<DataVector, double> pd_result;
			if(!collider.isCollide(q)) 
			{
				pd_result.first = frames[i][1];
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
					new_animation_stream << R[j][k] << " ";
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
					new_animation_stream << R[j][k] << " ";
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

		std::string frame_file_name = "../data/models/Teeth/teeth0116.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		std::vector<C2A_Model*> P;
		std::vector<C2A_Model*> Q;

		readObjFiles(P, "../data/models/Teeth/lo_03de_0116_a1_convex2.obj");
		readObjFiles(Q, "../data/models/Teeth/up_03de_0116_a1_convex2.obj");

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);


			double R[3][3];
			double T[3];

			if(frames[i][0].dim() == 6)
			{
				ConfigEuler2RotTran(R, T, q);
			}
			else if(frames[i][0].dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, q);
			}

			boost::timer t;
			aTimer.Reset();
			aTimer.Start();
			double pd = Collider3D::PDt(P, Q, q);
			aTimer.Stop();
			PD_file << pd << " ";	
			//timing_file << t.elapsed() << " ";
			timing_file << aTimer.GetTime() * 1000 << " ";

			timing_file.flush();
			PD_file.flush();
		}

	}

	void playback()
	{
		double dynamic_weight = 10;

		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/Teeth/teeth0116.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_0116.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_0116.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		Collider3D collider(P, Q);

		double distance_weight_backup[7];

		{
			double Ix, Iy, Iz;
			inertia_weight(Q, Ix, Iy, Iz);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = 1;
			std::cout << Ix << " " << Iy << " " << Iz << std::endl;
			distance_weight[3] = Ix; distance_weight[4] = Iy; distance_weight[5] = Iz;

			distance_weight[5] *= (dynamic_weight * 10);

			for(int i = 0; i < 7; ++i)
				distance_weight_backup[i] = distance_weight[i];
		}

		AABB3D aabb = computeAABB(P);

		//ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		double scale_x = aabb.b_max[0] - aabb.b_min[0];
		double scale_y = aabb.b_max[1] - aabb.b_min[1];
		double scale_z = aabb.b_max[2] - aabb.b_min[2];
		ContactSpaceSE3Euler2 contactspace(P, Q, scale_x, scale_y, scale_z);
		
		std::ofstream scaler_file("scaler_3d_rotation_teeth.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		int n_samples = 50000;
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(n_samples);

		SVMLearner learner;
		learner.setDim(contactspace.active_data_dim());
		learner.setC(20);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50); 

		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");

		// flann::HierarchicalClusteringIndex<ContactSpaceSE3Euler::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();

		int knn_k = 50;
		std::vector<ContactSpaceSampleData> support_samples;
		learner.collectSupportVectors(support_samples);
		//ExtendedModel<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex> extended_model = 
		//	constructExtendedModelForModelDecisionBoundary<ContactSpaceSE3Euler, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(contactspace, learner, support_samples, 0.01, knn_k);
		ExtendedModel<ContactSpaceSE3Euler2, flann::HierarchicalClusteringIndex> extended_model = 
			constructExtendedModelForModelDecisionBoundary<ContactSpaceSE3Euler2, SVMLearner, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(contactspace, learner, support_samples, 0.01, knn_k);


		distance_weight[3] = distance_weight_backup[3];
		distance_weight[4] = distance_weight_backup[4] * dynamic_weight;
		flann::HierarchicalClusteringIndex<typename ContactSpaceSE3Euler::DistanceType>* index1 = constructIndexForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(extended_model.samples);


		distance_weight[3] = distance_weight_backup[3] * dynamic_weight;
		distance_weight[4] = distance_weight_backup[4];
		flann::HierarchicalClusteringIndex<typename ContactSpaceSE3Euler::DistanceType>* index2 = constructIndexForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>(extended_model.samples);




		std::ofstream timing_file("timing_APD.txt");
		std::ofstream PD_file("PD_APD.txt");

		std::vector<DataVector> result_qs;

		
		
		for(std::size_t i = 0; i < frames.size(); ++i)
		{

			flann::HierarchicalClusteringIndex<typename ContactSpaceSE3Euler::DistanceType>* index = NULL;
			//if(i < 102) 
			//{
			//	index = index1;
			//	distance_weight[3] = distance_weight_backup[3];
			//	distance_weight[4] = distance_weight_backup[4] * dynamic_weight;
			//}
			//else if(i < 251) 
			//{
			//	index = index2;
			//	distance_weight[3] = distance_weight_backup[3] * dynamic_weight;
			//	distance_weight[4] = distance_weight_backup[4];
			//}
			//else if(i < 300)
			{
				index = extended_model.index;
				distance_weight[3] = distance_weight_backup[3];
				distance_weight[4] = distance_weight_backup[4];
			}
			//else if(i < 351)
			//{
			//	index = index1;
			//	distance_weight[3] = distance_weight_backup[3];
			//	distance_weight[4] = distance_weight_backup[4] * dynamic_weight;
			//}
			//else if(i < 700)
			//{
			//	index = index2;
			//	distance_weight[3] = distance_weight_backup[3] * dynamic_weight;
			//	distance_weight[4] = distance_weight_backup[4];
			//}
			//else if(i < 751)
			//{
			//	index = index1;
			//	distance_weight[3] = distance_weight_backup[3];
			//	distance_weight[4] = distance_weight_backup[4] * dynamic_weight;
			//}
			//else if(i < 852)
			//{
			//	index = index2;
			//	distance_weight[3] = distance_weight_backup[3] * dynamic_weight;
			//	distance_weight[4] = distance_weight_backup[4];
			//}
			//else
			//{
			//	index = extended_model.index;
			//	distance_weight[3] = distance_weight_backup[3];
			//	distance_weight[4] = distance_weight_backup[4];
			//}


			DataVector q = relative3D(frames[i][0], frames[i][1]);

			boost::timer t;
			aTimer.Reset();
			aTimer.Start();
			if(!collider.isCollide(q)) 
			{
				result_qs.push_back(frames[i][1]);
				//PD_file << 0 << " ";
			}
			else
			{
				// QueryResult pd_result = PD_query(learner, contactspace, query_index, q);
				QueryResult pd_result = PD_query(learner, contactspace, index, extended_model.samples, q);
				DataVector q2 = unrelative3D(frames[i][0], pd_result.v);
				result_qs.push_back(q2);
				//PD_file << pd_result.PD << " ";
			}

			aTimer.Stop();
			//timing_file << t.elapsed() << " ";
			timing_file << aTimer.GetTime() * 1000 << " ";


			
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

	void checkEulerAngle()
	{
		std::vector<std::vector<DataVector> > frames;

		// std::string frame_file_name = "../data/models/Teeth/teeth0116.ani";
		std::string frame_file_name = "new_animation_APD.ani";

		std::ofstream file("test_euler_angle.txt");

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			for(std::size_t j = 0; j < frames[i][1].dim(); ++j)
				file << frames[i][1][j] << " ";
			file << std::endl;
		}
	}
}

void main()
{
	APDL::playback();
	// APDL::playback_local_PD();
	// APDL::playback_exact_CSpace();
	// APDL::checkEulerAngle();
}