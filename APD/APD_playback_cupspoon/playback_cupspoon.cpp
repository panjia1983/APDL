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
		//{
		//	Quaternion q(1, 3, 4, 5);
		//	for(int i = 0; i < 4; ++i)
		//		std::cout << q[i] << std::endl;
		//	double R[3][3];
		//	Quat2Rot(R, q);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}
		//	Rot2Quat(q, R);
		//	for(int i = 0; i < 4; ++i)
		//		std::cout << q[i] << std::endl;
		//	Quat2Rot(R, q);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}

		//	return;

		//}
		//{
		//	double a = 1, b = 2, c = 3;
		//	double R[3][3];
		//	Euler2Rot(R, a, b, c);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}
		//	Rot2Euler(a, b, c, R);
		//	std::cout << a << " " << b << " " << c << std::endl;

		//	Euler2Rot(R, a, b, c);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}
		//}

		//{
		//	Quaternion q(1, 2, 3, 4);
		//	double a, b, c;
		//	Quat2Euler(a, b, c, q);
		//	double R[3][3];
		//	Quat2Rot(R, q);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}

		//	Euler2Rot(R, a, b, c);
		//	for(int i = 0; i < 3; ++i)
		//	{
		//		for(int j = 0; j < 3; ++j)
		//			std::cout << R[i][j] << " ";
		//		std::cout << std::endl;
		//	}
		//}

		{
			double R[3][3];
			R[0][0] = 0.0;
			R[0][1] = -1.0;
			R[0][2] = 0.0;
			R[1][0] = 0.984808;
			R[1][1] = 0.0;
			R[1][2] = 0.173648;
			R[2][0] = -0.173648;
			R[2][1] = 0.0;
			R[2][2] = 0.984808;

			double a, b, c;
			Rot2Euler(a, b, c, R);
			Quaternion q;
			Euler2Quat(q, a, b, c);

			Euler2Rot(R, a, b, c);
			for(int j = 0; j < 3; ++j)
			{
				for(int k = 0; k < 3; ++k)
					std::cout << R[j][k] << " ";
				std::cout << std::endl;
			}

			Quat2Rot(R, q);
			for(int j = 0; j < 3; ++j)
			{
				for(int k = 0; k < 3; ++k)
					std::cout << R[j][k] << " ";
				std::cout << std::endl;
			}
		}
		return;
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
			if(use_euler)
			{
				distance_weight[3] = Ix; distance_weight[4] = Iy; distance_weight[5] = Iz;
			}
			else
			{
				distance_weight[3] = 1; distance_weight[4] = Ix; distance_weight[5] = Iy; distance_weight[6] = Iz;
			}
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

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon.ani";

		frames = readAnimationFile(frame_file_name, use_euler);

		std::vector<DataVector> result_qs;

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			double R[3][3];
			double T[3];

			if(q.dim() == 6)
			{
				ConfigEuler2RotTran(R, T, q);
			}
			else if(q.dim() == 7)
			{
				ConfigQuat2RotTrans(R, T, q);
			}

			for(int j = 0; j < 3; ++j)
				for(int k = 0; k < 3; ++k)
					std::cout << R[j][k] << " ";
			for(int j = 0; j < 3; ++j)
				std::cout << T[j] << " ";
			std::cout << std::endl;

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

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon.ani";

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

		std::string frame_file_name = "../data/models/CupSpoon/CupSpoon.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			for(int j = 0; j < q.dim(); ++j)
				std::cout << q[j] << " ";
			std::cout << std::endl;
		}

	}
}

void main()
{
	// APDL::playback();
	// APDL::playback_local_PD();
	APDL::playback_exact_CSpace();
}