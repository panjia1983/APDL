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


	void playback()
	{
		std::vector<std::vector<std::pair<std::string, DataVector> > > frames;

		std::string base_name = "../data/models/Box2D/angrybird_configs/dump_transform";

		// for(std::size_t i = 0; i < 999; ++i)
		for(std::size_t i = 0; i < 10; ++i)
		{
			std::stringstream ss;
			ss << i;
			std::string ret;
			ss >> ret;
			std::size_t len = ret.length();
			for(std::size_t i = 0; i < 4 - len; ++i)
				ret = "0" + ret;

			std::string filename = base_name + ret + ".txt";

			std::vector<std::vector<std::pair<std::string, DataVector> > > frames_ = readDumpFile(filename);

			for(std::size_t i = 0; i < frames_.size(); ++i)
				frames.push_back(frames_[i]);
		}

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/GreenPig32.polys");
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/BigRedBird17.polys");
		std::vector<Polygon> polys3 = readPolyFile("../data/models/Box2D/WhiteBird30.polys");


		ContactSpaceSE2 contactspace11(polys1, polys1, 0.2 * (getCircle(polys1).second + getCircle(polys1).second));
		ContactSpaceSE2 contactspace12(polys1, polys2, 0.2 * (getCircle(polys1).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace13(polys1, polys3, 0.2 * (getCircle(polys1).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace22(polys2, polys2, 0.2 * (getCircle(polys2).second + getCircle(polys2).second));
		ContactSpaceSE2 contactspace23(polys2, polys3, 0.2 * (getCircle(polys2).second + getCircle(polys3).second));
		ContactSpaceSE2 contactspace33(polys3, polys3, 0.2 * (getCircle(polys3).second + getCircle(polys3).second));

		// ofstream test_file("test.txt");

		for(std::size_t i = 0; i < frames.size(); ++i)
		{

			//test_file << frames[i][0].first << " ";
			//for(std::size_t j = 0; j < frames[i][0].second.dim(); ++j)
			//	test_file << frames[i][0].second[j] << " ";
			//test_file << std::endl;

			//test_file << frames[i][1].first << " ";
			//for(std::size_t j = 0; j < frames[i][1].second.dim(); ++j)
			//	test_file << frames[i][1].second[j] << " ";
			//test_file << std::endl;

			//test_file << std::endl;
			
			std::vector<Polygon>* obj1 = NULL;
			std::vector<Polygon>* obj2 = NULL;

			if(frames[i][0].first == "GreenPig") obj1 = &polys1;
			else if(frames[i][0].first == "RedBird") obj1 = &polys2;
			else if(frames[i][0].first == "WhiteBird") obj1 = &polys3;
			else 
			{
				std::cout << "error" << std::endl;
				return;
			}

			if(frames[i][1].first == "GreenPig") obj2 = &polys1;
			else if(frames[i][1].first == "RedBird") obj2 = &polys2;
			else if(frames[i][1].first == "WhiteBird") obj2 = &polys3;
			else 
			{
				std::cout << "error" << std::endl;
				return;
			}

			DataVector q = relative2D(frames[i][0].second, frames[i][1].second);

			// std::cout << q[0] << " " << q[1] << " " << q[2] << std::endl;

		}

	}
}

void main()
{
	// APDL::playback();
	APDL::playback_exact_CSpace();
	// APDL::playback_local_PD();
}