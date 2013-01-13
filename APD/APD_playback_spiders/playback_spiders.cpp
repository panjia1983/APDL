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


	void playback()
	{
		std::vector<std::vector<std::pair<std::string, DataVector> > > frames;

		std::string base_name = "../data/models/Box2D/spider_cinfigs/dump_transform";

		// for(std::size_t i = 1001; i < 9999; ++i)
		for(std::size_t i = 1001; i < 1101; ++i)
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

		std::vector<Polygon> polys1 = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);
		std::vector<Polygon> polys2 = readPolyFile("../data/models/Box2D/nazca_spider77.polys", 3);

		ContactSpaceSE2 contactspace(polys1, polys1, 0.2 * (getCircle(polys1).second + getCircle(polys1).second));

		//ofstream test_file("test.txt");

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			if(frames[i][0].first == "Wall" || frames[i][1].first == "Wall") continue;

			//test_file << frames[i][0].first << " ";
			//for(std::size_t j = 0; j < frames[i][0].second.dim(); ++j)
			//	test_file << frames[i][0].second[j] << " ";
			//test_file << std::endl;

			//test_file << frames[i][1].first << " ";
			//for(std::size_t j = 0; j < frames[i][1].second.dim(); ++j)
			//	test_file << frames[i][1].second[j] << " ";
			//test_file << std::endl;

			//test_file << std::endl;

			DataVector q = relative2D(frames[i][0].second, frames[i][1].second);

			//std::cout << q[0] << " " << q[1] << " " << q[2] << std::endl;

		}

	}
}

void main()
{
	// APDL::playback();

	// APDL::playback_local_PD();
	APDL::playback_exact_CSpace();
}