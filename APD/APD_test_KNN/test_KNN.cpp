#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>

namespace APDL
{	
	void test_KNN()
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
			
			ContactSpaceR2 contactspace(p1, p2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

			flann::Index<ContactSpaceR2::DistanceType>* index = NULL;

			generateIndex<ContactSpaceR2::DistanceType>(contactspace_samples, 
				contactspace.active_data_dim(), 
				index,
				flann::KDTreeIndexParams());

		    std::vector<std::vector<int> > indices;
		    std::vector<std::vector<double> > dists;

			knnSearch<ContactSpaceR2::DistanceType>(contactspace_samples, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   10,
				   flann::SearchParams());

			indices.clear();
			dists.clear();

			radiusSearch<ContactSpaceR2::DistanceType>(contactspace_samples, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   0.5,
				   flann::SearchParams());

			std::ofstream indices_file("indices.txt");
			std::ofstream dists_file("dists.txt");

			indices_file << indices; dists_file << dists;
			indices_file.close(); dists_file.close();

			delete index;
		}

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
			
			ContactSpaceSE2 contactspace(p1, p2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(100000);

			flann::Index<ContactSpaceR2::DistanceType>* index = NULL;

			generateIndex<ContactSpaceR2::DistanceType>(contactspace_samples, 
				contactspace.active_data_dim(), 
				index,
				flann::KDTreeIndexParams());

		    std::vector<std::vector<int> > indices;
		    std::vector<std::vector<double> > dists;

			knnSearch<ContactSpaceR2::DistanceType>(contactspace_samples, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   10,
				   flann::SearchParams());

			indices.clear();
			dists.clear();

			radiusSearch<ContactSpaceR2::DistanceType>(contactspace_samples, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   0.5,
				   flann::SearchParams());

			std::ofstream indices_file("indices2.txt");
			std::ofstream dists_file("dists2.txt");

			indices_file << indices; dists_file << dists;
			indices_file.close(); dists_file.close();

			delete index;
		}
	}
}

void main()
{
	APDL::test_KNN();
}