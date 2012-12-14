#include <APD/minkowski_cspace.h>

namespace APDL
{

	static void test_Minkowskiw_2D()
	{
		{
			Minkowski_Cspace_2D::Polygon_2 P;
			P.push_back(Minkowski_Cspace_2D::Point_2(0, 0));
			P.push_back(Minkowski_Cspace_2D::Point_2(6, 0));
			P.push_back(Minkowski_Cspace_2D::Point_2(3, 5));
			
			Minkowski_Cspace_2D::Polygon_2 Q;
			Q.push_back(Minkowski_Cspace_2D::Point_2(0, 0));
			Q.push_back(Minkowski_Cspace_2D::Point_2(2, -2));
			Q.push_back(Minkowski_Cspace_2D::Point_2(2, 2));
			
			Minkowski_Cspace_2D::Polygon_with_holes_2 cspace_R2 = Minkowski_Cspace_2D::Minkowski_Cobstacle_R2(P, Q);
			
			assert(cspace_R2.number_of_holes() == 0);
			
			std::cout << "P = ";
			print_polygon(P);
			std::cout << "Q = ";
			print_polygon(Q);
			std::cout << "P - Q = ";
			print_polygon(cspace_R2.outer_boundary());
			
			Minkowski_Cspace_2D::Point_2 O(0, 0);
			
			assert(cspace_R2.outer_boundary().oriented_side(O) == CGAL::ON_POSITIVE_SIDE);
			
			DataVector query(2);
			query[0] = 0;
			query[1] = 0;
			
			std::pair<DataVector, double> PD_result = Minkowski_Cspace_2D::Exact_PD_R2(query, cspace_R2);
			DataVector PD_point = PD_result.first;
			std::cout << PD_point[0] << " " << PD_point[1] << " " << PD_result.second << std::endl;
			
			std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> > cspace_SE2 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(P, Q, 30);
			
			DataVector query2(3);
			query2[0] = 0;
			query2[1] = 0;
			query2[2] = 0;
			std::pair<DataVector, double> PD_result2 = Minkowski_Cspace_2D::Exact_PD_SE2(query2, cspace_SE2, 1);
			DataVector PD_point2 = PD_result2.first;
			std::cout << PD_point2[0] << " " << PD_point2[1] << " " << PD_point2[2] << " " << PD_result2.second << std::endl;
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
			
			Minkowski_Cspace_2D::Polygon_with_holes_2 cspace_R2 = Minkowski_Cspace_2D::Minkowski_Cobstacle_R2(P, Q);
			
			std::cout << "P - Q = ";
			print_polygon_with_holes(cspace_R2);
			
			DataVector query(2);
			query[0] = 0;
			query[1] = 0;
			
			std::pair<DataVector, double> PD_result = Minkowski_Cspace_2D::Exact_PD_R2(query, cspace_R2);
			DataVector PD_point = PD_result.first;
			std::cout << PD_point[0] << " " << PD_point[1] << " " << PD_result.second << std::endl;
			
			std::vector<std::pair<Minkowski_Cspace_2D::Polygon_with_holes_2, double> > cspace_SE2 = Minkowski_Cspace_2D::Minkowski_CObstacle_SE2(P, Q, 30);
			
			DataVector query2(3);
			query2[0] = 0;
			query2[1] = 0;
			query2[2] = 0;
			std::pair<DataVector, double> PD_result2 = Minkowski_Cspace_2D::Exact_PD_SE2(query2, cspace_SE2, 1);
			DataVector PD_point2 = PD_result2.first;
			std::cout << PD_point2[0] << " " << PD_point2[1] << " " << PD_point2[2] << " " << PD_result2.second << std::endl;
		}
	}
}

void main()
{
	APDL::test_Minkowskiw_2D();
}