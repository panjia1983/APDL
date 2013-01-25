#include <APD/minkowski_cspace.h>
#include <CGAL/IO/Polyhedron_iostream.h>

namespace APDL
{
	static void test_Minkowskiw_3D()
	{
		{
			// make tetrahedron
			Minkowski_Cspace_3D::Polyhedron P;
			Minkowski_Cspace_3D::Point_3 p(1, 0, 0);
			Minkowski_Cspace_3D::Point_3 q(0, 1, 0);
			Minkowski_Cspace_3D::Point_3 r(0, 0, 1);
			Minkowski_Cspace_3D::Point_3 s(0, 0, 0);
			P.make_tetrahedron(p, q, r, s);
			
			// make cube
			std::ifstream cube_file("../data/cube.off");
			
			if(!cube_file.is_open())
			{
				std::cerr << "Failed to open the cube file." << std::endl;
				return;
			}
			
			Minkowski_Cspace_3D::Polyhedron Q;
			cube_file >> Q;
			cube_file.close();
			
			
			std::cout << P << std::endl;
			std::cout << Q << std::endl;
			
			Minkowski_Cspace_3D::Polyhedron CSpace_R3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_R3(P, Q);
			
			
			DataVector query(3);
			query[0] = 0;
			query[1] = 0;
			query[2] = 0;
			std::pair<DataVector, double> PD_result = Minkowski_Cspace_3D::Exact_PD_R3(query, CSpace_R3);
			DataVector PD_point = PD_result.first;
			std::cout << PD_point[0] << " " << PD_point[1] << " " << PD_point[2] << " " << PD_result.second << std::endl;
			
			std::vector<std::pair<Minkowski_Cspace_3D::Polyhedron, Quaternion> > CSpace_SE3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_SE3(P, Q, 5);
			
			DataVector query2(6);
			query2[0] = 0;
			query2[1] = 0;
			query2[2] = 0;
			query2[3] = 0;
			query2[4] = 0;
			query2[5] = 0;
			
			std::pair<DataVector, double> PD_result2 = Minkowski_Cspace_3D::Exact_PD_SE3(query2, CSpace_SE3);
			DataVector PD_point2 = PD_result2.first;
			std::cout << "(" << PD_point2[0] << " " << PD_point2[1] << " " << PD_point2[2] << " "
			          << PD_point2[3] << " " << PD_point2[4] << " " << PD_point2[5] << ") " << PD_result2.second << std::endl;
			          
			          
			DataVector query3(7);
			query3[0] = 0;
			query3[1] = 0;
			query3[2] = 0;
			query3[3] = 1;
			query3[4] = 0;
			query3[5] = 0;
			query3[6] = 0;
			
			std::pair<DataVector, double> PD_result3 = Minkowski_Cspace_3D::Exact_PD_SE3(query3, CSpace_SE3);
			DataVector PD_point3 = PD_result3.first;
			std::cout << "(" << PD_point3[0] << " " << PD_point3[1] << " " << PD_point3[2] << " "
			          << PD_point3[3] << " " << PD_point3[4] << " " << PD_point3[5] << " " << PD_point3[6] << ") " << PD_result3.second << std::endl;
		}
		
		

		
		{
			Minkowski_Cspace_3D::Polyhedron P, Q;
			
			Minkowski_Cspace_3D::Point_3 p(1, 0, 0);
			Minkowski_Cspace_3D::Point_3 q(0, 1, 0);
			Minkowski_Cspace_3D::Point_3 r(0, 0, 1);
			Minkowski_Cspace_3D::Point_3 s(0, 0, 0);
			P.make_tetrahedron(p, q, r, s);
			
			std::ifstream spoon_file("../data/spoon.off");
			
			if(!spoon_file.is_open())
			{
				std::cerr << "Failed to open the spoon file." << std::endl;
				return;
			}
			
			spoon_file >> Q;
			spoon_file.close();
			
			
			Minkowski_Cspace_3D::Polyhedron CSpace_R3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_R3(P, Q);
		}


        ////// Expensive
		//{
		//	Minkowski_Cspace_3D::Polyhedron P, Q;
		//	std::ifstream cup_file("../data/cup.off");
		
		//	if(!cup_file.is_open())
		//	{
		//		std::cerr << "Failed to open the cup file." << std::endl;
		//		return;
		//	}
		
		//	cup_file >> P;
		//	cup_file.close();
		
		//	std::ifstream spoon_file("../data/spoon.off");
		
		//	if(!spoon_file.is_open())
		//	{
		//		std::cerr << "Failed to open the spoon file." << std::endl;
		//		return;
		//	}
		
		//	spoon_file >> Q;
		//	spoon_file.close();
		
		//	std::cout << P.size_of_vertices() << " " << Q.size_of_vertices() << std::endl;
		
		//	Minkowski_Cspace_3D::Polyhedron CSpace_R3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_R3(P, Q);
		//}
		
	}

	static void test_Minkowskiw_3D_cupspoon()
	{
		std::ifstream cup_file("../data/cup.off");

		if(!cup_file.is_open())
		{
			std::cerr << "Failed to open the cup file." << std::endl;
			return;
		}

		std::ifstream spoon_file("../data/spoon.off");

		if(!spoon_file.is_open())
		{
			std::cerr << "Failed to open the spoon file." << std::endl;
			return;
		}

		Minkowski_Cspace_3D::Polyhedron P;
		Minkowski_Cspace_3D::Polyhedron Q;
		cup_file >> P;
		spoon_file >> Q;

		Minkowski_Cspace_3D::Polyhedron CSpace_R3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_R3(P, Q);

		std::ofstream cupspoon_file("cupspoon.off");
		cupspoon_file << CSpace_R3;
	}
}

void main()
{
	// APDL::test_Minkowskiw_3D();
	APDL::test_Minkowskiw_3D_cupspoon();
}
