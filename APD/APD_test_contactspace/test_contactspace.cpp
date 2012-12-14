#include <APD/minkowski_cspace.h>
#include <APD/contact_space.h>
#include <APD/mesh_io.h>

namespace APDL
{
	void test_ContactSpace()
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
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);
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
			
			std::pair<Vec2D, double> c2 = p2.getCircle();
			
			ContactSpaceSE2 contactspace(p1, p2, 0.1 * c2.second);
			for(int i = 0; i < 100000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d_2.txt");
			asciiWriter(out, contactspace);
		}
		
		{
			C2A_Model* model1 = new C2A_Model;
			C2A_Model* model2 = new C2A_Model;
			readOffFile(model1, "../data/cup.off");
			readOffFile(model2, "../data/spoon.off");
			ContactSpaceR3 contactspace(model1, model2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test.txt");
			asciiWriter(out, contactspace);
			
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			//std::cout << model1->radius << " " << model2->radius << std::endl;
			//std::cout << model1->com[0] << " " << model1->com[1] << " " << model1->com[2] << std::endl;
			//std::cout << model2->com[0] << " " << model2->com[1] << " " << model2->com[2] << std::endl;
		}
		
		
		{
			C2A_Model* model1 = new C2A_Model;
			C2A_Model* model2 = new C2A_Model;
			readObjFile(model1, "../data/cup.obj");
			readObjFile(model2, "../data/spoon.obj");
			ContactSpaceR3 contactspace(model1, model2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test.txt");
			asciiWriter(out, contactspace);
			
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			//std::cout << model1->radius << " " << model2->radius << std::endl;
			//std::cout << model1->com[0] << " " << model1->com[1] << " " << model1->com[2] << std::endl;
			//std::cout << model2->com[0] << " " << model2->com[1] << " " << model2->com[2] << std::endl;
		}
		
		{
			C2A_Model* model1 = new C2A_Model;
			C2A_Model* model2 = new C2A_Model;
			readTriFile(model1, "../data/cup.tri");
			readTriFile(model2, "../data/spoon.tri");
			ContactSpaceR3 contactspace(model1, model2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test.txt");
			asciiWriter(out, contactspace);
			
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			//std::cout << model1->radius << " " << model2->radius << std::endl;
			//std::cout << model1->com[0] << " " << model1->com[1] << " " << model1->com[2] << std::endl;
			//std::cout << model2->com[0] << " " << model2->com[1] << " " << model2->com[2] << std::endl;
		}
		
		{
			C2A_Model* model1 = new C2A_Model;
			C2A_Model* model2 = new C2A_Model;
			readTriFile(model1, "../data/cup.tri");
			readTriFile(model2, "../data/spoon.tri");
			ContactSpaceSE3Quat contactspace(model1, model2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test2.txt");
			asciiWriter(out, contactspace);
			
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			//std::cout << model1->radius << " " << model2->radius << std::endl;
			//std::cout << model1->com[0] << " " << model1->com[1] << " " << model1->com[2] << std::endl;
			//std::cout << model2->com[0] << " " << model2->com[1] << " " << model2->com[2] << std::endl;
		}
		
		{
			C2A_Model* model1 = new C2A_Model;
			C2A_Model* model2 = new C2A_Model;
			readTriFile(model1, "../data/cup.tri");
			readTriFile(model2, "../data/spoon.tri");
			ContactSpaceSE3Quat2 contactspace(model1, model2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test3.txt");
			asciiWriter(out, contactspace);
			
			model1->ComputeCenterOfMass();
			model1->ComputeRadius();
			model2->ComputeCenterOfMass();
			model2->ComputeRadius();
			//std::cout << model1->radius << " " << model2->radius << std::endl;
			//std::cout << model1->com[0] << " " << model1->com[1] << " " << model1->com[2] << std::endl;
			//std::cout << model2->com[0] << " " << model2->com[1] << " " << model2->com[2] << std::endl;
			
			out.close();
			
			std::ifstream in("space_test3.txt");
			contactspace.data.clear();
			asciiReader(in, contactspace);
			in.close();
			
			std::ofstream out2("space_test4.txt");
			asciiWriter(out2, contactspace);
			out2.close();
			
			std::ofstream out3("space_test5.txt", std::ios::binary);
			binaryWriter(out3, contactspace);
			out3.close();
			
			std::ifstream in2("space_test5.txt", std::ios::binary);
			contactspace.data.clear();
			binaryReader(in2, contactspace);
			in2.close();
			
			std::ofstream out4("space_test6.txt");
			asciiWriter(out4, contactspace);
			out4.close();
			
			
		}
		
	}

}

void main()
{
	APDL::test_ContactSpace();
}