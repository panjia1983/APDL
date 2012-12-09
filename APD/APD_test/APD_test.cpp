#include <APD/metric_learning.h>
#include <APD/minkowski_cspace.h>

#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_distance.h>
#include <APD/decision_boundary_sampler.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <APD/mesh_io.h>
#include <iostream>
#include <fstream>


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
			std::ifstream room_file("rooms_star.dat");
			
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
			std::ifstream cube_file("cube.off");
			
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
			
			std::pair<DataVector, double> PD_result2 = Minkowski_Cspace_3D::Exact_PD_SE3(query2, CSpace_SE3, 1);
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
			
			std::pair<DataVector, double> PD_result3 = Minkowski_Cspace_3D::Exact_PD_SE3(query3, CSpace_SE3, 1);
			DataVector PD_point3 = PD_result3.first;
			std::cout << "(" << PD_point3[0] << " " << PD_point3[1] << " " << PD_point3[2] << " "
			          << PD_point3[3] << " " << PD_point3[4] << " " << PD_point3[5] << " " << PD_point3[6] << ") " << PD_result3.second << std::endl;
		}
		
		
		//{
		//	Minkowski_Cspace_3D::Polyhedron P, Q;
		//	std::ifstream cup_file("cup.off");
		
		//	if(!cup_file.is_open())
		//	{
		//		std::cerr << "Failed to open the cup file." << std::endl;
		//		return;
		//	}
		
		//	cup_file >> P;
		//	cup_file.close();
		
		//	std::ifstream spoon_file("spoon.off");
		
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
		
		
		{
			Minkowski_Cspace_3D::Polyhedron P, Q;
			
			Minkowski_Cspace_3D::Point_3 p(1, 0, 0);
			Minkowski_Cspace_3D::Point_3 q(0, 1, 0);
			Minkowski_Cspace_3D::Point_3 r(0, 0, 1);
			Minkowski_Cspace_3D::Point_3 s(0, 0, 0);
			P.make_tetrahedron(p, q, r, s);
			
			std::ifstream spoon_file("spoon.off");
			
			if(!spoon_file.is_open())
			{
				std::cerr << "Failed to open the spoon file." << std::endl;
				return;
			}
			
			spoon_file >> Q;
			spoon_file.close();
			
			
			Minkowski_Cspace_3D::Polyhedron CSpace_R3 = Minkowski_Cspace_3D::Minkowski_Cobstacle_R3(P, Q);
		}
	}
	
	
	void test_angle_interpolation()
	{
		const double pi = boost::math::constants::pi<double>();
		{
			DataVector a(3), b(3);
			a[0] = 0;
			a[1] = 0;
			a[2] = pi / 3;
			b[0] = 0;
			b[1] = 0;
			b[2] = pi / 6;
			DataVector c = InterpConfig2D(a, b, 0.5);
			std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
		}
		
		{
			DataVector a(3), b(3);
			a[0] = 0;
			a[1] = 0;
			a[2] = -2 * pi / 3;
			b[0] = 0;
			b[1] = 0;
			b[2] = 2 * pi / 3;
			DataVector c = InterpConfig2D(a, b, 0.5);
			std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
		}
		
		{
			DataVector a(3), b(3);
			a[0] = 0;
			a[1] = 0;
			a[2] = 2 * pi / 3;
			b[0] = 0;
			b[1] = 0;
			b[2] = -2 * pi / 3;
			DataVector c = InterpConfig2D(a, b, 0.5);
			std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
		}
	}
	
	void test_GJK()
	{
		const double pi = boost::math::constants::pi<double>();
		{
			Mat2D M(1, 2, 3, 4);
			Vec2D v(1, 2);
		}
		{
			Polygon p1, p2;
			p1.points.resize(4);
			p1.points[0] = Vec2D(1, 1);
			p1.points[1] = Vec2D(-1, 1);
			p1.points[2] = Vec2D(-1, -1);
			p1.points[3] = Vec2D(1, -1);
			
			p2.points.resize(4);
			p2.points[0] = Vec2D(1, 1);
			p2.points[1] = Vec2D(-1, 1);
			p2.points[2] = Vec2D(-1, -1);
			p2.points[3] = Vec2D(1, -1);
			
			Collider2D collider(&p1, &p2);
			
			DataVector q(3);
			q[0] = 5;
			q[1] = 0.5;
			q[2] = pi / 4;
			Collider2D::DistanceResult res1 = collider.distance(q);
			std::cout << res1.distance << std::endl;
			
			q[0] = 2.0;
			q[1] = 0.5;
			q[2] = pi / 4;
			Collider2D::DistanceResult res2 = collider.distance(q);
			std::cout << res2.distance << std::endl;
		}
	}
	
	void test_ContactSpace()
	{
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			std::ifstream room_file("rooms_star.dat");
			
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
			readOffFile(model1, "cup.off");
			readOffFile(model2, "spoon.off");
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
			readObjFile(model1, "cup.obj");
			readObjFile(model2, "spoon.obj");
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
			readTriFile(model1, "cup.tri");
			readTriFile(model2, "spoon.tri");
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
			readTriFile(model1, "cup.tri");
			readTriFile(model2, "spoon.tri");
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
			readTriFile(model1, "cup.tri");
			readTriFile(model2, "spoon.tri");
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

	void test_conlitron_learner_2d()
	{
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceR2 contactspace(p1, p2, 2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.learn(contactspace.data, 2);

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			}


			int error_num = 0;
			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);
		}
	}

	void test_conlitron_learner_2d_rotation()
	{
		// scaled
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceSE2 contactspace(p1, p2, 2);
			for(int i = 0; i < 10000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace);

			DataVector w(3);
			w[0] = 1; w[1] = 1; w[2] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.learn(contactspace.data, 3);

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			}

			int error_num = 0;

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;

			learner.saveVisualizeData("conlitron_2d_rotation_vis.txt", contactspace.getScaler(), 100);
		}
	}
	
	void test_svm_learner_2d_rotation()
	{
		// scaled
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceSE2 contactspace(p1, p2, 2);
			for(int i = 0; i < 10000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			/// we need to change gamma!!! for scaled version
			learner.setGamma(0.5 * contactspace.getScaler().getScale() * 10);


			std::ofstream scaler_file("scaler_2d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model_2d_rotation.txt");

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			std::size_t error_num = 0;

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;


		}


		return;

		// no scaled
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceSE2 contactspace(p1, p2, 2);
			for(int i = 0; i < 10000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d_rotation.txt");
			asciiWriter(out, contactspace);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			//learner.setUseScaler(true);

			std::ofstream scaler_file("scaler_2d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model_2d_rotation.txt");

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			std::size_t error_num = 0;

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;
		}
	}


	void test_svm_learner()
	{

		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceR2 contactspace(p1, p2, 2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			// learner.setUseScaler(true);


			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			

			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			std::cout << std::endl;
		}

		return;

		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			SVMLearner learner;
			learner.setLinearClassifier();
			learner.setC(10);
			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model_2d.txt");

			HyperPlane hp = learner.getLinearModel();
			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				DataVector v(contactspace.active_data_dim());
				for(std::size_t j = 0; j < contactspace.active_data_dim(); ++j)
					v[j] = contactspace.data[i].v[j];
				double pred = hp.evaluate(v);
				if(pred > 0) 
					std::cout << 1 << " ";
				else 
					std::cout << 0 << " ";
			}
			std::cout << std::endl;

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
				std::cout << results[i].label << " ";
			std::cout << std::endl;
		}

	}

	void test_boosting_learner()
	{

		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			AdaBoostLearner learner;
			learner.learn(contactspace.data, contactspace.active_data_dim());

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			//for(std::size_t i = 0; i < results.size(); ++i)
			//{
			//	std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			//}
			//std::cout << std::endl;
		}
	}

	void test_incremental_svm()
	{
		test_incsvm(1.0, "data_kernel.txt");
		test_incsvm(1.0, "diabetes_kernel.txt");
	}

	void test_KNN()
	{
		{
			std::ifstream room_file("rooms_star.dat");
			
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

			flann::Index<ContactSpaceR2::DistanceType>* index = NULL;

			generateIndex<ContactSpaceR2::DistanceType>(contactspace.data, 
				contactspace.active_data_dim(), 
				index,
				flann::KDTreeIndexParams());

		    std::vector<std::vector<int> > indices;
		    std::vector<std::vector<double> > dists;

			knnSearch<ContactSpaceR2::DistanceType>(contactspace.data, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   10,
				   flann::SearchParams());

			indices.clear();
			dists.clear();

			radiusSearch<ContactSpaceR2::DistanceType>(contactspace.data, 
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
			std::ifstream room_file("rooms_star.dat");
			
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
			for(int i = 0; i < 100000; ++i)
				contactspace.random_sample();

			flann::Index<ContactSpaceR2::DistanceType>* index = NULL;

			generateIndex<ContactSpaceR2::DistanceType>(contactspace.data, 
				contactspace.active_data_dim(), 
				index,
				flann::KDTreeIndexParams());

		    std::vector<std::vector<int> > indices;
		    std::vector<std::vector<double> > dists;

			knnSearch<ContactSpaceR2::DistanceType>(contactspace.data, 
		           contactspace.active_data_dim(), 
		           index,
		           indices,
				   dists,
				   10,
				   flann::SearchParams());

			indices.clear();
			dists.clear();

			radiusSearch<ContactSpaceR2::DistanceType>(contactspace.data, 
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

	void test_distance_to_decision_boundary_conlitron()
	{
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceR2 contactspace(p1, p2, 2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.learn(contactspace.data, 2);

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			}


			int error_num = 0;
			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				if(results[i].label != contactspace.data[i].col) error_num++;
			}
			std::cout << "error ratio: " << error_num / (double)contactspace.data.size() << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);

			MulticonlitronDistanceToDecisionBoundary_BruteForce distancer1(learner);
			MulticonlitronDistanceToDecisionBoundary_KNN distancer2(learner);
			MulticonlitronDistanceToDecisionBoundary_EmbedKNN distancer3(learner);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				const DataVector& v = contactspace.data[i].v;
				double d1 = distancer1.distance(v);
				double d2 = distancer2.distance(v);
				double d3 = distancer3.distance(v);
				if(d1 > 5 || d2 > 5 || d3 > 5) 
					std::cout << d1 << " " << d2 << " " << d3 << std::endl;
			}

		}
	}


	void test_distance_to_decision_boundary_svm()
	{
		{
			std::ifstream room_file("rooms_star.dat");
			
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
			
			ContactSpaceR2 contactspace(p1, p2, 2);
			for(int i = 0; i < 1000; ++i)
				contactspace.random_sample();
				
			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());

			learner.setUseScaler(true);

			learner.learn(contactspace.data, contactspace.active_data_dim());
			learner.save("model.txt");

			std::vector<PredictResult> results = learner.predict(contactspace.data);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
				std::cout << "(" << results[i].label << "," << contactspace.data[i].col << ")";
			std::cout << std::endl;

			std::cout << learner.hyperw_normsqr << std::endl;


			SVMDistanceToDecisionBoundary_Bruteforce distancer1(learner);
			SVMDistanceToDecisionBoundary_OptimizationGradient distancer2(learner);
		    SVMDistanceToDecisionBoundary_Projection distancer3(learner);
			SVMDistanceToDecisionBoundary_Optimization distancer4(learner);
			SVMDistanceToDecisionBoundary_RoughLowerBound distancer5(learner);

			for(std::size_t i = 0; i < contactspace.data.size(); ++i)
			{
				DataVector v = (learner.scaler && learner.use_scaler) ? contactspace.getScaler().scale(contactspace.data[i].v) : contactspace.data[i].v;
				double d1 = distancer1.distance(v);
				double d2 = distancer2.distance(v);
				double d3 = distancer3.distance(v);
				double d4 = distancer4.distance(v);
				double d5 = distancer5.distance(v);
				std::cout << d1 << " " << d2 << " " << d3 << " " << d4 << " " << d5 << std::endl;
			}

			return;


			//std::vector<DataVector> samples;
			////sample_decision_boundary_interpolation(learner, samples);
			//sample_decision_boundary_hierarchial_tree<SVMDistanceToDecisionBoundary_OptimizationGradient>(learner, samples);
			//

			//std::cout << samples.size() << std::endl;

			//SVMDistanceToDecisionBoundary_Bruteforce distancer1(learner);
			//SVMDistanceToDecisionBoundary_OptimizationGradient distancer2(learner);
			//for(std::size_t i = 0; i < samples.size(); ++i)
			//	std::cout << "(" << distancer1.distance(samples[i]) << " " << distancer2.distance(samples[i]) << ") ";
			//std::cout << std::endl;
		}

	}
	
}

void main()
{
	// APDL::test_Minkowskiw_2D();
	
	// APDL::test_Minkowskiw_3D();
	
	// APDL::test_angle_interpolation();
	
	// APDL::test_GJK();
	
	// APDL::test_ContactSpace();
	
	APDL::test_svm_learner();

	// APDL::test_svm_learner_2d_rotation();

	// APDL::test_incremental_svm();

	// APDL::test_KNN();

	// APDL::test_boosting_learner();

	// APDL::test_distance_to_decision_boundary_svm();

	APDL::test_distance_to_decision_boundary_conlitron();

	// APDL::test_conlitron_learner_2d_rotation();

	// APDL::test_conlitron_learner_2d();

	// APDL::test_conlitron_learner_2d_rotation();
}

