#include <APD/contact_space.h>
#include <APD/minkowski_cspace.h>
#include <APD/mesh_io.h>

#include <APD/profile.h>

namespace APDL
{
	static void test_collide_box()
	{
		Polygon p1;
		p1.points.push_back(Vec2D(1, -1));
		p1.points.push_back(Vec2D(1, 1));
		p1.points.push_back(Vec2D(-1, 1));
		p1.points.push_back(Vec2D(-1, -1));
		Polygon p2 = p1;
		Transform2D tf1;
		Transform2D tf2;
		tf2.t = Vec2D(1, 0);

		EPAResult res = doGJKEPA(p1, tf1, p2, tf2);
		for(std::size_t i = 0; i < res.contacts.size(); ++i)
		{
			std::cout << res.contacts[i].penetration << " " << res.contacts[i].normal.x << " " << res.contacts[i].normal.y 
				<< " " << res.contacts[i].point.x << " " << res.contacts[i].point.y << std::endl;
		}
		
	}

	static void test_collide_2d()
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

		ContactSpaceR2 contactspace(p1, p2, 2);

		std::vector<DataVector> samples = contactspace.uniform_sample0(1000);

		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			const DataVector& v = samples[i];
			Collider2D::CollisionResult res = contactspace.collider.collide(v);
			
			if(res.is_collide != contactspace.collider.isCollide(v)) std::cout << "error" << std::endl;
			//std::cout << res.is_collide << " " << res.contacts.size() << std::endl;
			//if(res.contacts.size() > 0)
			//{
			//	for(std::size_t j = 0; j < res.contacts.size(); ++j)
			//	{
			//		std::cout << res.contacts[j].penetration_depth << " ";
			//		std::cout << res.contacts[j].contact_point.x << " " << res.contacts[j].contact_point.y << " ";
			//		std::cout << res.contacts[j].normal.x << " " << res.contacts[j].normal.y << std::endl;
			//	}
			//}
		}
	}

	static void test_collide_2d_rotation()
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

		ContactSpaceSE2 contactspace(p1, p2, 2);

		std::vector<DataVector> samples = contactspace.uniform_sample0(1000);

		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			const DataVector& v = samples[i];

			Collider2D::CollisionResult res = contactspace.collider.collide(v);
			if(res.is_collide != contactspace.collider.isCollide(v)) std::cout << "error" << std::endl;
			//std::cout << res.is_collide << " " << res.contacts.size() << std::endl;
			//if(res.contacts.size() > 0)
			//{
			//	for(std::size_t j = 0; j < res.contacts.size(); ++j)
			//	{
			//		std::cout << res.contacts[j].penetration_depth << " ";
			//		std::cout << res.contacts[j].contact_point.x << " " << res.contacts[j].contact_point.y << " ";
			//		std::cout << res.contacts[j].normal.x << " " << res.contacts[j].normal.y << std::endl;
			//	}
			//}
		}
	}

	static void test_continuous_collide_2d()
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

		ContactSpaceR2 contactspace(p1, p2, 2);

		DataVector qt(3);
		while(1)
		{
			qt = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qt);
			if(res)
				break;
		}


		std::size_t id = 0;
		while(id < 1000)
		{
			DataVector qs = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qs);
			if(res) continue;

			id++;

			tools::Profiler::Begin("dcd");
			std::pair<bool, double> dcd_result = contactspace.collider.isCCDCollide(qs, qt, 1000);
			tools::Profiler::End("dcd");
			tools::Profiler::Begin("ccd");
			std::pair<bool, double> ccd_result = contactspace.collider.isCCDCollide(qs, qt);
			tools::Profiler::End("ccd");

			if(!ccd_result.first || ccd_result.second > dcd_result.second)
			{
				std::cout << "(" << dcd_result.first << "," << dcd_result.second << ")(" << ccd_result.first << "," << ccd_result.second << ")" << std::endl;
			}
		}
	}

	static void test_continuous_collide_2d_rotation()
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

		ContactSpaceSE2 contactspace(p1, p2, 2);

		DataVector qt(3);
		while(1)
		{
			qt = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qt);
			if(res)
				break;
		}


		std::size_t id = 0;
		while(id < 1000)
		{
			DataVector qs = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qs);
			if(res) continue;

			id++;

			tools::Profiler::Begin("dcd");
			std::pair<bool, double> dcd_result = contactspace.collider.isCCDCollide(qs, qt, 1000);
			tools::Profiler::End("dcd");
			tools::Profiler::Begin("ccd");
			std::pair<bool, double> ccd_result = contactspace.collider.isCCDCollide(qs, qt);
			tools::Profiler::End("ccd");

			if(!ccd_result.first || ccd_result.second > dcd_result.second)
			{
				std::cout << "(" << dcd_result.first << "," << dcd_result.second << ")(" << ccd_result.first << "," << ccd_result.second << ")" << std::endl;
			}
		}
	}


	static void test_collide_3d()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		ContactSpaceR3 contactspace(P, Q);

		std::vector<DataVector> samples = contactspace.uniform_sample0(1000);

		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			const DataVector& v = samples[i];

			Collider3D::CollisionResult res = contactspace.collider.collide(v);
			if(res.is_collide != contactspace.collider.isCollide(v)) std::cout << "error" << std::endl;
			//std::cout << res.is_collide << " " << res.contacts.size() << std::endl;
			//if(res.contacts.size() > 0)
			//{
			//	for(std::size_t j = 0; j < res.contacts.size(); ++j)
			//	{
			//		std::cout << res.contacts[j].penetration_depth << " ";
			//		std::cout << res.contacts[j].contact_point[0] << " " << res.contacts[j].contact_point[1] << " " << res.contacts[j].contact_point[2] << " ";
			//		std::cout << res.contacts[j].normal[0] << " " << res.contacts[j].normal[1] << " " << res.contacts[j].normal[2] << std::endl;
			//	}
			//}
		}

		delete P;
		delete Q;
	}

	static void test_collide_3d_rotation()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		ContactSpaceSE3Quat2 contactspace(P, Q);

		std::cout << P->com[0] << " " << P->com[1] << " " << P->com[2] << std::endl;
		std::cout << Q->com[0] << " " << Q->com[1] << " " << Q->com[2] << std::endl;

		std::vector<DataVector> samples = contactspace.uniform_sample0(1000);

		int sum = 0;
		for(std::size_t i = 0; i < samples.size(); ++i)
		{
		const DataVector& v = samples[i];

			Collider3D::CollisionResult res = contactspace.collider.collide(v);
			if(res.is_collide != contactspace.collider.isCollide(v)) std::cout << "error" << std::endl;
			if(res.is_collide) sum++;
			//std::cout << res.is_collide << " " << res.contacts.size() << std::endl;
			//if(res.contacts.size() > 0)
			//{
			//	for(std::size_t j = 0; j < res.contacts.size(); ++j)
			//	{
			//		std::cout << res.contacts[j].penetration_depth << " ";
			//		std::cout << res.contacts[j].contact_point[0] << " " << res.contacts[j].contact_point[1] << " " << res.contacts[j].contact_point[2] << " ";
			//		std::cout << res.contacts[j].normal[0] << " " << res.contacts[j].normal[1] << " " << res.contacts[j].normal[2] << std::endl;
			//	}
			//}
		}

		std::cout << sum / (double)samples.size() << std::endl;

		delete P;
		delete Q;
	}

	static void test_convex_decomposition_PD()
	{
		std::vector<C2A_Model*> model1;
		std::vector<C2A_Model*> model2;

		readObjFiles(model1, "../data/models/Bullet/ringz_convex.obj");
		readObjFiles(model2, "../data/models/Bullet/ringz_convex.obj");

		C2A_Model* P;
		C2A_Model* Q;
		readObjFile(P, "../data/models/Bullet/ringz.obj");
		readObjFile(Q, "../data/models/Bullet/ringz.obj");

		ContactSpaceSE3Euler2 contactspace(P, Q, 0.1);

		std::vector<DataVector> samples = contactspace.uniform_sample0(1000);

		for(std::size_t i = 0; i < samples.size(); ++i)
		{
			double PD = Collider3D::PDt(model1, model2, samples[i]);
			bool col = contactspace.collider.isCollide(samples[i]);

			bool col2 = false;
			for(std::size_t j = 0; j < model1.size(); ++j)
			{
				for(std::size_t k = 0; k < model2.size(); ++k)
				{
					//std::cout << j << " " << k  << std::endl;
					//for(std::size_t n = 0; n < model1[j]->num_tris; ++n)
					//{
					//	std::cout << model1[j]->num_tris << " " << n << std::endl;
					//	C2A_Tri* tri = model1[j]->GetTriangle(n);
					//	std::cout << tri->id << " " << tri->p1[0] << " " << tri->p1[1] << " " << tri->p1[2] << " ";
					//	std::cout << tri->p2[0] << " " << tri->p2[1] << " " << tri->p2[2] << " ";
					//	std::cout << tri->p3[0] << " " << tri->p3[1] << " " << tri->p3[2] << std::endl;
					//}

					//for(std::size_t n = 0; n < model2[k]->num_tris; ++n)
					//{
					//	std::cout << model1[j]->num_tris << " " << n << std::endl;
					//	C2A_Tri* tri = model2[k]->GetTriangle(n);
					//	std::cout << tri->id << " " << tri->p1[0] << " " << tri->p1[1] << " " << tri->p1[2] << " ";
					//	std::cout << tri->p2[0] << " " << tri->p2[1] << " " << tri->p2[2] << " ";
					//	std::cout << tri->p3[0] << " " << tri->p3[1] << " " << tri->p3[2] << std::endl;
					//}

					Collider3D collider(model1[j], model2[k]);
					if(collider.isCollide(samples[i]))
						col2 = true;
				}
			}
			std::cout << PD << " " << col2 << " " << col << std::endl;
		}
	}

	static void test_continuous_collide_3d()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		ContactSpaceR3 contactspace(P, Q);

	

		DataVector qt(6);
		while(1)
		{
			qt = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qt);
			if(res)
				break;
		}

		std::size_t id = 0;
		while(id < 1000)
		{
			DataVector qs = contactspace.uniform_sample0(1)[0];

			bool res = contactspace.collider.isCollide(qs);
			if(res) continue;

			id++;

			tools::Profiler::Begin("dcd");
			std::pair<bool, double> dcd_result = contactspace.collider.isCCDCollide(qs, qt, 1000);
			tools::Profiler::End("dcd");
			tools::Profiler::Begin("ccd");
			std::pair<bool, double> ccd_result = contactspace.collider.isCCDCollide(qs, qt);
			tools::Profiler::End("ccd");

			//if(!ccd_result.first || ccd_result.second > dcd_result.second)
			{
				std::cout << "(" << dcd_result.first << "," << dcd_result.second << ")(" << ccd_result.first << "," << ccd_result.second << ")" << std::endl;
			}
		}

		delete P;
		delete Q;
	}
}

void main()
{
	APDL::tools::Profiler::Start();
	// APDL::test_collide_box();
	// APDL::test_collide_2d();
	// APDL::test_collide_2d_rotation();
	// APDL::test_continuous_collide_2d();
	// APDL::test_continuous_collide_2d_rotation();
	
	// APDL::test_collide_3d();
	// APDL::test_collide_3d_rotation();
	// APDL::test_continuous_collide_3d();

	APDL::test_convex_decomposition_PD();

	APDL::tools::Profiler::Stop();
	APDL::tools::Profiler::Status();
	
}