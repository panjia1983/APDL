#include <APD/contact_space.h>
#include <APD/minkowski_cspace.h>

namespace APDL
{
	void test_collide()
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

		for(std::size_t i = 0; i < 1000; ++i)
		{
			DataVector v(3); 

			DataVector v_ = contactspace.sampler.sample();

			v[0] = v_[0];
			v[1] = v_[1];

			Collider2D::CollisionResult res = contactspace.collider.collide(v);
			std::cout << res.is_collide << " " << res.contacts.size() << std::endl;
		}

	}
}

void main()
{
	APDL::test_collide();
}