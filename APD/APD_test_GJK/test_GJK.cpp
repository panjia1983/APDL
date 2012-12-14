#include <APD/collider.h>

namespace APDL
{
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
}

void main()
{
	APDL::test_GJK();
}