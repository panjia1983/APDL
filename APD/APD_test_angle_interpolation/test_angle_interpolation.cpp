#include <APD/math_utility.h>
#include <iostream>

namespace APDL
{
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
}

void main()
{
	APDL::test_angle_interpolation();
}