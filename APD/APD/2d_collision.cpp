#include "2d_collision.h"

namespace APDL
{
	Vec2D Simplex::getSearchDirection() const
	{
		switch(count)
		{
		case 1:
			return -verts[0].point;
		case 2:
		{
			Vec2D edge_AB = verts[1].point - verts[0].point;
			double sgn = cross(edge_AB, -verts[0].point);
			if(sgn > 0)
			{
				return cross(1.0, edge_AB);
			}
			else
			{
				return cross(edge_AB, 1.0);
			}
		}
		default:
			assert(0);
			return Vec2D(0, 0);
		}
	}
	
	Vec2D Simplex::getClosestPoint() const
	{
		switch(count)
		{
		case 1:
			return verts[0].point;
			
		case 2:
		{
			double s = 1.0 / divisor;
			return (s * verts[0].u) * verts[0].point + (s * verts[1].u) * verts[1].point;
		}
		
		case 3:
			return Vec2D(0, 0);
			
		default:
			assert(0);
			return Vec2D(0, 0);
		}
	}
	
	void Simplex::getWitnessPoints(Vec2D& point1, Vec2D& point2) const
	{
		double factor = 1.0 / divisor;
		
		switch(count)
		{
		case 1:
			point1 = verts[0].point1;
			point2 = verts[0].point2;
			break;
			
		case 2:
		{
			double s = 1.0 / divisor;
			point1 = (s * verts[0].u) * verts[0].point1 + (s * verts[1].u) * verts[1].point1;
			point2 = (s * verts[0].u) * verts[0].point2 + (s * verts[1].u) * verts[1].point2;
		}
		break;
		
		case 3:
		{
			double s = 1.0 / divisor;
			point1 = (s * verts[0].u) * verts[0].point1 + (s * verts[1].u) * verts[1].point1 + (s * verts[2].u) * verts[2].point1;
			point2 = point1;
		}
		break;
		
		default:
			assert(0);
			break;
		}
	}
	
	
	// Closest point on line segment to Q.
	// Voronoi regions: A, B, AB
	void Simplex::solve2(const Vec2D& Q)
	{
		Vec2D A = verts[0].point;
		Vec2D B = verts[1].point;
		
		// Compute barycentric coordinates (pre-division).
		double u = dot(Q - B, A - B);
		double v = dot(Q - A, B - A);
		
		// Region A
		if(v <= 0.0)
		{
			// Simplex is reduced to just vertex A.
			verts[0].u = 1.0;
			divisor = 1.0;
			count = 1;
			return;
		}
		
		// Region B
		if(u <= 0.0)
		{
			// Simplex is reduced to just vertex B.
			// We move vertex B into vertex A and reduce the count.
			verts[0] = verts[1];
			verts[0].u = 1.0;
			divisor = 1.0;
			count = 1;
			return;
		}
		
		// Region AB. Due to the conditions above, we are
		// guaranteed the the edge has non-zero length and division
		// is safe.
		verts[0].u = u;
		verts[1].u = v;
		Vec2D e = B - A;
		divisor = dot(e, e);
		count = 2;
	}
	
	// Closest point on triangle to Q.
	// Voronoi regions: A, B, C, AB, BC, CA, ABC
	void Simplex::solve3(const Vec2D& Q)
	{
		Vec2D A = verts[0].point;
		Vec2D B = verts[1].point;
		Vec2D C = verts[2].point;
		
		// Compute edge barycentric coordinates (pre-division).
		double uAB = dot(Q - B, A - B);
		double vAB = dot(Q - A, B - A);
		
		double uBC = dot(Q - C, B - C);
		double vBC = dot(Q - B, C - B);
		
		double uCA = dot(Q - A, C - A);
		double vCA = dot(Q - C, A - C);
		
		// Region A
		if(vAB <= 0.0 && uCA <= 0.0)
		{
			verts[0].u = 1.0;
			divisor = 1.0;
			count = 1;
			return;
		}
		
		// Region B
		if(uAB <= 0.0 && vBC <= 0.0)
		{
			verts[0] = verts[1];
			verts[0].u = 1.0;
			divisor = 1.0;
			count = 1;
			return;
		}
		
		// Region C
		if(uBC <= 0.0 && vCA <= 0.0)
		{
			verts[0] = verts[2];
			verts[0].u = 1.0;
			divisor = 1.0;
			count = 1;
			return;
		}
		
		// Compute signed triangle area.
		double area = cross(B - A, C - A);
		
		// Compute triangle barycentric coordinates (pre-division).
		double uABC = cross(B - Q, C - Q);
		double vABC = cross(C - Q, A - Q);
		double wABC = cross(A - Q, B - Q);
		
		// Region AB
		if(uAB > 0.0 && vAB > 0.0 && wABC * area <= 0.0)
		{
			verts[0].u = uAB;
			verts[1].u = vAB;
			Vec2D e = B - A;
			divisor = dot(e, e);
			count = 2;
			return;
		}
		
		// Region BC
		if(uBC > 0.0 && vBC > 0.0 && uABC * area <= 0.0)
		{
			verts[0] = verts[1];
			verts[1] = verts[2];
			
			verts[0].u = uBC;
			verts[1].u = vBC;
			Vec2D e = C - B;
			divisor = dot(e, e);
			count = 2;
			return;
		}
		
		// Region CA
		if(uCA > 0.0 && vCA > 0.0 && vABC * area <= 0.0)
		{
			verts[1] = verts[0];
			verts[0] = verts[2];
			
			verts[0].u = uCA;
			verts[1].u = vCA;
			Vec2D e = A - C;
			divisor = dot(e, e);
			count = 2;
			return;
		}
		
		// Region ABC
		// The triangle area is guaranteed to be non-zero.
		assert(uABC > 0.0 && vABC > 0.0 && wABC > 0.0);
		verts[0].u = uABC;
		verts[1].u = vABC;
		verts[2].u = wABC;
		divisor = area;
		count = 3;
	}
	
	
}