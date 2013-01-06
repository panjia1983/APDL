#ifndef PQP_EXTENSION_H
#define PQP_EXTENSION_H

#include <C2A/C2A.h>
#include <C2A/LinearMath.h>
#include <PQP/MatVec.h>
#include <limits>

namespace APDL
{

inline PQP_REAL max2(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

inline PQP_REAL min2(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}

inline int project6(PQP_REAL *ax, 
         PQP_REAL *p1, PQP_REAL *p2, PQP_REAL *p3, 
         PQP_REAL *q1, PQP_REAL *q2, PQP_REAL *q3)
{
  PQP_REAL P1 = VdotV(ax, p1);
  PQP_REAL P2 = VdotV(ax, p2);
  PQP_REAL P3 = VdotV(ax, p3);
  PQP_REAL Q1 = VdotV(ax, q1);
  PQP_REAL Q2 = VdotV(ax, q2);
  PQP_REAL Q3 = VdotV(ax, q3);
  
  PQP_REAL mx1 = max2(P1, P2, P3);
  PQP_REAL mn1 = min2(P1, P2, P3);
  PQP_REAL mx2 = max2(Q1, Q2, Q3);
  PQP_REAL mn2 = min2(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}

inline bool buildTrianglePlane(PQP_REAL v1[3], PQP_REAL v2[3], PQP_REAL v3[3], PQP_REAL n[3], PQP_REAL& t)
{
	PQP_REAL v21[3];
	PQP_REAL v31[3];
	v21[0] = v2[0] - v1[0]; v21[1] = v2[1] - v1[1]; v21[2] = v2[2] - v1[2];
	v31[0] = v3[0] - v1[0]; v31[1] = v3[1] - v1[1]; v31[2] = v3[2] - v1[2];
	VcrossV(n, v21, v31);
	
	PQP_REAL len = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	if(len > 0)
	{
		len = sqrt(len);
		n[0] /= len; n[1] /= len; n[2] /= len;
		t = n[0] * v1[0] + n[1] * v1[1] + n[2] * v1[2];
		return true;
	}

	return false;
}


inline PQP_REAL distanceToPlane(PQP_REAL n[3], PQP_REAL t, PQP_REAL v[3])
{
	return n[0] * v[0] + n[1] * v[1] + n[2] * v[2] - t;
}
inline void computeDeepestPoints(PQP_REAL* clipped_points, int num_clipped_points, PQP_REAL n[3], PQP_REAL t, 
								 PQP_REAL& penetration_depth, PQP_REAL* deepest_points, int& num_deepest_points)
{
  num_deepest_points = 0;
  PQP_REAL max_depth = -DBL_MAX;
  int num_deepest_points_ = 0;
  int num_neg = 0;
  int num_pos = 0;
  int num_zero = 0;

  for(int i = 0; i < num_clipped_points; ++i)
  {
    PQP_REAL dist = -distanceToPlane(n, t, &clipped_points[i * 3]);
    if(dist > 1e-6) num_pos++;
    else if(dist < -1e-6) num_neg++;
    else num_zero++;
    if(dist > max_depth)
    {
      max_depth = dist;
      num_deepest_points_ = 1;
	  int idx = num_deepest_points_ - 1;
      deepest_points[idx * 3] = clipped_points[i * 3];
	  deepest_points[idx * 3 + 1] = clipped_points[i * 3 + 1];
	  deepest_points[idx * 3 + 2] = clipped_points[i * 3 + 2];
    }
    else if(dist + 1e-6 >= max_depth)
    {
      num_deepest_points_++;
	  int idx = num_deepest_points - 1;
      deepest_points[idx * 3] = clipped_points[i * 3];
	  deepest_points[idx * 3 + 1] = clipped_points[i * 3 + 1];
	  deepest_points[idx * 3 + 2] = clipped_points[i * 3 + 2];
    }
  }

  if(max_depth < -1e-6)
    num_deepest_points_ = 0;

  if(num_zero == 0 && ((num_neg == 0) || (num_pos == 0)))
    num_deepest_points_ = 0;

  penetration_depth = max_depth;
  num_deepest_points = num_deepest_points_;
}


inline int TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
					  PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3,
					  PQP_REAL* contact_points,
					  int& num_contact_points,
					  PQP_REAL& penetration_depth,
					  PQP_REAL normal[3]) 
{
	// One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
	// Edges are (e1,e2,e3) and (f1,f2,f3).
	// Normals are n1 and m1
	// Outwards are (g1,g2,g3) and (h1,h2,h3).
	//  
	// We assume that the triangle vertices are in the same coordinate system.
	//
	// First thing we do is establish a new c.s. so that p1 is at (0,0,0).

	PQP_REAL p1[3], p2[3], p3[3];
	PQP_REAL q1[3], q2[3], q3[3];
	PQP_REAL e1[3], e2[3], e3[3];
	PQP_REAL f1[3], f2[3], f3[3];
	PQP_REAL g1[3], g2[3], g3[3];
	PQP_REAL h1[3], h2[3], h3[3];
	PQP_REAL n1[3], m1[3];

	PQP_REAL ef11[3], ef12[3], ef13[3];
	PQP_REAL ef21[3], ef22[3], ef23[3];
	PQP_REAL ef31[3], ef32[3], ef33[3];

	p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
	p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
	p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];

	q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
	q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
	q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];

	e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
	e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
	e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

	f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
	f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
	f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];

	VcrossV(n1, e1, e2);
	VcrossV(m1, f1, f2);

	VcrossV(g1, e1, n1);
	VcrossV(g2, e2, n1);
	VcrossV(g3, e3, n1);
	VcrossV(h1, f1, m1);
	VcrossV(h2, f2, m1);
	VcrossV(h3, f3, m1);

	VcrossV(ef11, e1, f1);
	VcrossV(ef12, e1, f2);
	VcrossV(ef13, e1, f3);
	VcrossV(ef21, e2, f1);
	VcrossV(ef22, e2, f2);
	VcrossV(ef23, e2, f3);
	VcrossV(ef31, e3, f1);
	VcrossV(ef32, e3, f2);
	VcrossV(ef33, e3, f3);

	// now begin the series of tests

	if (!project6(n1, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(m1, p1, p2, p3, q1, q2, q3)) return 0;

	if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return 0;

	if (!project6(g1, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(g2, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(g3, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(h1, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(h2, p1, p2, p3, q1, q2, q3)) return 0;
	if (!project6(h3, p1, p2, p3, q1, q2, q3)) return 0;

	PQP_REAL normal1[3];
	PQP_REAL normal2[3];
	PQP_REAL t1, t2;

	buildTrianglePlane(P1, P2, P3, normal1, t1);
	buildTrianglePlane(Q1, Q2, Q3, normal2, t2);

	PQP_REAL deepest_points1[9];
	int num_deepest_points1 = 0;
	PQP_REAL deepest_points2[9];
	int num_deepest_points2 = 0;

	PQP_REAL penetration_depth1, penetration_depth2;

	PQP_REAL P[9];
	P[0] = P1[0]; P[1] = P1[1]; P[2] = P1[2];
	P[3] = P2[0]; P[4] = P2[1]; P[5] = P2[2];
	P[6] = P3[0]; P[7] = P3[1]; P[8] = P3[2];

	PQP_REAL Q[9];
	Q[0] = Q1[0]; Q[1] = Q1[1]; Q[2] = Q1[2];
	Q[3] = Q2[0]; Q[4] = Q2[1]; Q[5] = Q2[2];
	Q[6] = Q3[0]; Q[7] = Q3[1]; Q[8] = Q3[2];


	computeDeepestPoints(Q, 3, normal1, t1, penetration_depth2, deepest_points2, num_deepest_points2);
	computeDeepestPoints(P, 3, normal2, t2, penetration_depth1, deepest_points1, num_deepest_points1);


	if(penetration_depth1 > penetration_depth2)
	{
		num_contact_points = (std::min)(num_deepest_points2, 2);
		for(int i = 0; i < num_contact_points; ++i)
		{
			contact_points[3 * i] = deepest_points2[3 * i];
			contact_points[3 * i + 1] = deepest_points2[3 * i + 1];
			contact_points[3 * i + 2] = deepest_points2[3 * i + 2];
		}

		normal[0] = normal1[0]; normal[1] = normal1[1]; normal[2] = normal1[2]; 
		penetration_depth = penetration_depth2;
	}
	else
	{
		num_contact_points = (std::min)(num_deepest_points1, 2);
		for(int i = 0; i < num_contact_points; ++i)
		{
			contact_points[3 * i] = deepest_points1[3 * i];
			contact_points[3 * i + 1] = deepest_points1[3 * i + 1];
			contact_points[3 * i + 2] = deepest_points1[3 * i + 2];
		}

		normal[0] = -normal2[0]; normal[1] = -normal2[1]; normal[2] = -normal2[2]; 
		penetration_depth = penetration_depth1;
	}
	
	return 1;
}

}

#endif