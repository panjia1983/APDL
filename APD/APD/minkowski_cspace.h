#ifndef MINKOWSKI_CSPACE_H
#define MINKOWSKI_CSPACE_H

#include <CGAL/basic.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/aff_transformation_tags.h>
#include <iostream>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <limits>
#include "data_vector.h"
#include "math_utility.h"
#include "distance_proxy.h"
#include "2d_collision.h"

#include <C2A/C2A.h>
#include <C2A/LinearMath.h>

namespace APDL
{

	extern double distance_weight[7];
//-----------------------------------------------------------------------------
// Pretty-print a CGAL polygon.
//
	template<class Kernel, class Container>
	void print_polygon(const CGAL::Polygon_2<Kernel, Container>& P)
	{
		typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator  vit;
		
		std::cout << "[ " << P.size() << " vertices:";
		for(vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
			std::cout << " (" << *vit << ')';
		std::cout << " ]" << std::endl;
		
		return;
	}
	
//-----------------------------------------------------------------------------
// Pretty-print a polygon with holes.
//
	template<class Kernel, class Container>
	void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container>& pwh)
	{
		if(! pwh.is_unbounded())
		{
			std::cout << "{ Outer boundary = ";
			print_polygon(pwh.outer_boundary());
		}
		else
			std::cout << "{ Unbounded polygon." << std::endl;
			
		typename CGAL::Polygon_with_holes_2<Kernel, Container>::
		Hole_const_iterator  hit;
		unsigned int k = 1;
		
		std::cout << "  " << pwh.number_of_holes() << " holes:" << std::endl;
		for(hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k)
		{
			std::cout << "    Hole #" << k << " = ";
			print_polygon(*hit);
		}
		std::cout << " }" << std::endl;
		
		return;
	}
	
	
	
	class Minkowski_Cspace_2D
	{
	public:
		typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
		typedef Kernel::Point_2 Point_2;
		typedef Kernel::Segment_2 Segment_2;
		typedef CGAL::Polygon_2<Kernel> Polygon_2;
		typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;
		typedef CGAL::Aff_transformation_2<Kernel> Transformation;
		typedef CGAL::Cartesian<double> InExactKernel;
		typedef CGAL::Cartesian_converter<Kernel, InExactKernel> ToInexact;

		typedef Polygon_2::Container Container;
		typedef CGAL::Arr_segment_traits_2<Kernel>             Arr_segment_traits;
		typedef CGAL::Gps_segment_traits_2<Kernel,Container,Arr_segment_traits>  Traits_2;
		typedef CGAL::General_polygon_set_2<Traits_2>          General_polygon_set_2;
		
		static Polygon_2 rotate(const Polygon_2& in, double angle)
		{
			Transformation R(CGAL::ROTATION, sin(angle), cos(angle));
			return CGAL::transform(R, in);
		}
		
		static Polygon_with_holes_2 Minkowski_Cobstacle_R2(const Polygon_2& P, const Polygon_2& Q)
		{
			Transformation reflect(CGAL::SCALING, -1);
			return CGAL::minkowski_sum_2(P, CGAL::transform(reflect, Q));
		}

		static Polygon_with_holes_2 Minkowski_Cobstacle_R2(const std::vector<Polygon>& P_, const std::vector<Polygon>& Q_)
		{
			std::vector<Polygon_2> P, Q;

			for(std::size_t i = 0; i < P_.size(); ++i)
			{
				std::size_t num_p = P_[i].points.size(); 
				Point_2* p = new Point_2[num_p];
				for(std::size_t j = 0; j < num_p; ++j)
				{
					p[j] = Point_2(P_[i].points[j].x, P_[i].points[j].y);
				}

				Polygon_2 poly(p, p + num_p);
				delete [] p;

				P.push_back(poly);
			}

			for(std::size_t i = 0; i < Q_.size(); ++i)
			{
				std::size_t num_p = Q_[i].points.size(); 
				Point_2* p = new Point_2[num_p];
				for(std::size_t j = 0; j < num_p; ++j)
				{
					p[j] = Point_2(Q_[i].points[j].x, Q_[i].points[j].y);
				}

				Polygon_2 poly(p, p + num_p);
				delete [] p;

				Q.push_back(poly);
			}

			return Minkowski_Cobstacle_R2(P, Q);
		}

		static Polygon_with_holes_2 Minkowski_Cobstacle_R2(const std::vector<Polygon_2>& P, const std::vector<Polygon_2>& Q)
		{
			std::list<Polygon_2> sub_sum_polygons;

			for(std::size_t i = 0; i < P.size(); ++i)
			{
				for(std::size_t j = 0; j < Q.size(); ++j)
				{
					Polygon_with_holes_2 sub_sum = Minkowski_Cobstacle_R2(P[i], Q[j]);
					sub_sum_polygons.push_back(sub_sum.outer_boundary());
				}
			}

			General_polygon_set_2 gps;

			gps.join(sub_sum_polygons.begin(), sub_sum_polygons.end());

			std::list<Polygon_with_holes_2> sum;
			gps.polygons_with_holes(std::back_inserter(sum));

			return (*(sum.begin()));
		}
		
		static std::vector<std::pair<Polygon_with_holes_2, double> > Minkowski_CObstacle_SE2(const Polygon_2& P, const Polygon_2& Q, std::size_t n)
		{
			std::vector<std::pair<Polygon_with_holes_2, double> > ms;
			for(std::size_t i = 0; i < n; ++i)
			{
				if(i == 0)
					ms.push_back(std::make_pair(Minkowski_Cobstacle_R2(P, Q), 0));
				else
				{
					double angle = 2 * boost::math::constants::pi<double>() * i / (double)n;
					ms.push_back(std::make_pair(Minkowski_Cobstacle_R2(P, rotate(Q, angle)), angle));
				}
			}
			
			return ms;
		}

		static std::vector<std::pair<Polygon_with_holes_2, double> > Minkowski_CObstacle_SE2(const std::vector<Polygon_2>& P, const std::vector<Polygon_2>& Q, std::size_t n)
		{
			std::vector<std::pair<Polygon_with_holes_2, double> > ms;
			for(std::size_t i = 0; i < n; ++i)
			{
				if(i == 0)
					ms.push_back(std::make_pair(Minkowski_Cobstacle_R2(P, Q), 0));
				else
				{
					double angle = 2 * boost::math::constants::pi<double>() * i / (double)n;
					std::vector<Polygon_2> rotated_Q;
					for(std::size_t j = 0; j < Q.size(); ++j)
						rotated_Q.push_back(rotate(Q[j], angle));
					ms.push_back(std::make_pair(Minkowski_Cobstacle_R2(P, rotated_Q), angle));
				}
			}
			
			return ms;
		}

		static std::vector<std::pair<Polygon_with_holes_2, double> > Minkowski_CObstacle_SE2(const std::vector<Polygon>& P_, const std::vector<Polygon>& Q_, std::size_t n)
		{
			std::vector<Polygon_2> P, Q; 

			for(std::size_t i = 0; i < P_.size(); ++i)
			{
				std::size_t num_p = P_[i].points.size(); 
				Point_2* p = new Point_2[num_p];
				for(std::size_t j = 0; j < num_p; ++j)
				{
					p[j] = Point_2(P_[i].points[j].x, P_[i].points[j].y);
				}

				Polygon_2 poly(p, p + num_p);
				delete [] p;

				P.push_back(poly);
			}

			for(std::size_t i = 0; i < Q_.size(); ++i)
			{
				std::size_t num_p = Q_[i].points.size(); 
				Point_2* p = new Point_2[num_p];
				for(std::size_t j = 0; j < num_p; ++j)
				{
					p[j] = Point_2(Q_[i].points[j].x, Q_[i].points[j].y);
				}

				Polygon_2 poly(p, p + num_p);
				delete [] p;

				Q.push_back(poly);
			}

			return Minkowski_CObstacle_SE2(P, Q, n);
		}
		
		static bool inside(const Point_2& p, const Polygon_with_holes_2& CSpace)
		{
			for(Polygon_with_holes_2::Hole_const_iterator h_it = CSpace.holes_begin(); h_it != CSpace.holes_end(); ++h_it)
			{
				if(h_it->oriented_side(p) != CGAL::ON_POSITIVE_SIDE)
					return false;
			}
			
			return true;
		}
		
		static std::pair<DataVector, double> Exact_PD_R2(const DataVector& q, const Polygon_with_holes_2& CSpace)
		{
			ToInexact to_inexact;
			Point_2 query(q[0], q[1]);
			if(!inside(query, CSpace)) return std::make_pair(q, -1);
			
			Point_2 best_p;
			double best_d = (std::numeric_limits<double>::max)();
			
			std::vector<Polygon_2> polygons;
			if(!CSpace.is_unbounded())
			{
				polygons.push_back(CSpace.outer_boundary());
			}
			
			for(Polygon_with_holes_2::Hole_const_iterator h_it = CSpace.holes_begin(); h_it != CSpace.holes_end(); ++h_it)
				polygons.push_back(*h_it);
				
			for(std::size_t i = 0; i < polygons.size(); ++i)
			{
				for(Polygon_2::Vertex_const_iterator v_it = polygons[i].vertices_begin(); v_it != polygons[i].vertices_end(); ++v_it)
				{
					double d = to_inexact((*v_it - query).squared_length());
					if(d < best_d)
					{
						best_d = d;
						best_p = *v_it;
						
						if(best_d == 0) return std::make_pair(q, 0);
					}
				}
				
				for(Polygon_2::Edge_const_iterator e_it = polygons[i].edges_begin(); e_it != polygons[i].edges_end(); ++e_it)
				{
					const Segment_2& seg = *e_it;
					
					double l2 = to_inexact(seg.squared_length());
					if(l2 == 0) return std::make_pair(q, 0);
					
					double t = to_inexact((query - seg[0]) * (seg[1] - seg[0])) / l2;
					if(t > 0 && t < 1)
					{
						Point_2 proj = seg[0] + t * (seg[1] - seg[0]);
						double d = to_inexact((proj - query).squared_length());
						if(d < best_d)
						{
							best_d = d;
							best_p = proj;
						}
					}
				}
			}
			
			DataVector PD_p(2);
			PD_p[0] = to_inexact(best_p[0]);
			PD_p[1] = to_inexact(best_p[1]);
			return std::make_pair(PD_p, sqrt(best_d));
		}
		
		static std::pair<DataVector, double> Exact_PD_SE2(const DataVector& q, const std::vector<std::pair<Polygon_with_holes_2, double> >& CSpace)
		{
			DistanceProxySE2 metric;
			DataVector q2(2);
			q2[0] = q[0];
			q2[1] = q[1];
			DataVector best_p(3);
			double best_d = (std::numeric_limits<double>::max)();
			for(std::size_t i = 0; i < CSpace.size(); ++i)
			{
				std::pair<DataVector, double> res = Exact_PD_R2(q2, CSpace[i].first);
				double angle = CSpace[i].second;
				DataVector q_cur(3);
				q_cur[0] = res.first[0];
				q_cur[1] = res.first[1];
				q_cur[2] = angle;

				double d = metric.sqrDistance(q, q_cur);
				
				if((i == 0) || (d < best_d))
				{
					best_d = d;
					best_p[0] = q_cur[0];
					best_p[1] = q_cur[1];
					best_p[2] = q_cur[2];
				}
			}
			
			return std::make_pair(best_p, sqrt(best_d));
		}
	};
	
}

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/minkowski_sum_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

namespace APDL
{

	class Minkowski_Cspace_3D
	{
	public:
		typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
		typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
		typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
		typedef Kernel::Point_3 Point_3;
		typedef CGAL::Aff_transformation_3<Kernel> Transformation;
		
		typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel, Polyhedron> Primitive;
		typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
		typedef CGAL::AABB_tree<Traits> Tree;
		
		typedef CGAL::Cartesian<double> InExactKernel;
		typedef CGAL::Cartesian_converter<Kernel, InExactKernel> ToInexact;
		
		
		
		static Nef_polyhedron rotate(const Nef_polyhedron& P, double M[3][3])
		{
			Transformation R(M[0][0], M[0][1], M[0][2],
			                 M[1][0], M[1][1], M[1][2],
			                 M[2][0], M[2][1], M[2][2]);
			                 
			Nef_polyhedron P_(P);
			P_.transform(R);
			return P_;
		}
		
		static Nef_polyhedron Minkowski_Cobstacle_R3(const Nef_polyhedron& P, const Nef_polyhedron& Q)
		{
			Transformation reflection(CGAL::SCALING, -1);
			
			// compute the reflection of Q, cannot directly apply reflection to Q because it will cause assertion violation problem in CGAL
			// so we first reflect the polyhedron in Q and then convert it to Nef_polyhedron.
			Polyhedron Q_;
			Q.convert_to_polyhedron(Q_);
			std::transform(Q_.points_begin(), Q_.points_end(), Q_.points_begin(), reflection);
			
			return CGAL::minkowski_sum_3(const_cast<Nef_polyhedron&>(P), Nef_polyhedron(Q_));
		}
		
		static Polyhedron Minkowski_Cobstacle_R3(const Polyhedron& P, const Polyhedron& Q)
		{
		
			Nef_polyhedron result_ = Minkowski_Cobstacle_R3(Nef_polyhedron(const_cast<Polyhedron&>(P)), Nef_polyhedron(const_cast<Polyhedron&>(Q)));
			
			if(result_.is_simple())
			{
				Polyhedron result;
				result_.convert_to_polyhedron(result);
				return result;
			}
			else
			{
				std::cout << "result is not simple" << std::endl;
				return Polyhedron();
			}
		}
		
		
		
		static std::vector<std::pair<Nef_polyhedron, Quaternion> > Minkowski_Cobstacle_SE3(const Nef_polyhedron& P, const Nef_polyhedron& Q, std::size_t n)
		{
			std::vector<std::pair<Nef_polyhedron, Quaternion> > ms;
			for(std::size_t i = 0; i < 2 * n; ++i)
			{
				double psi = boost::math::constants::pi<double>() / n * i;
				for(std::size_t j = 0; j < n; ++j)
				{
					double theta = boost::math::constants::pi<double>() / n * i;
					for(std::size_t k = 0; k < 2 * n; ++k)
					{
						std::cout << i << " " << j << " " << k << std::endl;
						double phi = boost::math::constants::pi<double>() / n * i;
						
						double c = cos(theta * 0.5);
						double s = sin(theta * 0.5);
						double x = c * cos(psi * 0.5);
						double y = c * sin(psi * 0.5);
						double z = s * cos(psi * 0.5 + phi);
						double w = s * sin(psi * 0.5 + phi);
						
						double R[3][3];
						Quat2Rot(R, Quaternion(x, y, z, w));
						
						ms.push_back(std::make_pair(Minkowski_Cobstacle_R3(P, rotate(Q, R)), Quaternion(x, y, z, w)));
					}
				}
			}
			
			return ms;
		}
		
		static std::vector<std::pair<Polyhedron, Quaternion> > Minkowski_Cobstacle_SE3(const Polyhedron& P, const Polyhedron& Q, std::size_t n)
		{
			std::vector<std::pair<Polyhedron, Quaternion> > results;
			
			std::vector<std::pair<Nef_polyhedron, Quaternion> > results_ = Minkowski_Cobstacle_SE3(Nef_polyhedron(const_cast<Polyhedron&>(P)), Nef_polyhedron(const_cast<Polyhedron&>(Q)), n);
			for(std::size_t i = 0; i < results_.size(); ++i)
			{
				Polyhedron result;
				results_[i].first.convert_to_polyhedron(result);
				results.push_back(std::make_pair(result, results_[i].second));
			}
			
			return results;
		}
		
		
		static std::pair<DataVector, double> Exact_PD_R3(const DataVector& q, const Polyhedron& CSpace)
		{
			Point_3 query(q[0], q[1], q[2]);
			
			bool is_inside = true;
			// check whether query is inside the CSpace
			for(Polyhedron::Plane_const_iterator p_it = CSpace.planes_begin(); p_it != CSpace.planes_end(); ++p_it)
			{
				//if(!Kernel::HasOnNegativeSide_3(*p_it, query))
				if(p_it->oriented_side(query) == CGAL::ON_POSITIVE_SIDE)
				{
					is_inside = false;
					break;
				}
			}
			
			if(!is_inside)
				return std::make_pair(q, -1);
				
				
			Tree tree(const_cast<Polyhedron&>(CSpace).facets_begin(), const_cast<Polyhedron&>(CSpace).facets_end());
			tree.accelerate_distance_queries();
			
			ToInexact to_inexact;
			
			double best_d = to_inexact(tree.squared_distance(query));
			Point_3 best_p = tree.closest_point(query);
			
			DataVector PD_p(3);
			PD_p[0] = to_inexact(best_p[0]);
			PD_p[1] = to_inexact(best_p[1]);
			PD_p[2] = to_inexact(best_p[2]);
			return std::make_pair(PD_p, sqrt(best_d));
		}
		
		static std::pair<DataVector, double> Exact_PD_SE3(const DataVector& q, const std::vector<std::pair<Polyhedron, Quaternion> >& CSpace)
		{
			if(q.dim() == 6)
			{
				DistanceProxySE3EulerAngle metric;
				DataVector q_trans(3);
				q_trans[0] = q[0];
				q_trans[1] = q[1];
				q_trans[2] = q[2];
				
				DataVector best_p(6);
				double best_d = (std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < CSpace.size(); ++i)
				{
					std::pair<DataVector, double> res_trans = Exact_PD_R3(q_trans, CSpace[i].first);
					double a, b, c;
					Quat2Euler(a, b, c, CSpace[i].second);
					
					DataVector res(6);
					res[0] = res_trans.first[0];
					res[1] = res_trans.first[1];
					res[2] = res_trans.first[2];
					res[3] = a;
					res[4] = b;
					res[5] = c;
					
					double d = metric.sqrDistance(q, res);
					if((i == 0) || (d < best_d))
					{
						best_d = d;
						best_p = res;
					}
				}
				
				return std::make_pair(best_p, sqrt(best_d));
			}
			else if(q.dim() == 7)
			{
				DistanceProxySE3Quaternion metric;
				DataVector q_trans(3);
				q_trans[0] = q[0];
				q_trans[1] = q[1];
				q_trans[2] = q[2];
				
				DataVector best_p(7);
				double best_d = (std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < CSpace.size(); ++i)
				{
					std::pair<DataVector, double> res_trans = Exact_PD_R3(q_trans, CSpace[i].first);
					
					DataVector res(7);
					res[0] = res_trans.first[0];
					res[1] = res_trans.first[1];
					res[2] = res_trans.first[2];
					
					const Quaternion& r = CSpace[i].second;
					res[3] = r[0];
					res[4] = r[1];
					res[5] = r[2];
					res[6] = r[3];
					
					double d = metric.sqrDistance(q, res);
					
					if((i == 0) || (d < best_d))
					{
						best_d = d;
						best_p = res;
					}
				}
				
				return std::make_pair(best_p, sqrt(best_d));
			}
			else
				return std::make_pair(q, -1);
		}
		
		static std::pair<DataVector, double> Exact_PD_R3(const DataVector& q, C2A_Model* CSpace)
		{
			typedef CGAL::Simple_cartesian<double> K;

			typedef K::FT FT;
			typedef K::Ray_3 Ray;
			typedef K::Line_3 Line;
			typedef K::Point_3 Point;
			typedef K::Triangle_3 Triangle;

			typedef std::list<Triangle>::iterator Iterator;
			typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
			typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
			typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;			 
			
			std::list<Triangle> triangles;
			for(std::size_t i = 0; i < CSpace->num_tris; ++i)
			{
				C2A_Tri* tri = CSpace->GetTriangle(i);
				Point a(tri->p1[0], tri->p1[1], tri->p1[2]);
				Point b(tri->p2[0], tri->p2[1], tri->p2[2]);
				Point c(tri->p3[0], tri->p3[1], tri->p3[2]);
				triangles.push_back(Triangle(a, b, c));
			}

			Tree tree(triangles.begin(), triangles.end());
			Point query(q[0], q[1], q[2]);
			Point closest_point = tree.closest_point(query);
			
			DataVector PD_p(3);
			PD_p[0] = closest_point.x();
			PD_p[1] = closest_point.y();
			PD_p[2] = closest_point.z();

			double dist = tree.squared_distance(query);
			return std::make_pair(PD_p, sqrt(dist));
		}

		static std::pair<DataVector, double> Exact_PD_SE3(const DataVector& q, const std::vector<std::pair<C2A_Model*, Quaternion> >& CSpace)
		{
			if(q.dim() == 6)
			{
				DistanceProxySE3EulerAngle metric;
				DataVector q_trans(3);
				q_trans[0] = q[0];
				q_trans[1] = q[1];
				q_trans[2] = q[2];
				
				DataVector best_p(6);
				double best_d = (std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < CSpace.size(); ++i)
				{
					std::pair<DataVector, double> res_trans = Exact_PD_R3(q_trans, CSpace[i].first);
					double a, b, c;
					Quat2Euler(a, b, c, CSpace[i].second);
					
					DataVector res(6);
					res[0] = res_trans.first[0];
					res[1] = res_trans.first[1];
					res[2] = res_trans.first[2];
					res[3] = a;
					res[4] = b;
					res[5] = c;
					
					double d = metric.sqrDistance(q, res);

					if((i == 0) || (d < best_d))
					{
						best_d = d;
						best_p = res;
					}
				}
				
				return std::make_pair(best_p, sqrt(best_d));
			}
			else if(q.dim() == 7)
			{
				DistanceProxySE3Quaternion metric;
				DataVector q_trans(3);
				q_trans[0] = q[0];
				q_trans[1] = q[1];
				q_trans[2] = q[2];
				
				DataVector best_p(7);
				double best_d = (std::numeric_limits<double>::max)();
				for(std::size_t i = 0; i < CSpace.size(); ++i)
				{
					std::pair<DataVector, double> res_trans = Exact_PD_R3(q_trans, CSpace[i].first);
					
					DataVector res(7);
					res[0] = res_trans.first[0];
					res[1] = res_trans.first[1];
					res[2] = res_trans.first[2];
					
					const Quaternion& r = CSpace[i].second;
					res[3] = r[0];
					res[4] = r[1];
					res[5] = r[2];
					res[6] = r[3];
					
					double d = metric.sqrDistance(q, res);
					
					if((i == 0) || (d < best_d))
					{
						best_d = d;
						best_p = res;
					}
				}
				
				return std::make_pair(best_p, sqrt(best_d));
			}
			else
				return std::make_pair(q, -1);
		}
	};
	
	
	static std::ofstream& operator << (std::ofstream& os, const Minkowski_Cspace_3D::Polyhedron& P)
	{
		CGAL::set_ascii_mode(os);
		os << "OFF" << std::endl << P.size_of_vertices() << ' '
		   << P.size_of_facets() << " 0" << std::endl;
		std::copy(P.points_begin(), P.points_end(), std::ostream_iterator<Minkowski_Cspace_3D::Point_3>(os, "\n"));
		for(Minkowski_Cspace_3D::Polyhedron::Facet_const_iterator i = P.facets_begin(); i != P.facets_end(); ++i)
		{
			Minkowski_Cspace_3D::Polyhedron::Halfedge_around_facet_const_circulator j = i->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			os << CGAL::circulator_size(j) << ' ';
			do
			{
				os << ' ' << std::distance(P.vertices_begin(), j->vertex());
			}
			while(++j != i->facet_begin());
			os << std::endl;
		}
		
		return os;
	}
	
	
}


#endif