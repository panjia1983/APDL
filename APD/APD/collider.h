#ifndef COLLIDER_H
#define COLLIDER_H

#include <C2A/C2A.h>
#include <C2A/LinearMath.h>
#include "data_vector.h"
#include "math_utility.h"

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <cassert>
#include <list>
#include <limits>

#include <boost/shared_ptr.hpp>

#include "2d_collision.h"
#include "PQP_extension.h"

namespace APDL
{

	template<typename CGAL_Polygon, typename Kernel>
	Polygon toPolygon(const CGAL_Polygon& in)
	{
		typedef CGAL::Cartesian<double>                             InExactKernel;
		typedef CGAL::Cartesian_converter<Kernel, InExactKernel>    ToInexact;
		
		ToInexact to_inexact;
		Polygon out;
		for(CGAL_Polygon::Vertex_iterator it = in.vertices_begin(); it != in.vertices_end(); ++it)
		{
			double x = to_inexact(it->x());
			double y = to_inexact(it->y());
			out.points.push_back(Vec2D(x, y));
		}
		
		return out;
	}
	
	class Collider2D
	{
	public:
		typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
		typedef CGAL::Partition_traits_2<Kernel>                    Traits;
		typedef Traits::Point_2                                     Point_2;
		typedef Traits::Polygon_2                                   Polygon_2;
		typedef Polygon_2::Vertex_iterator                          Vertex_iterator;
		typedef std::list<Polygon_2>                                Polygon_list;
		typedef CGAL::Aff_transformation_2<Kernel> Transformation;
		typedef CGAL::Cartesian<double>                             InExactKernel;
		typedef CGAL::Cartesian_converter<Kernel, InExactKernel>    ToInexact;
		
		struct DistanceResult
		{
			Vec2D point1;
			Vec2D point2;
			double distance;
		};

		struct Contact
		{
			Vec2D normal;
			Vec2D contact_point;
			double penetration_depth;

			Contact(const Vec2D& point_, const Vec2D& normal_, double penetration_)
			{
				contact_point = point_;
				normal = normal_;
				penetration_depth = penetration_;
			}
		};

		struct CollisionResult
		{
			bool is_collide;

			std::vector<Contact> contacts;

			CollisionResult()
			{
				is_collide = false;
			}
			
		};
		
		Collider2D(const Polygon_2* model1_, const Polygon_2* model2_)
		{
			set(model1_, model2_);
		}
		
		Collider2D(const Polygon* model1_, const Polygon* model2_)
		{
			Polygon_2 m1, m2;
			for(std::size_t i = 0; i < model1_->size(); ++i)
			{
				m1.push_back(Point_2(model1_->points[i].x, model1_->points[i].y));
			}
			
			for(std::size_t i = 0; i < model2_->size(); ++i)
			{
				m2.push_back(Point_2(model2_->points[i].x, model2_->points[i].y));
			}
			
			set(&m1, &m2);
		}
		
		void set(const Polygon_2* model1_, const Polygon_2* model2_)
		{
			ToInexact to_inexact;
			model1 = *model1_;
			model2 = *model2_;
			
			if(!model1_->is_convex())
			{
				Polygon_list convex_list;
				CGAL::approx_convex_partition_2(model1_->vertices_begin(),
				                                model1_->vertices_end(),
				                                std::back_inserter(convex_list));
				                                
				for(Polygon_list::const_iterator it = convex_list.begin();
				        it != convex_list.end();
				        ++it)
				{
					const Polygon_2& cgal_poly = *it;
					Polygon poly;
					poly.points.resize(cgal_poly.size());
					
					int i = 0;
					for(Polygon_2::Vertex_const_iterator v_it = cgal_poly.vertices_begin();
					        v_it != cgal_poly.vertices_end();
					        ++v_it)
					{
						double x = to_inexact(v_it->x());
						double y = to_inexact(v_it->y());
						poly.points[i++] = Vec2D(x, y);
					}
					
					model1_convex.push_back(poly);
				}
			}
			else
			{
				Polygon poly;
				poly.points.resize(model1_->size());
				
				int i = 0;
				for(Polygon_2::Vertex_const_iterator v_it = model1_->vertices_begin();
				        v_it != model1_->vertices_end();
				        ++v_it)
				{
					double x = to_inexact(v_it->x());
					double y = to_inexact(v_it->y());
					poly.points[i++] = Vec2D(x, y);
				}
				
				model1_convex.push_back(poly);
			}
			
			
			if(!model2_->is_convex())
			{
				Polygon_list convex_list;
				CGAL::approx_convex_partition_2(model2_->vertices_begin(),
				                                model2_->vertices_end(),
				                                std::back_inserter(convex_list));
				                                
				for(Polygon_list::const_iterator it = convex_list.begin();
				        it != convex_list.end();
				        ++it)
				{
					const Polygon_2& cgal_poly = *it;
					Polygon poly;
					poly.points.resize(cgal_poly.size());
					
					int i = 0;
					for(Polygon_2::Vertex_const_iterator v_it = cgal_poly.vertices_begin();
					        v_it != cgal_poly.vertices_end();
					        ++v_it)
					{
						double x = to_inexact(v_it->x());
						double y = to_inexact(v_it->y());
						poly.points[i++] = Vec2D(x, y);
					}
					
					model2_convex.push_back(poly);
				}
			}
			else
			{
				Polygon poly;
				poly.points.resize(model2_->size());
				
				int i = 0;
				for(Polygon_2::Vertex_const_iterator v_it = model2_->vertices_begin();
				        v_it != model2_->vertices_end();
				        ++v_it)
				{
					double x = to_inexact(v_it->x());
					double y = to_inexact(v_it->y());
					poly.points[i++] = Vec2D(x, y);
				}
				
				model2_convex.push_back(poly);
			}

			max_distance_to_origin_model2 = 0;
			for(Polygon_2::Vertex_const_iterator v_it = model2_->vertices_begin();
				        v_it != model2_->vertices_end();
				        ++v_it)
			{
				double x = to_inexact(v_it->x());
				double y = to_inexact(v_it->y());
				double dist = x * x + y * y;
				if(dist > max_distance_to_origin_model2) max_distance_to_origin_model2 = dist;
			}
			max_distance_to_origin_model2 = sqrt(max_distance_to_origin_model2);
		}
		
		DistanceResult distance(const DataVector& q) const
		{
			GJKResult global_res;
			global_res.distance = (std::numeric_limits<double>::max)();
			double c = cos(q[2]);
			double s = sin(q[2]);
			
			Transform2D tf1;
			Transform2D tf2(Mat2D(c, -s, s, c), Vec2D(q[0], q[1]));
			
			for(std::size_t i = 0; i < model1_convex.size(); ++i)
			{
				for(std::size_t j = 0; j < model2_convex.size(); ++j)
				{
				
					GJKResult res = doGJK(model1_convex[i], tf1,
					                      model2_convex[j], tf2);
					                      
					if(res.distance < global_res.distance)
						global_res = res;
				}
			}
			
			DistanceResult dist_res;
			
			dist_res.point1 = global_res.point1;
			dist_res.point2 = global_res.point2;
			dist_res.distance = global_res.distance;
			
			return dist_res;
		}
		
		CollisionResult collide(const DataVector& q) const
		{
			double c = cos(q[2]);
			double s = sin(q[2]);
			Transform2D tf1;
			Transform2D tf2(Mat2D(c, -s, s, c), Vec2D(q[0], q[1]));

			CollisionResult res;
			res.is_collide = false;

			for(std::size_t i = 0; i < model1_convex.size(); ++i)
			{
				for(std::size_t j = 0; j < model2_convex.size(); ++j)
				{
					EPAResult epa_res = doGJKEPA(model1_convex[i], tf1,
											     model2_convex[j], tf2);

					if(epa_res.distance <= 0)
					{
						res.is_collide = true;

						for(std::size_t k = 0; k < epa_res.contacts.size(); ++k)
						{
							res.contacts.push_back(Contact(epa_res.contacts[k].point, epa_res.contacts[k].normal, epa_res.contacts[k].penetration));
						}
					}
				}
			}

			return res;
		}

		bool isCollideCGAL(const DataVector& q) const
		{
			double c = cos(q[2]);
			double s = sin(q[2]);
			Transformation R(c, -s, q[0], s, c, q[1]);
			Polygon_2 model2_new = CGAL::transform(R, model2);
			
			return CGAL::do_intersect(model1, model2_new);
		}
		
		bool isCollide(const DataVector& q) const
		{
			double c = cos(q[2]);
			double s = sin(q[2]);
			Transform2D tf1;
			Transform2D tf2(Mat2D(c, -s, s, c), Vec2D(q[0], q[1]));
			
			for(std::size_t i = 0; i < model1_convex.size(); ++i)
			{
				for(std::size_t j = 0; j < model2_convex.size(); ++j)
				{
				
					GJKResult res = doGJK(model1_convex[i], tf1,
					                      model2_convex[j], tf2);
					                      
					if(res.distance <= 0)
						return true;
				}
			}
			
			return false;
		}


		static DataVector interpolate(const DataVector& qs, const DataVector& qt, double t)
		{
			return InterpConfig2D(qs, qt, t);
		}
		
		std::pair<bool, double> isCCDCollide(const DataVector& qs, const DataVector& qt, std::size_t n_dcd) const
		{
			// perform n_dcd + 1 collisions
			for(std::size_t i = 0; i <= n_dcd; ++i)
			{
				if(i == 0)
				{
					if(isCollide(qs)) return std::make_pair(true, 0);
				}
				else if(i == n_dcd)
				{
					if(isCollide(qt)) return std::make_pair(true, 1);
				}
				else
				{
					double t = i / (double)n_dcd;
					
					DataVector q = InterpConfig2D(qs, qt, t);
					if(isCollide(q)) return std::make_pair(true, t);
				}
			}
			
			return std::make_pair(false, 1);
		}
		
		std::pair<bool, double> isCCDCollide(const DataVector& qs, const DataVector& qt) const
		{
			double toc = 0;
			double t = 0;
			double d;
			double u;
			Vec2D n;
			Vec2D p1, p2;
			
			double t_threshold = 1e-6;
			double d_threshold = 1e-6;
			
			Transform2D tf1;
			double c = cos(qs[2]), s = sin(qs[2]);
			Transform2D tf2(Mat2D(c, -s, s, c), Vec2D(qs[0], qs[1]));
			Vec2D v(qt[0] - qs[0], qt[1] - qs[1]);
			double w = angleTruncate(angleTruncate(qt[2]) - angleTruncate(qs[2]));

			double v_len = v.length();

			DataVector q(3);
			Vec2D cached_n(2, 1);
			
			do
			{
				q[0] = qs[0] + v.x * toc;
				q[1] = qs[1] + v.y * toc;
				q[2] = qs[2] + w * toc;
				
				DistanceResult res = distance(q);
				n = res.point2 - res.point1;
				d = n.length();
				if(d > 0)
				{
					n *= 1 / d;
					cached_n = n;
				}
				else 
					n = cached_n;

				//u = max(0, dot(v, n) + abs(q[2]) * max_distance_to_origin_model2); /// can not use this, because this is for convex shapes
				u = v_len + abs(q[2]) * max_distance_to_origin_model2;

				if(d > d_threshold && u > 0)
					t = d / u;
				else break;
				
				toc += t;
			}
			while(t > t_threshold && toc < 1);

			if(toc < 1)
				return std::make_pair(true, toc);
			else
				return std::make_pair(false, 1);
		}
		
	protected:
	
		std::vector<Polygon> model1_convex;
		std::vector<Polygon> model2_convex;
		
		Polygon_2 model1;
		Polygon_2 model2;
		double max_distance_to_origin_model2;  // for ccd

		
	};
	
	class Collider3D
	{
	public:
	
		struct DistanceResult
		{
			double point1[3];
			double point2[3];
			double distance;
		};

		struct Contact
		{
			double normal[3];
			double contact_point[3];
			double penetration_depth;

			Contact(double normal_[3], double contact_point_[3], double penetration_depth_)
			{
				for(std::size_t i = 0; i < 3; ++i)
				{
					normal[i] = normal_[i];
					contact_point[i] = contact_point_[i];
				}

				penetration_depth = penetration_depth_;
			}
		};

		struct CollisionResult
		{
			bool is_collide;

			std::vector<Contact> contacts;

			CollisionResult()
			{
				is_collide = false;
			}
		};
		
		Collider3D(C2A_Model* model1_, C2A_Model* model2_) : model1(model1_), model2(model2_)
		{}
		
		DistanceResult distance(const DataVector& q) const
		{
			double R1[3][3];
			double T1[3];
			
			IdentityRotTrans(R1, T1);
			
			double R2[3][3];
			double T2[3];
			
			if(q.dim() == 6)
				ConfigEuler2RotTran(R2, T2, q);
			else if(q.dim() == 7)
				ConfigQuat2RotTrans(R2, T2, q);
			else
				return DistanceResult();
				
			C2A_DistanceResult result;
			C2A_Distance(&result, R1, T1, model1.get(), R2, T2, model2.get(), 0.001, 0.001);
			
			DistanceResult res;
			res.distance = result.Distance();
			res.point1[0] = result.p1[0];
			res.point1[1] = result.p1[1];
			res.point1[2] = result.p1[2];
			
			for(int i = 0; i < 3; ++i)
			{
				res.point2[i] = 0;
				for(int j = 0; j < 3; ++j)
					res.point2[i] += R2[i][j] * result.p2[j];
					
				res.point2[i] += T2[i];
			}
			
			return res;
		}

		CollisionResult collide(const DataVector& q) const
		{
			CollisionResult res;
			res.is_collide = false;

			double R1[3][3];
			double T1[3];
			
			IdentityRotTrans(R1, T1);
			
			double R2[3][3];
			double T2[3];
			
			if(q.dim() == 6)
				ConfigEuler2RotTran(R2, T2, q);
			else if(q.dim() == 7)
				ConfigQuat2RotTrans(R2, T2, q);
			else
				return res;
				
			PQP_CollideResult pqp_res;
			
			C2A_Collide(&pqp_res, R1, T1, model1.get(), R2, T2, model2.get(), C2A_FIRST_CONTACT);

			std::cout << pqp_res.num_pairs << std::endl;
			if(pqp_res.NumPairs() == 0)
			{
				res.is_collide = false;
			}
			else
			{
				res.is_collide = true;

				for(std::size_t i = 0; i < pqp_res.NumPairs(); ++i)
				{
					int id1 = pqp_res.pairs[i].id1;
					int id2 = pqp_res.pairs[i].id2;

					PQP_REAL P1[3], P2[3], P3[3];
					PQP_REAL Q1[3], Q2[3], Q3[3];
					PQP_REAL Q_tmp[3];

					C2A_Tri* tri1 = model1.get()->trisConst + id1;
					C2A_Tri* tri2 = model2.get()->trisConst + id2;

					P1[0] = tri1->p1[0]; P1[1] = tri1->p1[1]; P1[2] = tri1->p1[2];
					P2[0] = tri1->p2[0]; P2[1] = tri1->p2[1]; P2[2] = tri1->p2[2];
					P3[0] = tri1->p3[0]; P3[1] = tri1->p3[1]; P3[2] = tri1->p3[2];

					Q_tmp[0] = tri2->p1[0]; Q_tmp[1] = tri2->p1[1]; Q_tmp[2] = tri2->p1[2];
					MxVpV(Q1, R2, Q_tmp, T2);
					Q_tmp[0] = tri2->p2[0]; Q_tmp[1] = tri2->p2[1]; Q_tmp[2] = tri2->p2[2];
					MxVpV(Q2, R2, Q_tmp, T2);
					Q_tmp[0] = tri2->p3[0]; Q_tmp[1] = tri2->p3[1]; Q_tmp[2] = tri2->p3[2];
					MxVpV(Q3, R2, Q_tmp, T2);

					PQP_REAL contact_points[6];
					int num_contact_points = 0;
					PQP_REAL penetration_depth = 0;
					PQP_REAL normal[3];

					TriContact(P1, P2, P3, Q1, Q2, Q3, contact_points, num_contact_points, penetration_depth, normal);
	

					for(int j = 0; j < num_contact_points; ++j)
						res.contacts.push_back(Contact(normal, &contact_points[3 * j], penetration_depth));
				}
			}

			return res;
		}
		
		// q is the relative transform of B to A
		bool isCollide(const DataVector& q) const
		{
			double R1[3][3];
			double T1[3];
			
			IdentityRotTrans(R1, T1);
			
			double R2[3][3];
			double T2[3];
			
			if(q.dim() == 6)
				ConfigEuler2RotTran(R2, T2, q);
			else if(q.dim() == 7)
				ConfigQuat2RotTrans(R2, T2, q);
			else
				return false;
				
			PQP_CollideResult result;
			
			C2A_Collide(&result, R1, T1, model1.get(), R2, T2, model2.get(), C2A_FIRST_CONTACT);
			
			return result.NumPairs() > 0;
		}


		static DataVector interpolate(const DataVector& qs, const DataVector& qt, double t)
		{
			if(qs.dim() == 6)
				return InterpConfigEuler(qs, qt, t);
			else if(qs.dim() == 7)
				return InterpConfigQuat(qs, qt, t);
			else
				return DataVector();
		}
		
		std::pair<bool, double> isCCDCollide(const DataVector& qs, const DataVector& qt, std::size_t n_dcd) const
		{
			// perform n_dcd + 1 collisions
			double R1[3][3];
			double T1[3];
			
			IdentityRotTrans(R1, T1);
			
			double R2[3][3];
			double T2[3];

			//for(int j = 0; j < qs.dim(); ++j)
			//	std::cout << qs[j] << " ";
			//std::cout << std::endl;
			//for(int j = 0; j < qt.dim(); ++j)
			//	std::cout << qt[j] << " ";
			//std::cout << std::endl;
			
			for(std::size_t i = 0; i <= n_dcd; ++i)
			{
				if(i == 0)
				{
					if(qs.dim() == 6)
						ConfigEuler2RotTran(R2, T2, qs);
					else if(qs.dim() == 7)
						ConfigQuat2RotTrans(R2, T2, qs);
					else
						return std::make_pair(false, 0);
						
					PQP_CollideResult result;
					C2A_Collide(&result, R1, T1, model1.get(), R2, T2, model2.get(), C2A_FIRST_CONTACT);
					
					if(result.NumPairs() > 0)
						return std::make_pair(true, 0);
				}
				else if(i == n_dcd)
				{
					if(qt.dim() == 6)
						ConfigEuler2RotTran(R2, T2, qt);
					else if(qt.dim() == 7)
						ConfigQuat2RotTrans(R2, T2, qt);
					else
						return std::make_pair(false, 0);
						
					PQP_CollideResult result;
					C2A_Collide(&result, R1, T1, model1.get(), R2, T2, model2.get(), C2A_FIRST_CONTACT);
					
					if(result.NumPairs() > 0)
						return std::make_pair(true, 1);
				}
				else
				{
					double t = i / (double)n_dcd;
					
					
					if(qs.dim() == 6)
					{
						DataVector q = InterpConfigEuler(qs, qt, t);
						ConfigEuler2RotTran(R2, T2, qt);
					}
					else if(qs.dim() == 7)
					{
						DataVector q = InterpConfigQuat(qs, qt, t);
						ConfigQuat2RotTrans(R2, T2, qt);
					}
					else
					{
						std::cout << "should not happen" << std::endl;
						return std::make_pair(false, 0);
					}

					//for(std::size_t j = 0; j < 3; ++j)
					//	std::cout << T2[j] << " ";
					//std::cout << ", ";
					//for(std::size_t j = 0; j < 3; ++j)
					//{
					//	for(std::size_t k = 0; k < 3; ++k)
					//		std::cout << R2[j][k] << " ";
					//	std::cout << " ";
					//}
					//std::cout << std::endl;
						
					PQP_CollideResult result;
					C2A_Collide(&result, R1, T1, model1.get(), R2, T2, model2.get(), C2A_FIRST_CONTACT);
					
					if(result.NumPairs() > 0)
						return std::make_pair(true, t);
				}
			}
			
			return std::make_pair(false, 1);
		}
		
		std::pair<bool, double> isCCDCollide(const DataVector& qs, const DataVector& qt) const
		{
			Transform identity;
			identity.Set_Translation(::Coord3D(0, 0, 0));
			identity.Set_Rotation(::Quaternion(1, 0, 0, 0));
			
			Transform trans1, trans2;
			trans1.Set_Translation(Coord3D(qs[0], qs[1], qs[2]));
			trans2.Set_Translation(Coord3D(qt[0], qt[1], qt[2]));
			
			if(qs.dim() == 6)
			{
				DataVector qs2 = ConfigEuler2Quat(qs);
				DataVector qt2 = ConfigEuler2Quat(qt);
				trans1.Set_Rotation(::Quaternion(qs2[3], qs2[4], qs2[5], qs2[6]));
				trans2.Set_Rotation(::Quaternion(qt2[3], qt2[4], qt2[5], qt2[6]));
			}
			else if(qs.dim() == 7)
			{
				trans1.Set_Rotation(::Quaternion(qs[3], qs[4], qs[5], qs[6]));
				trans2.Set_Rotation(::Quaternion(qt[3], qt[4], qt[5], qt[6]));
			}
			
			Transform first_contact_trans;
			double first_time_of_contact = 1;
			bool is_collide = false;
			
			Transform contact_transA;
			Transform contact_transB;
			double time_of_contact;
			int num_of_iteration;
			int num_of_contact;
			C2A_TimeOfContactResult dres;
			dres.last_triA = model1.get()->last_tri;
			dres.last_triB = model2.get()->last_tri;

			C2A_Solve(&identity, &identity, model1.get(), &trans1, &trans2, model2.get(),
			          contact_transA, contact_transB, time_of_contact, num_of_iteration, num_of_contact, 0.0, dres);

			if(dres.toc < 1)
			{
				is_collide = true;
				if(dres.toc < first_time_of_contact)
				{
					first_time_of_contact = dres.toc;
					first_contact_trans = contact_transB;
				}
			}
			
			return std::make_pair(is_collide, first_time_of_contact);
		}
		
	protected:
		boost::shared_ptr<C2A_Model> model1;
		boost::shared_ptr<C2A_Model> model2;
	};
	
	
	struct AABB3D
	{
		double b_min[3];
		double b_max[3];
		
		AABB3D()
		{
			b_min[0] = b_min[1] = b_min[2] = (std::numeric_limits<double>::max)();
			b_max[0] = b_max[1] = b_max[2] = -(std::numeric_limits<double>::max)();
		}
	};
	
	inline AABB3D rotate(const AABB3D& aabb, double R[3][3])
	{
		AABB3D new_aabb;
		
		
		Matrix3x3 R_(R);
		
		Coord3D v[8];
		v[0] = Coord3D(aabb.b_min[0], aabb.b_min[1], aabb.b_min[2]);
		v[1] = Coord3D(aabb.b_min[0], aabb.b_max[1], aabb.b_min[2]);
		v[2] = Coord3D(aabb.b_min[0], aabb.b_min[1], aabb.b_max[2]);
		v[3] = Coord3D(aabb.b_min[0], aabb.b_max[1], aabb.b_max[2]);
		
		v[4] = Coord3D(aabb.b_max[0], aabb.b_min[1], aabb.b_min[2]);
		v[5] = Coord3D(aabb.b_max[0], aabb.b_max[1], aabb.b_min[2]);
		v[6] = Coord3D(aabb.b_max[0], aabb.b_min[1], aabb.b_max[2]);
		v[7] = Coord3D(aabb.b_max[0], aabb.b_max[1], aabb.b_max[2]);
		
		for(int i = 0; i < 8; ++i)
			v[i] = R_ * v[i];
			
			
		for(int i = 0; i < 3; ++i)
		{
			double min_ = v[0][i];
			double max_ = v[0][i];
			
			for(int j = 1; j < 8; ++j)
			{
				if(v[j][i] < min_) min_ = v[j][i];
				else if(v[j][i] > max_) max_ = v[j][i];
			}
			
			new_aabb.b_min[i] = min_;
			new_aabb.b_max[i] = max_;
		}
		
		return new_aabb;
		
	}
	
	inline AABB3D computeAABB(C2A_Model* model)
	{
		AABB3D aabb;
		for(int i = 0; i < model->num_tris; ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				double t = model->GetTriangle(i)->p1[j];
				if(t > aabb.b_max[j]) aabb.b_max[j] = t;
				if(t < aabb.b_min[j]) aabb.b_min[j] = t;
			}
			
			for(int j = 0; j < 3; ++j)
			{
				double t = model->GetTriangle(i)->p2[j];
				if(t > aabb.b_max[j]) aabb.b_max[j] = t;
				if(t < aabb.b_min[j]) aabb.b_min[j] = t;
			}
			
			for(int j = 0; j < 3; ++j)
			{
				double t = model->GetTriangle(i)->p3[j];
				if(t > aabb.b_max[j]) aabb.b_max[j] = t;
				if(t < aabb.b_min[j]) aabb.b_min[j] = t;
			}
		}
		
		return aabb;
	}
	
}


#endif