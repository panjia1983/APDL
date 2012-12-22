#ifndef COLLISION_2D_H
#define COLLISION_2D_H

#include <cmath>
#include <cassert>
#include <utility>
#include <iostream>
#include <vector>
#include <limits>

namespace APDL
{
	class Vec2D
	{
	public:
		Vec2D()
		{
			x = 0;
			y = 0;
		}
		Vec2D(double x_, double y_) : x(x_), y(y_) {}
		
		void set(double x_, double y_)
		{
			x = x_;
			y = y_;
		}
		
		Vec2D operator - () const
		{
			return Vec2D(-x, -y);
		}
		
		Vec2D& operator += (const Vec2D& v)
		{
			x += v.x;
			y += v.y;
			return *this;
		}
		
		Vec2D& operator -= (const Vec2D& v)
		{
			x -= v.x;
			y -= v.y;
			return *this;
		}
		
		Vec2D& operator *= (double d)
		{
			x *= d;
			y *= d;
			return *this;
		}
		
		Vec2D operator * (double d) const
		{
			return Vec2D(x * d, y * d);
		}
		
		Vec2D operator + (const Vec2D& v) const
		{
			return Vec2D(x + v.x, y + v.y);
		}
		
		Vec2D operator - (const Vec2D& v) const
		{
			return Vec2D(x - v.x, y - v.y);
		}
		
		double sqrLength() const
		{
			return x * x + y * y;
		}
		
		double length() const
		{
			return sqrt(x * x + y * y);
		}
		
		Vec2D& normalize()
		{
			double len = sqrt(x * x + y * y);
			if(len > 0)
			{
				x /= len;
				y /= len;
			}
			
			return *this;
		}
		
		double x, y;
	};
	
	inline Vec2D normalize(const Vec2D& v)
	{
		Vec2D v_(v);
		v_.normalize();
		return v_;
	}
	
	inline Vec2D operator * (double s, const Vec2D& v)
	{
		return Vec2D(v.x * s, v.y * s);
	}
	
	inline double cross(const Vec2D& a, const Vec2D& b)
	{
		return a.x * b.y - a.y * b.x;
	}
	
	inline Vec2D cross(const Vec2D& a, double s)
	{
		return Vec2D(s * a.y, -s * a.x);
	}
	
	inline Vec2D cross(double s, const Vec2D& a)
	{
		return Vec2D(-s * a.y, s * a.x);
	}
	
	inline double dot(const Vec2D& a, const Vec2D& b)
	{
		return a.x * b.x + a.y * b.y;
	}
	
	inline double sqrDistance(const Vec2D& a, const Vec2D& b)
	{
		return (b - a).sqrLength();
	}
	
	inline double distance(const Vec2D& a, const Vec2D& b)
	{
		return (b - a).length();
	}
	
	inline std::ostream& operator << (std::ostream& os, const Vec2D& v)
	{
		os << v.x << " " << v.y;
		return os;
	}
	
	
	struct Mat2D
	{
		Mat2D()
		{
			col[0] = Vec2D(1, 0);
			col[1] = Vec2D(0, 1);
		}
		
		Mat2D(double angle)
		{
			double c = cos(angle);
			double s = sin(angle);
			col[0] = Vec2D(c, s);
			col[1] = Vec2D(-s, c);
		}
		
		Mat2D(const Vec2D& col1, const Vec2D& col2)
		{
			col[0] = col1;
			col[1] = col2;
		}
		
		Mat2D(double a, double b, double c, double d)
		{
			col[0] = Vec2D(a, c);
			col[1] = Vec2D(b, d);
		}
		
		Mat2D& transpose()
		{
			double temp = col[0].y;
			col[0].y = col[1].x;
			col[1].x = temp;
			return *this;
		}
		
		Mat2D& inverse()
		{
			double a = col[0].x;
			double b = col[1].x;
			double c = col[0].y;
			double d = col[1].y;
			
			double det = a * d - b * c;
			assert(det != 0);
			det = 1.0 / det;
			col[0].x = det * d;
			col[0].y = -det * c;
			col[1].x = -det * b;
			col[1].y = det * a;
			
			return *this;
		}
		
		void set(double angle)
		{
			double c = cos(angle);
			double s = sin(angle);
			col[0].x = c;
			col[0].y = s;
			col[1].x = -s;
			col[1].y = c;
		}
		
		Vec2D operator * (const Vec2D& v) const
		{
			return Vec2D(col[0].x * v.x + col[1].x * v.y, col[0].y * v.x + col[1].y * v.y);
		}
		
		Vec2D col[2];
	};
	
	inline std::ostream& operator << (std::ostream& os, const Mat2D& m)
	{
		os << m.col[0].x << " " << m.col[1].x << std::endl;
		os << m.col[0].y << " " << m.col[1].y;
		return os;
	}
	
	
	inline Mat2D transpose(const Mat2D& m)
	{
		return Mat2D(m.col[0].x, m.col[0].y, m.col[1].x, m.col[1].y);
	}
	
	inline Mat2D inverse(const Mat2D& m)
	{
		Mat2D r(m);
		r.inverse();
		return r;
	}
	
	
	struct Transform2D
	{
		Mat2D R;
		Vec2D t;
		
		Transform2D()
		{
			R = Mat2D(1, 0, 0, 1);
			t = Vec2D(0, 0);
		}
		
		Transform2D(const Mat2D& R_, const Vec2D& t_) : R(R_), t(t_)
		{}
		
		Vec2D transform(const Vec2D& v) const
		{
			return R * v + t;
		}
		
		Vec2D untransform(const Vec2D& v) const
		{
			return inverse(R) * (v - t);
		}
		
		Vec2D transform_dir(const Vec2D& v) const
		{
			return R * v;
		}
		
		Vec2D untransform_dir(const Vec2D& v) const
		{
			return inverse(R) * v;
		}
	};
	
	
	
	struct SimplexVertex
	{
		Vec2D point1; // support point in shape 1
		Vec2D point2; // support point in shape 2
		Vec2D point;  // point2 - point1
		double u;     // unnormalized barycentric coordinate for closest point
		
		int index1;   // point1 index
		int index2;   // point2 index
	};
	
	struct Simplex
	{
		SimplexVertex verts[3];
		int count;
		double divisor; // denominator to normalize barycentric coordinates
		
		// Find direction toward origin
		Vec2D getSearchDirection() const;
		
		// Compute closest simplex point to origin
		Vec2D getClosestPoint() const;
		
		// Compute witness points
		void getWitnessPoints(Vec2D& point1, Vec2D& point2) const;
		void solve2(const Vec2D& P);
		void solve3(const Vec2D& P);
		
		Simplex()
		{
			count = 0;
			divisor = 1;
		}
	};
	
	
	struct AABB2D
	{
		double x_min, x_max;
		double y_min, y_max;
		
		AABB2D(double x1, double x2, double y1, double y2)
		{
			x_min = x1;
			x_max = x2;
			y_min = y1;
			y_max = y2;
		}
	};
	
	inline AABB2D rotate(const AABB2D& aabb, double angle)
	{
		Vec2D v[4];
		v[0] = Vec2D(aabb.x_min, aabb.y_min);
		v[1] = Vec2D(aabb.x_min, aabb.y_max);
		v[2] = Vec2D(aabb.x_max, aabb.y_min);
		v[3] = Vec2D(aabb.x_max, aabb.y_max);
		
		Mat2D mat(angle);
		for(int i = 0; i < 4; ++i)
			v[i] = mat * v[i];
			
		double x_min, x_max, y_min, y_max;
		x_min = x_max = v[0].x;
		y_min = y_max = v[0].y;
		for(int i = 1; i < 4; ++i)
		{
			if(v[i].x < x_min) x_min = v[i].x;
			else if(v[i].x > x_max) x_max = v[i].x;
			if(v[i].y < y_min) y_min = v[i].y;
			else if(v[i].y > y_max) y_max = v[i].y;
		}
		
		return AABB2D(x_min, x_max, y_min, y_max);
	}
	
	class Shape2D
	{
	public:
		virtual ~Shape2D() {}
		
		virtual std::pair<Vec2D, int> getSupport(const Vec2D& d) const = 0;
		
		AABB2D getAABB() const;
	};
	
	
	class Polygon : public Shape2D
	{
	public:
		Polygon()
		{
		}
		
		std::size_t size() const
		{
			return points.size();
		}
		
		std::pair<Vec2D, int> getSupport(const Vec2D& d) const
		{
			int best_id = 0;
			double best_v = dot(points[0], d);
			for(std::size_t i = 1; i < points.size(); ++i)
			{
				double v = dot(points[i], d);
				if(v > best_v)
				{
					best_id = i;
					best_v = v;
				}
			}
			
			return std::make_pair(points[best_id], best_id);
		}
		
		AABB2D getAABB() const
		{
			double x_min, x_max, y_min, y_max;
			x_min = x_max = points[0].x;
			y_min = y_max = points[0].y;
			
			for(std::size_t i = 1; i < points.size(); ++i)
			{
				double x = points[i].x;
				double y = points[i].y;
				
				if(x < x_min) x_min = x;
				else if(x > x_max) x_max = x;
				if(y < y_min) y_min = y;
				else if(y > y_max) y_max = y;
			}
			
			return AABB2D(x_min, x_max, y_min, y_max);
		}
		
		std::pair<Vec2D, double> getCircle() const
		{
			Vec2D c = points[0];
			for(std::size_t i = 1 ; i < points.size(); ++i)
				c += points[i];
				
			c *= (1.0 / points.size());
			
			double r_max = 0;
			
			for(std::size_t i = 0; i < points.size(); ++i)
			{
				double r = (points[i] - c).sqrLength();
				if(r > r_max) r_max = r;
			}
			
			return std::make_pair(c, std::sqrt(r_max));
		}
		
		std::vector<Vec2D> points;
	};
	
	
	struct GJKResult
	{
		Vec2D point1; // closest point on shape1
		Vec2D point2; // closest point on shape2
		double distance;
		int iterations;
	};
	
	
	inline GJKResult doGJK(const Shape2D& s1, const Transform2D& tf1,
	                       const Shape2D& s2, const Transform2D& tf2)
	{
		const double epsilon = 1e-4;
		Simplex simplex;
		std::pair<Vec2D, int> sup1 = s1.getSupport(tf1.untransform_dir(Vec2D(-1, 0)));
		std::pair<Vec2D, int> sup2 = s2.getSupport(tf2.untransform_dir(Vec2D(1, 0)));
		
		simplex.verts[0].point1 = tf1.transform(sup1.first);
		simplex.verts[0].point2 = tf2.transform(sup2.first);
		simplex.verts[0].point = simplex.verts[0].point2 - simplex.verts[0].point1;
		simplex.verts[0].u = 1;
		simplex.count = 1;
		
		// Main iteration loop
		int k_max_iters = 20;
		int iter = 0;
		while(iter < k_max_iters)
		{
			// determine the closest point on the simplex and remove unused vertices
			switch(simplex.count)
			{
			case 1:
				break;
			case 2:
				simplex.solve2(Vec2D(0, 0));
				break;
			case 3:
				simplex.solve3(Vec2D(0, 0));
				break;
			default:
				assert(0);
			}
			
			// if we have 3 points, then the origin is in the corresponding triangle
			if(simplex.count == 3)
				break;
				
			// get search direction
			Vec2D d = simplex.getSearchDirection();
			
			// ensure the seach direction is non-zero
			if(dot(d, d) == 0)
				break;
				
			// compute a tentative new simplex using support points
			SimplexVertex* vertex = simplex.verts + simplex.count;
			
			std::pair<Vec2D, int> sup1 = s1.getSupport(tf1.untransform_dir(-d));
			std::pair<Vec2D, int> sup2 = s2.getSupport(tf2.untransform_dir(d));
			
			vertex->point1 = tf1.transform(sup1.first);
			vertex->point2 = tf2.transform(sup2.first);
			vertex->point = vertex->point2 - vertex->point1;
			
			
			++iter;
			
			// check for duplicate support points, this is the main terminate criteria
			bool duplicate = false;
			for(int i = 0; i < simplex.count; ++i)
			{
				if(distance(vertex->point1, simplex.verts[i].point1) < epsilon && distance(vertex->point2, simplex.verts[i].point2) < epsilon)
				{
					duplicate = true;
					break;
				}
			}
			
			if(duplicate)
				break;
				
			++simplex.count;
		}
		
		GJKResult result;
		
		simplex.getWitnessPoints(result.point1, result.point2);
		result.distance = distance(result.point1, result.point2);
		result.iterations = iter;
		
		return result;
	}
	
	
	inline GJKResult doGJK(const Polygon& s1, const Transform2D& tf1,
	                       const Polygon& s2, const Transform2D& tf2,
						   Simplex* last_simplex = NULL)
	{
		Simplex simplex;
		std::pair<Vec2D, int> sup1 = s1.getSupport(tf1.untransform_dir(Vec2D(-1, 0)));
		std::pair<Vec2D, int> sup2 = s2.getSupport(tf2.untransform_dir(Vec2D(1, 0)));
		
		simplex.verts[0].index1 = sup1.second;
		simplex.verts[0].index2 = sup2.second;
		
		simplex.verts[0].point1 = tf1.transform(sup1.first);
		simplex.verts[0].point2 = tf2.transform(sup2.first);
		simplex.verts[0].point = simplex.verts[0].point2 - simplex.verts[0].point1;
		simplex.verts[0].u = 1;
		simplex.count = 1;
		
		// These store the vertices of the last simplex so that we can check for duplicates and prevent cycling.
		int save1[3], save2[3];
		int save_count = 0;
		
		// Main iteration loop
		const int k_max_iters = 20;
		int iter = 0;
		while(iter < k_max_iters)
		{
			save_count = simplex.count;
			for(int i = 0; i < save_count; ++i)
			{
				save1[i] = simplex.verts[i].index1;
				save2[i] = simplex.verts[i].index2;
			}
			
			// determine the closest point on the simplex and remove unused vertices
			switch(simplex.count)
			{
			case 1:
				break;
			case 2:
				simplex.solve2(Vec2D(0, 0));
				break;
			case 3:
				simplex.solve3(Vec2D(0, 0));
				break;
			default:
				assert(0);
			}
			
			// if we have 3 points, then the origin is in the corresponding triangle
			if(simplex.count == 3)
				break;
				
			// get search direction
			Vec2D d = simplex.getSearchDirection();
			
			// ensure the seach direction is non-zero
			if(dot(d, d) == 0)
				break;
				
			// compute a tentative new simplex using support points
			SimplexVertex* vertex = simplex.verts + simplex.count;
			
			std::pair<Vec2D, int> sup1 = s1.getSupport(tf1.untransform_dir(-d));
			std::pair<Vec2D, int> sup2 = s2.getSupport(tf2.untransform_dir(d));
			vertex->index1 = sup1.second;
			vertex->index2 = sup2.second;
			vertex->point1 = tf1.transform(sup1.first);
			vertex->point2 = tf2.transform(sup2.first);
			vertex->point = vertex->point2 - vertex->point1;
			
			++iter;
			
			// check for duplicate support points, this is the main terminate criteria
			bool duplicate = false;
			for(int i = 0; i < save_count; ++i)
			{
				if(vertex->index1 == save1[i] && vertex->index2 == save2[i])
				{
					duplicate = true;
					break;
				}
			}
			
			if(duplicate)
				break;
				
			++simplex.count;
		}
		
		GJKResult result;
		
		simplex.getWitnessPoints(result.point1, result.point2);
		result.distance = distance(result.point1, result.point2);
		
		result.iterations = iter;

		if(last_simplex) *last_simplex = simplex;
		
		return result;
	}
	
	struct PolytopeEdge
	{
		int index1, index2;
		PolytopeEdge* next, *prev;
		
		double distsq;
		Vec2D dir;
		
		PolytopeEdge(int id1, int id2)
		{
			next = prev = NULL;
			distsq = (std::numeric_limits<double>::max)();
			index1 = id1;
			index2 = id2;
		}
		
	};
	
	struct Polytope
	{
		std::vector<SimplexVertex> verts;
		int count;
		
		PolytopeEdge* head, *tail;
		
		Polytope(const Simplex& simplex)
		{
			count = simplex.count;
			verts.resize(count);
			for(int i = 0; i < count; ++i)
				verts[i] = simplex.verts[i];
				
			head = tail = NULL;
			
			if(count == 2)
			{
				insertEdge(tail, new PolytopeEdge(0, 1));
				insertEdge(tail, new PolytopeEdge(1, 0));
			}
			else if(count == 3)
			{
				Vec2D a = simplex.verts[0].point;
				Vec2D b = simplex.verts[1].point;
				Vec2D c = simplex.verts[2].point;
				Vec2D ab = b - a;
				Vec2D bc = c - b;
				
				if(cross(ab, bc) > 0)
				{
					insertEdge(tail, new PolytopeEdge(0, 1));
					insertEdge(tail, new PolytopeEdge(1, 2));
					insertEdge(tail, new PolytopeEdge(2, 0));
				}
				else
				{
					insertEdge(tail, new PolytopeEdge(0, 2));
					insertEdge(tail, new PolytopeEdge(2, 1));
					insertEdge(tail, new PolytopeEdge(1, 0));
				}
			}
		}

		~Polytope()
		{
			PolytopeEdge* p = head;
			while(p != tail)
			{
				PolytopeEdge* t = p;
				p = p->next;
				delete t;
			}

			delete tail;
		}
		
		void insertEdge(PolytopeEdge* pre_e,
		                PolytopeEdge* new_e)
		{
			if(head == NULL)
			{
				head = new_e;
				tail = new_e;
				new_e->next = new_e;
				new_e->next = new_e;
			}
			else
			{
				new_e->prev = pre_e;
				new_e->next = pre_e->next;
				new_e->next->prev = new_e;
				pre_e->next = new_e;
				
				if(pre_e == tail)
					tail = new_e;
			}
		}
		
		void deleteEdge(PolytopeEdge* e)
		{
			if(e == head)
				head = e->next;
				
			if(e == tail)
				tail = e->prev;
				
			e->prev->next = e->next;
			e->next->prev = e->prev;
		}
		
		PolytopeEdge* getClosestEdge()
		{
			PolytopeEdge* first_edge = head;
			
			if(first_edge->distsq == (std::numeric_limits<double>::max)())
			{
				Vec2D a = verts[first_edge->index1].point;
				Vec2D b = verts[first_edge->index2].point;
				Vec2D ab = b - a;
				
				double v = -dot(ab, a);
				if(v <= 0)
				{
					Vec2D cp(a.x, a.y);
					first_edge->distsq = cp.sqrLength();
					first_edge->dir = cp;
				}
				else
				{
					double u = dot(ab, b);
					if(u <= 0)
					{
						Vec2D cp(b.x, b.y);
						first_edge->distsq = cp.sqrLength();
						first_edge->dir = cp;
					}
					else
					{
						double s = 1 / ab.sqrLength();
						double t = v * s;
						Vec2D cp = (1 - t) * a + t * b;
						first_edge->distsq = cp.sqrLength();
						first_edge->dir = Vec2D(ab.y, -ab.x);
					}
				}
			}
			
			
			PolytopeEdge* closest_edge = first_edge;
			
			
			for(PolytopeEdge* edge = first_edge->next;
			        edge != head;
			        edge = edge->next)
			{
				if(edge->distsq == (std::numeric_limits<double>::max)())
				{
					Vec2D a = verts[edge->index1].point;
					Vec2D b = verts[edge->index2].point;
					Vec2D ab = b - a;
					
					double v = -dot(ab, a);
					if(v <= 0)
					{
						Vec2D cp(a);
						edge->distsq = cp.sqrLength();
						edge->dir = cp;
					}
					else
					{
						double u = dot(ab, b);
						if(u <= 0)
						{
							Vec2D cp(b);
							edge->distsq = cp.sqrLength();
							edge->dir = cp;
						}
						else
						{
							double s = 1 / ab.sqrLength();
							double t = v * s;
							Vec2D cp = (1 - t) * a + t * b;
							edge->distsq = cp.sqrLength();
							edge->dir = Vec2D(ab.y, -ab.x);
						}
					}
				}
				
				if(edge->distsq > 1e-4 && edge->distsq < closest_edge->distsq)
				{
					closest_edge = edge;
				}
			}
			
			return closest_edge;
		}
	};
	
	
	inline Polytope doEPA(const Polygon& s1, const Transform2D& tf1,
	                      const Polygon& s2, const Transform2D& tf2,
	                      const Simplex& simplex)
	{
		Polytope polytope(simplex);
		
		std::vector<int> save1, save2;
		
		int max_iters = 20;
		int iter = 0;
		while(iter < max_iters)
		{
			int save_count = polytope.verts.size();
			save1.resize(save_count);
			save2.resize(save_count);
			for(int i = 0; i < save_count; i++)
			{
				save1[i] = polytope.verts[i].index1;
				save2[i] = polytope.verts[i].index2;
			}
			
			PolytopeEdge* edge = polytope.getClosestEdge();
			
			Vec2D d = edge->dir;
			
			if(dot(d, d) == 0)
				break;
				
			std::pair<Vec2D, int> sup1 = s1.getSupport(tf1.untransform_dir(-d));
			Vec2D p1 = tf1.transform(sup1.first);
			std::pair<Vec2D, int> sup2 = s2.getSupport(tf2.untransform_dir(d));
			Vec2D p2 = tf2.transform(sup2.first);
			Vec2D p = p2 - p1;
			
			SimplexVertex v1 = polytope.verts[edge->index1];
			SimplexVertex v2 = polytope.verts[edge->index2];
			
			if((v1.index1 == sup1.second && v1.index2 == sup2.second) || (v2.index1 == sup1.second && v2.index2 == sup2.second))
				break;
				
			SimplexVertex new_v;
			new_v.index1 = sup1.second;
			new_v.index2 = sup2.second;
			new_v.point1 = p1;
			new_v.point2 = p2;
			new_v.point = p;
			
			polytope.verts.push_back(new_v);
			int new_index = polytope.verts.size() - 1;
			PolytopeEdge* prev_e = edge->prev;
			PolytopeEdge* next_e = edge->next;
			polytope.deleteEdge(edge);
			delete edge;
			edge = NULL;
			
			polytope.insertEdge(prev_e, new PolytopeEdge(prev_e->index2, new_index));
			polytope.insertEdge(prev_e->next, new PolytopeEdge(new_index, next_e->index1));
			
			bool duplicate = false;
			for(int i = 0; i < save_count; ++i)
			{
				if(new_v.index1 == save1[i] && new_v.index2 == save2[i])
				{
					duplicate = true;
					break;
				}
			}
			
			if(duplicate)
				break;
		}
		
		return polytope;
	}

	struct Edge
	{
		Vec2D v1, v2;
	};

	inline Edge findSeparateEdge(const Polygon& s, const Transform2D& tf, const Vec2D& n)
	{
		Vec2D n_local = tf.untransform_dir(n);
		int index = s.getSupport(n_local).second;
		int index_prev = (index + s.size() - 1) % (s.size());
		int index_next = (index + 1) % (s.size());

		const Vec2D& v = s.points[index];
		const Vec2D& v_prev = s.points[index_prev];
		const Vec2D& v_next = s.points[index_next];
		Vec2D l = v - v_next;
		Vec2D r = v - v_prev;

		Edge edge;
		
		if(dot(r, n_local) <= dot(l, n_local))
		{
			edge.v1 = tf.transform(v_prev);
			edge.v2 = tf.transform(v);
			return edge;
		}

		edge.v1 = tf.transform(v);
		edge.v2 = tf.transform(v_next);
		return edge;
	}

	inline std::vector<Vec2D> clipLineSegment(const Vec2D& v1, const Vec2D& v2, const Vec2D& n, double o)
	{
		double d1 = dot(n, v1) - o;
		double d2 = dot(n, v2) - o;
		
		std::vector<Vec2D> cp;
		if(d1 >= 0) cp.push_back(v1);
		if(d2 >= 0) cp.push_back(v2);
		if(d1 * d2 < 0)
		{
			Vec2D delta = v2 - v1;
			Vec2D p = v1 + delta * (d1 / (d1 - d2));
			cp.push_back(p);
		}

		return cp;
	}

	struct EPAContact
	{
		Vec2D normal;
		Vec2D point;
		double penetration;

		EPAContact(const Vec2D& point_, const Vec2D& normal_, double penetration_)
		{
			point = point_;
			normal = normal_;
			penetration = penetration_;
		}
	};

	inline std::vector<EPAContact> computeContactPoints(const Polygon& s1, const Transform2D& tf1,
														const Polygon& s2, const Transform2D& tf2,
														const Vec2D & n)
	{
		std::vector<EPAContact> cp;
		Edge e1 = findSeparateEdge(s1, tf1, n);
		Edge e2 = findSeparateEdge(s2, tf2, -n);

		Vec2D e1d = e1.v2 - e1.v1;
		Vec2D e2d = e2.v2 - e2.v1;

		Edge ref, inc;
		Vec2D ref_n;
		bool flip;


		double en1 = abs(dot(e1d, n));
		double en2 = abs(dot(e2d, n));
		if(en1 <= en2)
		{
			ref = e1;
			ref_n = normalize(e1d);
			inc = e2;
			flip = true;
		}
		else
		{
			ref = e2;
			ref_n = normalize(e2d);
			inc = e1;
			flip = false;
		}

		std::vector<Vec2D> v;

		double o1 = dot(ref_n, ref.v1);
		v = clipLineSegment(inc.v1, inc.v2, ref_n, o1);
		if(v.size() < 2)
			return cp;

		double o2 = -dot(ref_n, ref.v2);
		v = clipLineSegment(v[0], v[1], -ref_n, o2);
		if(v.size() < 2)
			return cp;

		Vec2D ref_perp(-ref_n.y, ref_n.x);


		double o3 = dot(ref_perp, ref.v1);
		double depth0 = dot(ref_perp, v[0]) - o3;
		double depth1 = dot(ref_perp, v[1]) - o3;

		
		if(depth0 > 0)
			cp.push_back(EPAContact(v[0], n, -depth0));

		if(depth1 > 0)
			cp.push_back(EPAContact(v[1], n, -depth1));

		return cp;
	}
	
	struct EPAResult
	{
		bool distance;

		std::vector<EPAContact> contacts;
	};

	inline EPAResult doGJKEPA(const Polygon& s1, const Transform2D& tf1,
							  const Polygon& s2, const Transform2D& tf2)
	{
		Simplex last_simplex;
		GJKResult gjk_res = doGJK(s1, tf1, s2, tf2, &last_simplex);

		EPAResult epa_res;
		if(gjk_res.distance > 0)
		{
			epa_res.distance = gjk_res.distance;
			return epa_res;
		}

		Polytope polytope = doEPA(s1, tf1, s2, tf2, last_simplex);
		Vec2D normal = polytope.getClosestEdge()->dir;
		normal.normalize();

		std::vector<EPAContact> info = computeContactPoints(s1, tf1, s2, tf2, -normal);

		epa_res.distance = gjk_res.distance;
		epa_res.contacts = info;

		return epa_res;

	}
	
}

#endif