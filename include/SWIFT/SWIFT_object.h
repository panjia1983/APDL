/*************************************************************************\

  Copyright 2001 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify OR distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             S. Ehmann, M. Lin
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919) 962-1749

  EMail:               geom@cs.unc.edu
                       ehmann@cs.unc.edu
                       lin@cs.unc.edu

\**************************************************************************/


//////////////////////////////////////////////////////////////////////////////
//
// SWIFT_object.h
//
// Description:
//      Classes to manage objects in the scene.
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _SWIFT_OBJECT_H_
#define _SWIFT_OBJECT_H_

#include <iostream>
using namespace std;
#include <SWIFT_config.h>
#include <SWIFT_common.h>
#include <SWIFT_linalg.h>
#include <SWIFT_mesh.h>
#include <SWIFT_pair.h>
#include <SWIFT_boxnode.h>

//////////////////////////////////////////////////////////////////////////////
// SWIFT_Object
//
// Description:
//      Class to manage an object.  The object cannot be used until it has
//  been initialized by calling the Initialize function.
//////////////////////////////////////////////////////////////////////////////
class SWIFT_Object
{
public:
	SWIFT_Object();
	~SWIFT_Object();
	
	// Get functions
	const SWIFT_Transformation& Transformation() const
	{
		return transform;
	}
	const SWIFT_Transformation& Transformation_Source() const
	{
		return transform_s;    //zhangxy, Feb 22, 2006
	}
	const SWIFT_Transformation& Transformation_Target() const
	{
		return transform_t;    //zhangxy, Feb 22, 2006
	}
	
	SWIFT_Tri_Mesh* Mesh()
	{
		return mesh;
	}
	SWIFT_Triple& Center_Of_Mass()
	{
		return mesh->Center_Of_Mass();
	}
	SWIFT_Real Radius()
	{
		return mesh->Radius();
	}
	SWIFT_Real AngularRadius()
	{
		return mesh->AngularRadius();    //zhangxy 2005-12-06
	}
	void Velocity()
	{
		Velocities(transform_s, transform_t, 1.0, v, w);    //zhangxy 2006-02-11
	}
	SWIFT_Triple AngularVel()
	{
		return w;    //zhangxy 2006-02-11
	}
	SWIFT_Triple LinearVel()
	{
		return v;    //zhangxy 2006-02-11
	}
	SWIFT_Triple RootAngularVel()
	{
		return rootw;    //zhangxy 2006-09-13
	}
	SWIFT_Real RootAngularRadius()
	{
		return rootr;    //zhangxy 2006-30-30
	}
	SWIFT_Triple RootLinearVel()
	{
		return rootv;    //zhangxy 2006-09-13
	}
	SWIFT_Real AccumAngularVelLen()
	{
		return accumwlen;    //zhangxy 2006-09-13
	}
	SWIFT_Real AngularVelLen()
	{
		return wlen;    //zhangxy 2006-10-30
	}
	SWIFT_Real AccumLinearVelLen()
	{
		return accumvlen;    //zhangxy 2006-09-13
	}
	SWIFT_Real LinearVelLen()
	{
		return vlen;    //zhangxy 2006-10-30
	}
	SWIFT_Real AccumAngularRadius()
	{
		return accumar;    //zhangxy 2006-09-13
	}
	SWIFT_Real AccumAngularVelLenDotRadius()
	{
		return accumwlendotr;    //zhangxy 2006-10-30
	}
	
	bool Fixed()
	{
		return fixed;
	}
	bool Use_Cube()
	{
		return cube;
	}
	SWIFT_Box_Node* Min_Box_Node(int axis)
	{
		return min_bns + axis;
	}
	SWIFT_Box_Node* Max_Box_Node(int axis)
	{
		return max_bns + axis;
	}
	void Get_Box_Nodes(int i,
	                   SWIFT_Box_Node** min_0, SWIFT_Box_Node** max_0,
	                   SWIFT_Box_Node** min_1, SWIFT_Box_Node** max_1,
	                   SWIFT_Box_Node** min_2, SWIFT_Box_Node** max_2);
	SWIFT_Array<SWIFT_Pair>& Pairs()
	{
		return pairs;
	}
	int Num_Pairs()
	{
		return pairs.Length();
	}
	
	// Set functions
	void Set_Id(int i);
	int oId; //zhangxy, Oct 25, 2006
	bool isRoot;
	
	// Initialization functions.  Should only be called once after this
	// object has been constructed.
	
	// Single piece version
	void Initialize(SWIFT_Tri_Mesh* m, bool is_fixed, bool use_cube,
	                SWIFT_Real box_enl_rel, SWIFT_Real box_enl_abs,
	                bool copy);
	                
	// Multiple piece version.  If box_enl_rel or box_enl_abs are NULL then the
	// default of zero enlargement is done (for rel or abs).
	void Initialize(SWIFT_Tri_Mesh** ms, bool is_fixed, const bool* use_cube,
	                const SWIFT_Real* box_enl_rel,
	                const SWIFT_Real* box_enl_abs, const bool* copy);
	                
	// Object update functions
	inline void Set_Transformation_No_Boxes(const SWIFT_Real* R,
	                                        const SWIFT_Real* T);
	inline void Set_Transformation_No_Boxes(const SWIFT_Real* RT);
	inline void Set_Transformation(const SWIFT_Real* R, const SWIFT_Real* T);
	inline void Set_Transformation(const SWIFT_Real* RT);
	
	//zhangxy Feb 10, 2006-----------
	inline void Set_Source_Transformation(const SWIFT_Real* R, const SWIFT_Real* T);
	inline void Set_Source_Transformation(const SWIFT_Real* RT);
	inline void Set_Target_Transformation(const SWIFT_Real* R, const SWIFT_Real* T);
	inline void Set_Target_Transformation(const SWIFT_Real* RT);
	inline void Set_RootAngularVelocity(const SWIFT_Real* rav);   //September 12, 2006
	inline void Set_RootAngularRadius(const SWIFT_Real ar);   //September 12, 2006
	inline void Set_RootTranslationalVelocity(const SWIFT_Real* rtv);   //September 12, 2006
	inline void Set_AccumAngularVelocityLen(const SWIFT_Real aavl);   //September 12, 2006
	inline void Set_LinearVelocityLen(const SWIFT_Real lvl);   //October 30, 2006
	inline void Set_AngularVelocityLen(const SWIFT_Real avl);   //October 30, 2006
	inline void Set_AccumTranslationalVelocityLen(const SWIFT_Real atvl);   //September 12, 2006
	inline void Set_AccumAngularRadius(const SWIFT_Real aar);   //September 12, 2006
	inline void Set_AccumAngularVelocityLenDotRadius(const SWIFT_Real aavlr); //October 30, 2006
	inline void Integrate(const SWIFT_Real t);
	inline void SkewMotion(const SWIFT_Real dt);
	inline void SkewMotionSetup();
	inline void Reset_Transformation();
	//-----------zhangxy Feb 10, 2006
	
private:
	inline void Update_Boxes();
	
	bool fixed;
	int id;
	SWIFT_Transformation transform;
	SWIFT_Transformation transform_s; //zhangxy Feb 10, 2006
	SWIFT_Transformation transform_t; //zhangxy Feb 10, 2006
	SWIFT_Triple v; //zhangxy Feb 11, 2006
	SWIFT_Triple w; //zhangxy Feb 11, 2006
	SWIFT_Matrix33 Ra, Rb, Rc; //zhangxy August 07, 2006
	SWIFT_Triple rootv; //zhangxy Sept 13, 2006
	SWIFT_Triple rootw; //zhangxy Sept 13, 2006
	SWIFT_Real rootr; //zhangxy Oct 30, 2006
	SWIFT_Real accumvlen; //zhangxy Sept 13, 2006
	SWIFT_Real accumwlen; //zhangxy Sept 13, 2006
	SWIFT_Real vlen; //zhangxy Oct 30, 2006
	SWIFT_Real wlen; //zhangxy Oct 30, 2006
	SWIFT_Real accumar; //zhangxy Sept 13, 2006
	SWIFT_Real accumwlendotr; //zhangxy Oct 30, 2006
	
	SWIFT_Tri_Mesh* mesh;
	
	// AABB nodes
	SWIFT_Box_Node min_bns[3];
	SWIFT_Box_Node max_bns[3];
	
	// AABB parameters
	bool cube;
	SWIFT_Real enlargement;
	
	//   cube parameters
	SWIFT_Real radius;
	
	//   dynamic box parameters
	SWIFT_Tri_Edge* min_es[3];
	SWIFT_Tri_Edge* max_es[3];
	
	// Pairs
	SWIFT_Array<SWIFT_Pair> pairs;
};


///////////////////////////////////////////////////////////////////////////////
// Inline functions
///////////////////////////////////////////////////////////////////////////////

inline void SWIFT_Object::Set_Transformation_No_Boxes(const SWIFT_Real* R,
        const SWIFT_Real* T)
{
	transform.Set_Value(R, T);
}

inline void SWIFT_Object::Set_Transformation_No_Boxes(const SWIFT_Real* RT)
{
	transform.Set_Value(RT);
}

inline void SWIFT_Object::Set_Transformation(const SWIFT_Real* R,
        const SWIFT_Real* T)
{
	Set_Transformation_No_Boxes(R, T);
	
	// Update the bounding boxes
	Update_Boxes();
}

inline void SWIFT_Object::Set_Transformation(const SWIFT_Real* RT)
{
	Set_Transformation_No_Boxes(RT);
	
	// Update the bounding boxes
	Update_Boxes();
}

inline void SWIFT_Object::Set_Source_Transformation(const SWIFT_Real* R,
        const SWIFT_Real* T)
{
	transform_s.Set_Value(R, T);
}

inline void SWIFT_Object::Set_Source_Transformation(const SWIFT_Real* RT)
{
	transform_s.Set_Value(RT);
}

inline void SWIFT_Object::Set_Target_Transformation(const SWIFT_Real* R,
        const SWIFT_Real* T)
{
	transform_t.Set_Value(R, T);
}

inline void SWIFT_Object::Set_Target_Transformation(const SWIFT_Real* RT)
{
	transform_t.Set_Value(RT);
}


inline void SWIFT_Object::Set_RootAngularVelocity(const SWIFT_Real* rav)
{
	rootw.Set_Value(rav);
}

inline void SWIFT_Object::Set_RootAngularRadius(const SWIFT_Real ar)
{
	rootr = ar;
}

inline void SWIFT_Object::Set_RootTranslationalVelocity(const SWIFT_Real* rtv)
{
	rootv.Set_Value(rtv);
}

inline void SWIFT_Object::Set_AccumAngularVelocityLen(const SWIFT_Real aavl)
{
	accumwlen = aavl;
}

inline void SWIFT_Object::Set_LinearVelocityLen(const SWIFT_Real lvl)
{
	vlen = lvl;
}

inline void SWIFT_Object::Set_AngularVelocityLen(const SWIFT_Real avl)
{
	wlen = avl;
}

inline void SWIFT_Object::Set_AccumAngularVelocityLenDotRadius(const SWIFT_Real aavlr)
{
	accumwlendotr = aavlr;
}

inline void SWIFT_Object::Set_AccumTranslationalVelocityLen(const SWIFT_Real atvl)
{
	accumvlen = atvl;
}

inline void SWIFT_Object::Set_AccumAngularRadius(const SWIFT_Real aar)
{
	accumar = aar;
}

inline void SWIFT_Object::Reset_Transformation()
{
	transform.Set_Value(transform_s);
}

#define ANGULARMOTIONTRESHOLD 0.78539815 //PI/4
inline void SWIFT_Object::Integrate(const SWIFT_Real dt)
{
	if(dt <= 1.0)
	{
		transform.Set_Translation(transform_s.Translation() + v * dt);
		
#ifdef QUATERNION_DERIVATIVE
		Quaternion orn = getRotation();
		orn += ((*v) * orn) * (dt * 0.5f);
		orn.normalize();
#else
		//exponential map
		SWIFT_Triple axis;
		double	angle = w.Length();
		if((angle * dt) > ANGULARMOTIONTRESHOLD)
		{
			angle = (ANGULARMOTIONTRESHOLD) / dt;
		}
		
		//limit the angular motion
		if(angle < 0.001f)
		{
			// use Taylor's expansions of sync function
			axis   = w * (0.5f * dt - (dt * dt * dt) * (0.020833333333) * angle * angle);
		}
		else
		{
			// sync(fAngle) = sin(c*fAngle)/t
			axis   = w * (sin(0.5f * angle * dt) / angle);
		}
		SWIFT_Quaternion dorn(axis.X(), axis.Y(), axis.Z(), cos(angle * dt * 0.5f));
		SWIFT_Quaternion orn0 = transform_s.Quaternion();
		
		SWIFT_Quaternion predictedOrn = dorn % orn0;
#endif
		transform.Set_Rotation(predictedOrn);
	}
	else
	{
		transform = transform_t;
	}
	
}

inline void SWIFT_Object::SkewMotion(const SWIFT_Real dt)
{
	SWIFT_Real angle = w.Length() * dt;
	SWIFT_Real T = tan(0.5 * angle);
	SWIFT_Real T2 = T * T;
	SWIFT_Real denom = 1.0 / (1 + T2);
	SWIFT_Real C = (1 - T2) * denom;
	SWIFT_Real S = 2.0 * T * denom;
	
	SWIFT_Matrix33 Rt;
	Rt = Ra * C + Rb * S + Rc;
	transform.Set_Rotation(Rt);
	transform.Set_Translation(transform_s.Translation() + v * dt);
}


inline void SWIFT_Object::SkewMotionSetup()
{
	// Init for the rotation matrix of a link in the reference frame of the parent link
	//
	// u[3] contains the unit rotation axis
	// r0 contains the rotation matrix at time 0
	// In the end, the rotation matrix has the form : a.cos(wt)+b.sin(wt)+c
	//
	// a=(I-u.ut).r0
	// b=u*.r0 (u.star, a matrix such as u*.x=u^x for each x)
	// c=u.ut.r0
	
	// *********
	// c=u.ut.r0
	// *********
	
	// c=u.(r0t.u)t so we compute r0t.u first
	
	SWIFT_Triple u, R0tu;
	SWIFT_Matrix33 R0, R0t, ustar;
	u = w;
	u.Normalize();
	
	R0 = transform_s.Rotation();
	R0t = R0;
	R0t.Transpose();
	R0tu = R0t * u;
	Rc.Set_Value_Rows(u.X()*R0tu, u.Y()*R0tu, u.Z()*R0tu);
	
	// a=r0-u.ut.r0
	
	Ra = R0 - Rc;
	
	// b=u*.r0
	//
	//     ( 0    -u[2] u[1]  )
	// u*= ( u[2]   0   -u[0] )
	//     ( -u[1]  u[0]   0  )
	
	ustar.Set_Value_Rows(SWIFT_Triple(0,		-u.Z(),	u.Y()),
	                     SWIFT_Triple(u.Z(), 0,		-u.X()),
	                     SWIFT_Triple(-u.Y(), u.X(),	0));
	Rb = ustar * R0;
	//Rb.Value[0]=	-u.Z()*R0.Value[3] + u.Y()*R0.Value[6];
	//Rb.Value[3]=	 u.Z()*R0.Value[0] - u.X()*R0.Value[6];
	//Rb.Value[6]=	-u.Y()*R0.Value[0] + u.X()*R0.Value[3];
	//Rb.Value[1]=	-u.Z()*R0.Value[4] + u.Y()*R0.Value[7];
	//Rb.Value[4]=	 u.Z()*R0.Value[1] - u.X()*R0.Value[7];
	//Rb.Value[7]=	-u.Y()*R0.Value[1] + u.X()*R0.Value[4];
	//Rb.Value[2]=	-u.Z()*R0.Value[5] + u.Y()*R0.Value[8];
	//Rb.Value[5]=	 u.Z()*R0.Value[2] - u.X()*R0.Value[8];
	//Rb.Value[8]=	-u.Y()*R0.Value[2] + u.X()*R0.Value[5];
	
	
}

inline void SWIFT_Object::Update_Boxes()
{
	SWIFT_Real vals[9];
	SWIFT_Triple trans_center;
	
	if(cube)
	{
		// To update a cube, simply transform the center and then add and
		// subtract the enlarged radius.
		trans_center = transform * Center_Of_Mass();
		min_bns[0].Set_Value(trans_center.X() - radius);
		min_bns[1].Set_Value(trans_center.Y() - radius);
		min_bns[2].Set_Value(trans_center.Z() - radius);
		max_bns[0].Set_Value(trans_center.X() + radius);
		max_bns[1].Set_Value(trans_center.Y() + radius);
		max_bns[2].Set_Value(trans_center.Z() + radius);
	}
	else
	{
		// To update a dynamic bounding box, find new minimum and maximum
		// vertices and then add and subtract the enlargement factor.
		transform.Rotation().Get_Value(vals);
		min_bns[0].Set_Value(-Mesh()->Root()->
		                     Extremal_Vertex(-SWIFT_Triple(vals), 0, min_es[0]) -
		                     enlargement + transform.Translation().X());
		min_bns[1].Set_Value(-Mesh()->Root()->
		                     Extremal_Vertex(-SWIFT_Triple(vals + 3), 0, min_es[1]) -
		                     enlargement + transform.Translation().Y());
		min_bns[2].Set_Value(-Mesh()->Root()->
		                     Extremal_Vertex(-SWIFT_Triple(vals + 6), 0, min_es[2]) -
		                     enlargement + transform.Translation().Z());
		max_bns[0].Set_Value(Mesh()->Root()->
		                     Extremal_Vertex(SWIFT_Triple(vals), 0, max_es[0]) +
		                     enlargement + transform.Translation().X());
		max_bns[1].Set_Value(Mesh()->Root()->
		                     Extremal_Vertex(SWIFT_Triple(vals + 3), 0, max_es[1]) +
		                     enlargement + transform.Translation().Y());
		max_bns[2].Set_Value(Mesh()->Root()->
		                     Extremal_Vertex(SWIFT_Triple(vals + 6), 0, max_es[2]) +
		                     enlargement + transform.Translation().Z());
	}
}

#endif


