/*************************************************************************\

  Copyright 2009	Computer Graphics Lab,
									Ewha Womans University,
									Seoul, Korea.
  All Rights Reserved.

  The authors may be contacted via:

  Mailing Address:     Xinyu Zhang
                       Min Tang
                       Department of Computer Science & Engineering
                       Ewha Womans University
                       Seoul, 120-750, Korea.

  EMail:               zhangxy@ewha.ac.kr
                       tangmin@ewha.ac.kr


\**************************************************************************/

#ifndef __C2A_Internal_h
#define __C2A_Internal_h

#include "C2A_Tri.h"
#include "C2A_BV.h"
#include <vector>


class PQP_Model;
class Coord3D;
class C2A_Model: public PQP_Model
{
public:

	C2A_Tri *trisConst;
	
	C2A_BV  *B;
	
	int maxdeep;
	PQP_REAL P[3];
	PQP_REAL Q[3];
	
	int **VertexFace;
	int *number_vertexFace;
	
	
	PQP_REAL **Vertex;
	
	int num_vertex;
	
	//use a single tri coherence
	bool bMotionCoherence;
	int motionTri;
	int motionBV;
	
	//use a local tree coherence, set local C2A_BV level typicall 3
	int level; //set local C2A_BV level typicall 3
	int levelPID;//the parent ID at level tree from the leaf
	// int numLocalTris; // pow(2, level)
	// int *localTris;//all the tris in the local C2A_BV tree
	
	C2A_Model();
	~C2A_Model();
	BV *child(int n)
	{
		return (BV*) & ((C2A_BV*)b)[n];
	}
	
	C2A_BV *getC2Achild(int n)
	{
		return &B[n];
	}
	
	
	int BeginModel(int n = 8);
	
	int AddTri(const PQP_REAL *p1, const PQP_REAL *p2, const PQP_REAL *p3,
	           int i, int i1, int i2, int i3);
	int AddTri(const PQP_REAL *p1, const PQP_REAL *p2, const PQP_REAL *p3,
	           int id);
	int EndModel();
	int MemUsage(int msg);  // returns model mem usage.
	// prints message to stderr if msg == TRUE
	
	PQP_REAL com[3]; // the center of Mass, huangxin
	PQP_REAL radius; // the radius of the model using com, huangxin
	
	std::vector<Coord3D> max_clear_conf_;
	
	
	PQP_REAL Max_Cross_Product(PQP_REAL dir[3]);
	void SetCenterOfMass(PQP_REAL c[3]);
	void ComputeCenterOfMass();
	void ComputeRadius();
	
	inline C2A_Tri* GetTriangle(int idx)
	{
		return &((C2A_Tri*)tris)[idx];
	}
	
};


// Add by Liangjun: to compute the aboslute transform
void ComputeAbsTra(PQP_REAL R[3][3], PQP_REAL T[3], C2A_Model *o, int bv_ind);



struct C2A_DistanceResult: public PQP_DistanceResult
{
	int t1;//the closest pair of t1
	int t2;//the closest pair of t2
};

#endif