#ifndef __C2A_h
#define __C2A_h

#include <PQP.h>
#include <list>

#include "C2A_Internal.h"
#define articular_used true

enum C2A_Result
{
	OK,						//Succeed
	TOCFound,				//A time of collision is found
	CollisionFound,			//A collision is found
	CollisionFree,			//A collision is found
	CollisionNotFound,		//No collision is found
};

class Transform;
struct C2A_TimeOfContactResult;

//
// Main C2A Algorithm
//
C2A_Result C2A_Solve(Transform *trans00, // begin transformation of obj1
                     Transform *trans01, // end transformation of obj1
                     C2A_Model * obj1_tested,
                     Transform *trans10, // begin transformation of obj2
                     Transform *trans11,  // end transformation of obj2
                     C2A_Model * obj2_tested,
                     Transform & trans0, // contact motion of obj1
                     Transform & trans1, // contact motion of obj2
                     PQP_REAL &time_of_contact,
                     int& number_of_iteration,
                     int & number_of_contact,
                     PQP_REAL th_ca,
                     C2A_TimeOfContactResult& dres);



struct C2A_ContactFeature
{
	int type;//0: vertex; 1: edge; 2: face
	int tri; //index of the triangle
	int fid[3]; //index of the vertex or edge according to the type
	int bid; //the C2A_BV of the triangle 		//huangxin 071126
	
	C2A_ContactFeature(): type(-1) {}
};

struct C2A_ContactPair
{
	C2A_ContactFeature f1;
	C2A_ContactFeature f2;
	double distance; // the diatnce between the two features
	PQP_REAL p1[3]; //the nearest point of f1 to f2
	PQP_REAL p2[3]; //the nearest point of f2 to f1
};

struct C2A_ContactResult
{
	void FreePairsList();
	void SizeTo(int n);
	void Add(const C2A_ContactFeature& f1, const C2A_ContactFeature& f2, PQP_REAL distance, PQP_REAL p[], PQP_REAL q[]);
	
	// statistics
	
	int NumBVTests()
	{
		return num_bv_tests;
	}
	int NumTriTests()
	{
		return num_tri_tests;
	}
	double QueryTimeSecs()
	{
		return query_time_secs;
	}
	
	// The following distance and points established the minimum distance
	// for the models, within the relative and absolute error bounds
	// specified.
	// Points are defined: PQP_REAL p1[3], p2[3];
	C2A_ContactResult();
	~C2A_ContactResult();
	
	
	
	
//	PQP_REAL Distance() { return distance; }
	const PQP_REAL *P1()
	{
		return p1;
	}
	const PQP_REAL *P2()
	{
		return p2;
	}
	
	
public:
	// stats
	
	int num_bv_tests;
	int num_tri_tests;
	double query_time_secs;
	
	// xform from model 1 to model 2
	
	PQP_REAL R[3][3];
	PQP_REAL T[3];
	
	PQP_REAL rel_err;
	PQP_REAL abs_err;
	
//	PQP_REAL distance;
	PQP_REAL p1[3];
	PQP_REAL p2[3];
	int qsize;
	
	//contact features
	C2A_ContactPair *cpairs;
	
	int num_cpairs_alloced;
	int num_cpairs;
	
	bool modelcollided; //whether the two models are colliding
	
	PQP_REAL ctolerance;
};

struct ContactF
{
public:
	void IcI(int a[3], int b[3]);
	void VcV_L(PQP_REAL a[3], PQP_REAL b[3]);
	
	ContactF(
	    int FeatureType_A,
	    int FeatureType_B,
	    int FeatureID_A[3],
	    int FeatureID_B[3],
	    int TriangleID_A,
	    int TriangleID_B,
	    PQP_REAL P_A[3],
	    PQP_REAL P_B[3],
	    PQP_REAL Distance,
	    PQP_REAL Normal[3]);
	    
public:

	int		FeatureType_A;
	int		FeatureType_B;
	
	int		FeatureID_A[3];
	int		FeatureID_B[3];
	
	int		TriangleID_A;
	int		TriangleID_B;
	
	PQP_REAL P_A[3];
	PQP_REAL P_B[3];
	
	PQP_REAL Distance;
	PQP_REAL Normal[3];
	
	friend bool operator > (const ContactF&, const ContactF&);
	friend bool operator < (const ContactF&, const ContactF&);
	
	
};




struct contact_Dist
{

public:
	PQP_REAL Distance;
	PQP_REAL normal[3];
	int BV1_id;
	int BV2_id;
public:

	contact_Dist(PQP_REAL Dis, PQP_REAL N[3], int id1, int id2);
	
};



typedef std::list<ContactF> ContactFList;

typedef std::list<ContactF>::iterator ContactFListIterator;

typedef std::list<contact_Dist> contact_list;
typedef std::list<contact_Dist>::iterator ContactDListIterator;



struct C2A_TimeOfContactResult
{
	// stats
	
	int num_bv_tests;
	int num_tri_tests;
	double query_time_secs;
	
	ContactFList cont_l;
	// xform from model 1 to model 2
	
	PQP_REAL R[3][3];
	PQP_REAL T[3];
	
	PQP_REAL rel_err;
	PQP_REAL abs_err;
	PQP_REAL UpboundTOC;
	
	PQP_REAL distance;
	
	PQP_REAL p1[3];
	PQP_REAL p2[3];
	int qsize;
	PQP_REAL mint; //the minimum time of C2A_BV
	
	
	//toc and collision flag
	PQP_REAL toc;//toc >= 1.0 if collision free, else toc < 1.0
	bool collisionfree;//collision free flag
	int numCA; //the number of CA performed
	PQP_REAL R_toc[3][3], T_toc[3];
	
	int num_contact;
	
	Tri *last_triA;
	Tri *last_triB;
	
	PQP_REAL lastpoint_A[3];
	PQP_REAL lastpoint_B[3];
	
#ifdef articular_used
	
	contact_list cont_dist;
	
	
	PQP_REAL rotAxis[2][3];
	
	PQP_REAL radius2P[2][128];
	
	
	PQP_REAL Velocity[2][128][3];
	
	
	PQP_REAL rotAngle[2][128];
	
	
	int num_link[2];
	
	bool flag_TriangleVisted;
	bool flag_motionbound[2];
	
	
	PQP_REAL accum_RotAngle[2];
	
	PQP_REAL Motion_Upbound[2];
	
	PQP_REAL R_reasolution0[1000][3][3];
	PQP_REAL T_reasolution0[1000][3];
	
	PQP_REAL R_reasolution1[1000][3][3];
	PQP_REAL T_reasolution1[1000][3];
	
	bool flag_intersecstion;
	
#endif
	
public:

	C2A_TimeOfContactResult(): cont_l() {}
	
	int NumBVTests()
	{
		return num_bv_tests;
	}
	
	int NumTriTests()
	{
		return num_tri_tests;
	}
	
	double QueryTimeSecs()
	{
		return query_time_secs;
	}
	
	// The following distance and points established the minimum distance
	// for the models, within the relative and absolute error bounds
	// specified.
	// Points are defined: PQP_REAL p1[3], p2[3];
	
	PQP_REAL Distance()
	{
		return distance;
	}
	
	const PQP_REAL* P1()
	{
		return p1;
	}
	
	const PQP_REAL* P2()
	{
		return p2;
	}
};




PQP_REAL
C2A_TriDistance(PQP_REAL R[3][3],
                PQP_REAL T[3],
                Tri *t1,
                Tri *t2,
                PQP_REAL p[3],
                PQP_REAL q[3],
                C2A_ContactFeature &f1,
                C2A_ContactFeature &f2,
                bool &bCollided);


int
C2A_Distance_At_ContactSpace(C2A_DistanceResult *result,
                             PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
                             PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
                             PQP_REAL upper_bound);

const int C2A_ALL_CONTACTS = 1;  // find all pairwise intersecting triangles
const int C2A_FIRST_CONTACT = 2; // report first intersecting tri pair found

int
C2A_Collide(PQP_CollideResult *result,
            PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
            PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
            int flag = C2A_ALL_CONTACTS);

int
C2A_Distance(C2A_DistanceResult *result,
             PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
             PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
             PQP_REAL rel_err, PQP_REAL abs_err,
             int qsize = 2);

int
C2A_Collide(C2A_DistanceResult *result,
            PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
            PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
            PQP_REAL rel_err, PQP_REAL abs_err,
            int qsize = 2);

class CInterpMotion ;


//CInterpMotion &objmotion1, CInterpMotion &objmotion2,

PQP_REAL C2A_QueryTimeOfContact(CInterpMotion *objmotion1,
                                CInterpMotion *objmotion2,
                                C2A_TimeOfContactResult *res,
                                C2A_Model *o1,
                                C2A_Model *o2,
                                PQP_REAL tolerance_d,
                                PQP_REAL tolerance_t,
                                int qsize = 2);


PQP_REAL C2A_QueryContact(CInterpMotion *objmotion1,
                          CInterpMotion *objmotion2,
                          C2A_TimeOfContactResult *res,
                          C2A_Model *o1,
                          C2A_Model *o2,
                          double threshold);

PQP_REAL
C2A_QueryContactOnly(C2A_TimeOfContactResult *res,
                     PQP_REAL R1[3][3],
                     PQP_REAL T1[3],
                     C2A_Model *o1,
                     PQP_REAL R2[3][3],
                     PQP_REAL T2[3],
                     C2A_Model *o2,
                     double threshold);

int
C2A_TOCStep_articular(C2A_TimeOfContactResult *res,
                      PQP_REAL R1[3][3], PQP_REAL T1[3], C2A_Model *o1,
                      PQP_REAL R2[3][3], PQP_REAL T2[3], C2A_Model *o2,
                      PQP_REAL tolerance, PQP_REAL tolerance_t,
                      int qsize = 2);


double
C2A_TOCStepRecurse_Dis_Deep(C2A_TimeOfContactResult *res,
                            C2A_Model *o1, int b1,
                            C2A_Model *o2, int b2,
                            PQP_REAL mint,
                            PQP_REAL dis_threshold, PQP_REAL m_toc_delta, int Max_Iter, int Max_Deep);



void
Contact_normal_new(C2A_TimeOfContactResult *res, C2A_Model *o1, C2A_Model *o2);

double
CAonNonAdjacentTriangles_New(C2A_TimeOfContactResult *res,
                             C2A_Tri *t1, C2A_Tri *t2,
                             C2A_BV *b1, C2A_BV *b2,
                             PQP_REAL *mint, PQP_REAL *distance, int* Inter,
                             PQP_REAL dis_threshold, PQP_REAL m_toc_delta, int Max_Iter, PQP_REAL* upperbound);


#endif

