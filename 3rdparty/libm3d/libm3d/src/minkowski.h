//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MINKOWSKI_H_
#define _BF_MINKOWSKI_H_

//define a small offset for minkowski point
#define MK_POINT_OFFSET 1e-5f

//compute minkowski sum from two set of points

#include "model.h"
#include <vector.h>

namespace libm3d 
{

struct mksum_pt
{
    mksum_pt(){ pid=qid=0; }
    unsigned int pid, qid;
};

typedef std::list<mksum_pt> MKPTS;
typedef MKPTS::iterator MKPTIT;


//
// This class maintains a list of points
// that are the minkowski sum of two point sets
// sampled from the boundaries of two models
// 

class mksum
{
public:

    //Constructor and Destructor
    mksum(){
        P=Q=NULL;
    }

    ~mksum(){ destroy(); }

    void destroy(){
        pts.clear();
        pt_size=0;
    }

    //An efficient way to generate points are
    //are likely to be points on the Minkowski sum  of
    //two models P&Q
    void build_with_filter(model * P, model * Q)
    {
        if(P==NULL || Q==NULL) //bad input
            return;

        this->P=P;
        this->Q=Q;

        //P.v and Q.facets
        buildVF(P,Q);

        //P.v and Q.facets
        buildVF(Q,P,false);

        //P.v and Q.v
        buildVV(P,Q);

        //P.v and Q.e
        buildVE(P,Q);

        //Q.v and P.e
        buildVE(Q,P,false);

        //P.e and Q.e
        buildEE(P,Q);

        pt_size=pts.size();
    }

    
    // a brute force way to generate minkowski sum points
    // this is mainly for experiments
    void build_without_filter(model * P, model * Q)
    {
        if(P==NULL || Q==NULL) //bad input
            return;

        this->P=P;
        this->Q=Q;
        buildAll(P,Q);
        pt_size=pts.size();
    }


    //----------------------------------------------------------------------------
    // 
    // Access functions
    //
    //----------------------------------------------------------------------------

    //Position of a minkowski sum point
    template<class T> void getMPos
    (const mksum_pt& mpt, T pos[3])
    {
        getMPos(mpt.pid,mpt.qid,pos);
    }

    //Offset position of a minkowski sum point
    template<class T> void getMPos_offset
    (const mksum_pt& mpt, T pos[3])
    {
        getMPos_offset(mpt.pid,mpt.qid,pos);
    }

    //Normal direction of a minkowski sum point
    template<class T>
    void getMNormal(const mksum_pt& mpt, T n[3])
    {
        getMNormal(mpt.pid,mpt.qid,n);
    }

private:

    //build minkowski points from a vertex from P and a facet from Q
    void buildVF(model * P, model * Q, bool pfirst=true);

    //build minkowski points from a pair of vertices
    void buildVV(model * P, model * Q);

    //build minkowski points from a vertex and an edge
    void buildVE(model* P, model * Q, bool pfirst=true);

    //build minkowski points from a pair of edges!
    void buildEE(model * P, model * Q);

    //brute force method, for test only
    void buildAll(model * P, model * Q);

    //test if a vertex and a facet can be a potiential candidate
    bool mksum_test(model * m, vertex& v, triangle& t);

    //test if a pair of vertices can be a potiential candidate
    bool mksum_test(model * P, model * Q, vertex& p, vertex& q);

    //test if a vertex and an edge can be a potiential candidate
    bool mksum_test(model * P, model * Q, vertex& vp, edge& eq);   

    //test if a pair of edges can be a potiential candidate
    bool mksum_test(edge& e1, edge& e2);

    //minkowski sum position without offset
    template<class T> void getMPos
    (const unsigned int& pid, const unsigned int& qid, T pos[3])
    {
        const point& p=P->pts[pid];
        const point& q=Q->pts[qid];
        for(int i=0;i<3;i++)
            pos[i]=q.p[i]+p.p[i];
    }

    //minkowski sum position with a tiny offset
    template<class T> void getMPos_offset
    (const unsigned int& pid, const unsigned int& qid, T pos[3])
    {
        const point& p=P->pts[pid];
        const point& q=Q->pts[qid];
        for(int i=0;i<3;i++)
            pos[i]=(q.p[i]+p.p[i])+MK_POINT_OFFSET*q.n[i];
    }

    // mksum point normal
    // Note: Not good enough and in fact slow computation of point normal
    //       need a better and more efficient method!!
    void mksum_normal
    (model * P, model * Q, vertex& p, vertex& q, mathtool::Vector3d& n);
    
    void mksum_normal_ve
    (model * P, model * Q, const point& p, const point& q, mathtool::Vector3d& n);

    void mksum_normal(edge& e1, edge& e2, mathtool::Vector3d& n);

    //functions that compute the position and the normal of M+ points
    template<class T> void getMNormal
    (const unsigned int& pid, const unsigned int& qid, T n[3]){
        const point& p=P->pts[pid];
        const point& q=Q->pts[qid];

        mathtool::Vector3d tmp;

        if(p.from=='f') 
        {
            for(int i=0;i<3;i++) n[i]=P->tris[p.from_id].n[i];
        }
        else if(q.from=='f') 
        {
            for(int i=0;i<3;i++) n[i]=Q->tris[q.from_id].n[i];
        }
        else if(p.from=='e' && q.from=='e')
        {
            mksum_normal(P->edges[p.from_id],Q->edges[q.from_id],tmp);
            for(int i=0;i<3;i++) n[i]=tmp[i];
        }
        else if(p.from=='v' && q.from=='v')
        {
            mksum_normal(P,Q,P->vertices[p.from_id],Q->vertices[q.from_id],tmp);
            for(int i=0;i<3;i++) n[i]=tmp[i];
        }
        else if(p.from=='v' && q.from=='e')
        {
            mksum_normal_ve(P,Q,p,q,tmp);
            for(int i=0;i<3;i++) n[i]=tmp[i];
        }
        else if(p.from=='e' && q.from=='v')
        {
            mksum_normal_ve(Q,P,q,p,tmp);
            for(int i=0;i<3;i++) n[i]=tmp[i];
        }
        else{
            cerr<<"! Error: getMNormal"<<endl;
        }
    }


public:

    model * P, * Q;   //source of the minkowski sum
    MKPTS pts;        //sum
    unsigned int pt_size;
};

//-----------------------------------------------------------------------------

}
#endif //_BF_MINKOWSKI_H_

