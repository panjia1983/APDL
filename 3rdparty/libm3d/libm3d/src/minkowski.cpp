//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

//compute minkowski sum from two set of points
#include "minkowski.h"

namespace libm3d 
{

extern double current_rot[3][3];

//build minkowski points from a vertex from P and a facet from Q
void mksum::buildVF(model * P, model * Q, bool pfirst)
{
    mksum_pt pt;
    typedef std::list<unsigned int>::iterator IT;
    //P.v and Q.facets
    for(unsigned int p=0;p<P->v_size;p++){
        vertex& vp=P->vertices[p];
        if(pfirst) pt.pid=vp.sample;
        else pt.qid=vp.sample;

        for(unsigned int q=0;q<Q->t_size;q++){
            triangle& tq=Q->tris[q];
            if( mksum_test(P,vp,tq) ){ 
                for(IT s=tq.samples.begin();s!=tq.samples.end();s++){
                    if(pfirst) pt.qid=*s;
                    else pt.pid=*s;
                    pts.push_back(pt);
                }//end for 
            }
        }//end q
    }//end p
}

void mksum::buildVV(model * P, model * Q)
{
    mksum_pt pt;
    for(unsigned int p=0;p<P->v_size;p++){

        vertex& vp=P->vertices[p];
        pt.pid=vp.sample;

        for(unsigned int q=0;q<Q->v_size;q++){
            vertex& vq=Q->vertices[q];
            if( mksum_test(P,Q,vp,vq) ){
                pt.qid=vq.sample;
                pts.push_back(pt);
            }
        }//end q

    }//end p
}

void mksum::buildVE(model* P, model * Q, bool pfirst)
{
    mksum_pt pt;
    typedef std::list<unsigned int>::iterator IT;
    for(unsigned int p=0;p<P->v_size;p++){
        vertex& vp=P->vertices[p];
        if(pfirst) pt.pid=vp.sample;
        else pt.qid=vp.sample;
        for(unsigned int q=0;q<Q->e_size;q++){
            edge& eq=Q->edges[q];
            if(eq.samples.empty()) continue;
            if( mksum_test(P,Q,vp,eq) ){
                for(IT i=eq.samples.begin();i!=eq.samples.end();i++){
                    if(pfirst) pt.qid=*i;
                    else pt.pid=*i;
                    pts.push_back(pt);
                }
            }
        }//end q
    }//end p
}

//build minkowski points from a pair of edges!
void mksum::buildEE(model * P, model * Q)
{
    mksum_pt pt;
    typedef std::list<unsigned int>::iterator IT;
    //P.e and Q.e
    for(unsigned int p=0;p<P->e_size;p++){
        edge & ep=P->edges[p];
        if(ep.samples.empty()) continue; //nothing to generate
        for(unsigned int q=0;q<Q->e_size;q++){
            edge & eq=Q->edges[q];
            if(eq.samples.empty()) continue; //nothing to generate
            if( mksum_test(ep,eq) ){ 
                for(IT i=ep.samples.begin();i!=ep.samples.end();i++){
                    pt.pid=*i;
                    for(IT j=eq.samples.begin();j!=eq.samples.end();j++){
                        pt.qid=*j;
                        pts.push_back(pt);
                    }//end for j
                }//end for i
            }
        }//end q
    }//end p
}

//brute force, bulld all point pairs
void mksum::buildAll(model * P, model * Q)
{
    mksum_pt pt;
    for(unsigned int q=0;q<Q->p_size;q++)
    {
        for(unsigned int p=0;p<P->p_size;p++){

            pt.pid=(unsigned int)p;
            pt.qid=(unsigned int)q;
            pts.push_back(pt);

        }//q
    }//p
}

//test if a vertex and a facet can be a potiential candidate
bool mksum::mksum_test(model * m, vertex& v, triangle& t)
{
    //return true;
    typedef std::list<unsigned int>::iterator  IT;
    for(IT i=v.m_e.begin();i!=v.m_e.end();i++){
        //unsigned int& id=*i;
        edge& e=m->edges[*i];
        bool pos=true;
        if( &(m->vertices[e.vid[0]])==&v )
            pos=e.v*t.n>0;
        else 
            pos=e.v*t.n<0;

        if(pos) return false; //cannot be a candidate
    }//end for
    return true; //yes it is a cadidate
}


bool mksum::mksum_test(model * P, model * Q, vertex& p, vertex& q)
{
    typedef std::list<unsigned int>::iterator  IT;

    //for each facet incident to q
    for(IT i=q.m_f.begin();i!=q.m_f.end();i++){
        if(mksum_test(P,p,Q->tris[*i])) return true;
    }

    //for each facet incident to p
    for(IT i=p.m_f.begin();i!=p.m_f.end();i++){
        if(mksum_test(Q,q,P->tris[*i])) return true;
    }

    for(IT i=p.m_e.begin();i!=p.m_e.end();i++){
        edge& ep=P->edges[*i];
        for(IT j=q.m_e.begin();j!=q.m_e.end();j++){
            edge& eq=Q->edges[*j];
            if( mksum_test(ep,eq) ) return true;
        }
    }

    return false;
}

bool mksum::mksum_test(model * P, model * Q, vertex& vp, edge& eq)
{
    if(mksum_test(P,vp,Q->tris[eq.fid[0]])) return true;
    if( eq.fid[1]>=0 )
        if(mksum_test(P,vp,Q->tris[eq.fid[1]])) return true;

    typedef std::list<unsigned int>::iterator  IT;
    for(IT i=vp.m_e.begin();i!=vp.m_e.end();i++)
        if( mksum_test(P->edges[*i],eq) ) return true;
    return false;
}

//test if a pair of edges can be a potiential candidate
bool mksum::mksum_test(edge& e1, edge& e2)
{
    mathtool::Vector3d n=e1.v%e2.v;
    if(n.normsqr()==0) return false; //parallel edges

    float d1=n*e1.in_n[0];
    float d2=n*e1.in_n[1];
    if(d1*d2<0) return false;

    d2=n*e2.in_n[0];
    if(d1*d2<0) return false;

    d2=n*e2.in_n[1];
    if(d1*d2<0) return false;

    return true; //yes it is a cadidate
}


//-----------------------------------------------------------------------------
//
//
// Compute Normals
//
//
//-----------------------------------------------------------------------------

void mksum::mksum_normal
(model * P, model * Q, vertex& p, vertex& q, mathtool::Vector3d& n)
{
    n=(p.n+q.n).normalize();
}

void mksum::mksum_normal_ve
(model * P, model * Q, const point& p, const point& q, mathtool::Vector3d& n)
{
    vertex& v=P->vertices[p.from_id];
    edge&   e=Q->edges[q.from_id];

    std::vector<mathtool::Vector3d> candidates;

    if(mksum_test(P,v,Q->tris[e.fid[0]]))
        candidates.push_back(Q->tris[e.fid[0]].n);

    if( e.fid[1]>=0 ){
        if(mksum_test(P,v,Q->tris[e.fid[1]]))
            candidates.push_back(Q->tris[e.fid[1]].n);
    }

    if(candidates.empty()){
        typedef std::list<unsigned int>::iterator  IT;
        for(IT i=v.m_e.begin();i!=v.m_e.end();i++){
            if( mksum_test(P->edges[*i],e) ){
                mathtool::Vector3d tmp;
                mksum_normal(P->edges[*i],e,tmp);
                candidates.push_back(tmp);
            }
        }
    }

    n.set(0,0,0);
    int size=candidates.size();

    if(size>1){
        bool found=false;
        for(int i=0;i<size;i++)
        {
            double pos[3];
            for(int k=0;k<3;k++)
                pos[k]=p.p[k]+q.p[k]+candidates[i][k]*1e-5;
            bool r=is_in_collision(this->P->cd,this->Q->cd,pos,current_rot);
            if(!r){
                n=n+candidates[i];
                found=true;
            }
        }
        if(found)
            n=n.normalize();
    }
    else{
        if(candidates.empty())
            n=(p.n+q.n).normalize();
        else
            n=candidates.front();
    }
}

void mksum::mksum_normal(edge& e1, edge& e2, mathtool::Vector3d& n)
{
    n=e1.v%e2.v;
    float d=n*e1.in_n[0];
    if(d!=0){ if(d>0) n=-n; }
    else{
        d=n*e1.in_n[1];
        if(d!=0){ if(d>0) n=-n; }
        else{
            d=n*e2.in_n[0];
            if(d!=0){ if(d>0) n=-n; }
            else{
                d=n*e2.in_n[1];
                if(d!=0){ if(d>0) n=-n; }
            }
        }
    }

    n=n.normalize();
}


}
