//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "model.h"
#include "ModelGraph.h"

//quasirandom
//#include "halton.H"
//#include "hammersley.H"

#ifdef WIN32
#pragma warning(disable:4244)
#endif

//-----------------------------------------------------------------------------
// two global scope models, P, Q
model P,Q;
model& getP(){ return P; }
model& getQ(){ return Q; }
//-----------------------------------------------------------------------------

//
// build the model from a file
// the file will contain an obj model
// once obj file is read, vertices and facets will be built
// other information associated with facets and vertices,
// as well as edges will be build also using model graph
//

bool model::build(const string & name, bool fixed)
{
    if(pts!=NULL) return false; //built

    {//build/copy model
        objReader reader(name);
        if( !reader.Read() )
            return false;
        objModel& data=reader.getModel();

        //allocate memory
        v_size=data.pts.size();
        t_size=data.polys.size();
        vertices=new vertex[v_size];   //
        tris=new triangle[t_size];     //
        assert(vertices&&tris);        //make sure enough memory

        //copy vertices
        for(unsigned int i=0;i<v_size;i++){
            vertices[i].p.set(&data.pts[i].x);
            vertices[i].bk_p=vertices[i].p;  //backup for modification
            vertices[i].n.set(&data.pts[i].nx);
            vertices[i].bk_n=vertices[i].n;  //backup for modification
        }

        //copy triangles
        int tid=0;
        for(list<polygon>::iterator i=data.polys.begin();i!=data.polys.end();i++){
            list<int>& ids=i->pts;
            //check if triangle
            if(ids.size()!=3){
                cerr<<"! Error: polygon "<<tid<<" is not a triangle."<<endl;
                return false;
            }
            int vid=0;
            for(list<int>::iterator j=ids.begin();j!=ids.end();j++){
                tris[tid].v[vid++]=*j;
                vertices[*j].m_f.push_back(tid);
            }
            tid++;
        }
    }//end build/copy model

    {//build model graph and collect informations
        CModelGraph G;
        G.doInit(this);

        //create edges
        e_size=G.getEdgeSize();
        CModelEdge * ptrE=G.getEdges();
        edges=new edge[e_size];
        assert(edges);
        for(unsigned int i=0;i<e_size;i++,ptrE=ptrE->getNext()){
            int v1=edges[i].vid[0]=ptrE->getStartPt();
            int v2=edges[i].vid[1]=ptrE->getEndPt();
            edges[i].fid[0]=ptrE->getFacet0();
            edges[i].fid[1]=ptrE->getFacet1();
            //compute parallel vector
            edges[i].v=edges[i].bk_v=ptrE->getV();
            edges[i].in_n[0]=edges[i].bk_in_n[0]=ptrE->getInNormal(0);
            edges[i].in_n[1]=edges[i].bk_in_n[1]=ptrE->getInNormal(1);
            vertices[v1].m_e.push_back(i);
            vertices[v2].m_e.push_back(i);
        }//end i
        
        //facets
        vector<CModelFacet>& facets=G.getFacets();
        for(unsigned int i=0;i<t_size;i++){
            tris[i].n=tris[i].bk_n=facets[i].n;
        }//end i
    }

    return true;
}

void model::build_collision_detection(bool reverse)
{
    if(cd==NULL){
        //build collision detection model here
        cd=new cd_model();
        assert(cd);
        cd->build(this,reverse);
    }
}

void model::sample(float d, int edgescale)
{
    list<point> samples;

    //sample from edges
    for(unsigned int i=0;i<e_size;i++){

        const triangle& f1=tris[edges[i].fid[0]];
        const triangle& f2=tris[edges[i].fid[1]];
        const vertex&   v1=vertices[edges[i].vid[0]];
        const vertex&   v2=vertices[edges[i].vid[1]];

        //check if edge is significant enough
        //if( f1.n*f2.n>0.9999 )
        if( f1.n*f2.n>0.9999 )
            continue;

        //compute sample size
        Vector3d vec=(v2.p-v1.p);
        float l=vec.norm();
        int size=edgescale*(int)floor(l/d);
        Vector3d delta=vec/size;
        Vector3d enormal=(f1.n+f2.n).normalize();

        //sample
        for(int j=1;j<size;j++){
            point newp;
            newp.p=newp.bk_p=v1.p+delta*j;
            newp.n=newp.bk_n=enormal;
            newp.from='e';
            newp.from_id=i;
            samples.push_back(newp);
        }
    }//end for i

    //sample from vertices
    for(unsigned int i=0;i<v_size;i++){
        point newp;
        newp.p=newp.bk_p=vertices[i].p;
        newp.n=newp.bk_n=vertices[i].n;
        newp.from='v';
        newp.from_id=i;
        samples.push_back(newp);
    }


    //sample from facets

    //a test on quasirandom number generator
    double r[2];
    //halton_dim_num_set(2);
    //hammersley_dim_num_set(2);

    for(unsigned int i=0;i<t_size;i++){
        triangle& t=tris[i];
        Point3d& p1=vertices[t.v[0]].p;
        Point3d& p2=vertices[t.v[1]].p;
        Point3d& p3=vertices[t.v[2]].p;

        //compute triangle area
        Vector3d v1=p2-p1;
        Vector3d v2=p3-p1;
        Vector3d v3=p3-p2;
        Vector3d v4=p1-p3;
        float area=(v1%v2).norm()/2;

        //compute number of points needed
        int size=(int)ceil(area/(d*d));

        //randomly generate these points
        while(size>0)
        {
            //halton(r);
            //hammersley_sequence(1,r);

            r[0]=drand48();
            r[1]=drand48();
            Point3d pos=p1+v1*r[0]+v2*r[1];

            //check if inside the triangle
            if( (v1%(pos-p1))*t.n<=0 ) continue;
            if( (v3%(pos-p2))*t.n<=0 ) continue;
            if( (v4%(pos-p3))*t.n<=0 ) continue;
            //in tri
            point newp;
            newp.p=newp.bk_p=pos;
            newp.n=newp.bk_n=t.n;
            newp.from='f';
            newp.from_id=i;
            samples.push_back(newp);
            size--;
        }//j
    }

    //copy to pts
    p_size=samples.size();
    pts=(point*)malloc(p_size*sizeof(point));
    assert(pts);
    copy(samples.begin(),samples.end(),pts);
    samples.clear();

    //assign to trianles, edges and vertices
    for(unsigned int i=0;i<p_size;i++){
        if( pts[i].from=='v') vertices[pts[i].from_id].sample =i;
        else if( pts[i].from=='e') edges[pts[i].from_id].samples.push_back(i);
        else if( pts[i].from=='f') tris[pts[i].from_id].samples.push_back(i);
    }

}


void model::resample()
{
    for(unsigned int i=0;i<p_size;i++)
    {
        point& pt=pts[i];
        if(pt.from=='v') continue;
        if(pt.from=='e') continue;
        else //pt.from=='f'
        {
            const triangle& t=tris[pt.from_id];
            Point3d& p1=vertices[t.v[0]].p;
            Point3d& p2=vertices[t.v[1]].p;
            Point3d& p3=vertices[t.v[2]].p;

            //compute triangle area
            Vector3d v1=p2-p1;
            Vector3d v2=p3-p1;
            Vector3d v3=p3-p2;
            Vector3d v4=p1-p3;

            //randomly generate these points
            while(true)
            {
                Point3d pos=p1+v1*drand48()+v2*drand48();
                //check if inside the triangle
                float d1=(v1%(pos-p1))*t.n;
                float d2=(v3%(pos-p2))*t.n;
                if(d1*d2<0) continue;
                d2=(v4%(pos-p3))*t.n;
                if(d1*d2<0) continue;
                pt.p=pos;
                break;
            }//while
        }
    }//end i
}


void model::rotate(const Matrix2x2& M)
{
    //Vector2d tmp;
 //   for(unsigned int i=0;i<p_size;i++){
 //       tmp.set(backup[i].p[0],backup[i].p[1]);
 //       tmp=m*tmp;
 //       pts[i].p.set(tmp[0],tmp[1],pts[i].p[2]);
 //       tmp.set(backup[i].n[0],backup[i].n[1]);
 //       tmp=m*tmp;
 //       pts[i].n.set(tmp[0],tmp[1],0);
 //   }
    Vector2d tmp;
    //rotate points
    for(unsigned int i=0;i<p_size;i++){
        tmp.set(pts[i].bk_p[0],pts[i].bk_p[1]);
        tmp=M*tmp;
        pts[i].p.set(tmp[0],tmp[1]);
        tmp.set(pts[i].bk_n[0],pts[i].bk_n[1]);
        tmp=M*tmp;
        pts[i].n.set(tmp[0],tmp[1]);
    }

    //rotate vertices
    for(unsigned int i=0;i<v_size;i++){
        tmp.set(vertices[i].bk_p[0],vertices[i].bk_p[1]);
        tmp=M*tmp;
        vertices[i].p.set(tmp[0],tmp[1]);
        tmp.set(vertices[i].bk_n[0],vertices[i].bk_n[1]);
        tmp=M*tmp;
        vertices[i].n.set(tmp[0],tmp[1]);
    }
    
    //rotate edges
    for(unsigned int i=0;i<e_size;i++){
        for(int j=0;j<2;j++){
            tmp.set(edges[i].bk_in_n[j][0],edges[i].bk_in_n[j][1]);
            tmp=M*tmp;
            edges[i].in_n[j].set(tmp[0],tmp[1]);
        }
        
        tmp.set(edges[i].v[0],edges[i].v[1]);
        tmp=M*tmp;
        edges[i].v.set(tmp[0],tmp[1]);
    }
    //rotate facets
    for(unsigned int i=0;i<t_size;i++){
        tmp.set(tris[i].n[0],tris[i].n[1]);
        tmp=M*tmp;
        tris[i].n.set(tmp[0],tmp[1]);
    }
}




void model::rotate(const Matrix3x3& M)
{
    Vector3d tmp;
    //rotate points
    for(unsigned int i=0;i<p_size;i++){
        tmp.set(pts[i].bk_p.get());
        pts[i].p=M*tmp;
        pts[i].n=(M*pts[i].bk_n).normalize();
    }

    //rotate vertices
    for(unsigned int i=0;i<v_size;i++){
        tmp.set(vertices[i].bk_p.get());
        vertices[i].p=M*tmp;
        vertices[i].n=(M*vertices[i].bk_n).normalize();
    }
    
    //rotate edges
    for(unsigned int i=0;i<e_size;i++){
        edges[i].in_n[0]=(M*edges[i].bk_in_n[0]).normalize();
        edges[i].in_n[1]=(M*edges[i].bk_in_n[1]).normalize();
        edges[i].v=(M*edges[i].bk_v).normalize();
    }
    //rotate facets
    for(unsigned int i=0;i<t_size;i++)
        tris[i].n=(M*tris[i].bk_n).normalize();
}


void model::negate()
{
    for(unsigned int i=0;i<p_size;i++){
        Point3d& pt=pts[i].p;
        pt.set(-pt[0],-pt[1],-pt[2]);

        Point3d& pt2=pts[i].bk_p;
        pt2.set(-pt2[0],-pt2[1],-pt2[2]);

        pts[i].n=-pts[i].n;
        pts[i].bk_n=-pts[i].bk_n;
    }

    for(unsigned int i=0;i<v_size;i++){
        Point3d& pt=vertices[i].p;
        pt.set(-pt[0],-pt[1],-pt[2]);

        Point3d& pt2=vertices[i].bk_p;
        pt2.set(-pt2[0],-pt2[1],-pt2[2]);

        vertices[i].n=-vertices[i].n;
        vertices[i].bk_n=-vertices[i].bk_n;
    }

    for(unsigned int i=0;i<t_size;i++){
        tris[i].n=-tris[i].n;
        tris[i].bk_n=-tris[i].bk_n;
    }

    for(unsigned int i=0;i<e_size;i++){
        edge& e=edges[i];
        e.v=-e.v;
        e.in_n[0]=-e.in_n[0];
        e.in_n[1]=-e.in_n[1];
        e.bk_v=-e.bk_v;
        e.bk_in_n[0]=-e.bk_in_n[0];
        e.bk_in_n[1]=-e.bk_in_n[1];
    }
}
