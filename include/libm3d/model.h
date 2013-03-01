//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <Point.h>
#include <Vector.h>
#include <Matrix.h>

#include <string>
#include <cassert>

#include "objReader.h"
#include "cd.h"

namespace libm3d 
{

//these two functions are defined in model.cpp
struct model; 
model& getP();
model& getQ();


//this is a sampled point from the model
struct point
{
    mathtool::Point3d p;  //position
    mathtool::Vector3d n; //normal
    char from;  //source, v, f, or e
    int from_id;
    
    //backups
    mathtool::Point3d bk_p; 
    mathtool::Vector3d bk_n;
};

//a triangle of the model
struct triangle
{
    unsigned int v[3]; //vertex id
    mathtool::Vector3d n;
    std::list<unsigned int> samples;

    //backups
    mathtool::Vector3d bk_n;
};

//a vertex of the model
struct vertex
{
    mathtool::Point3d p;  //position
    mathtool::Vector3d n; //normal
    std::list<unsigned int> m_f;
    std::list<unsigned int> m_e; //a list of edges
    unsigned int sample;

    //backups
    mathtool::Point3d bk_p; 
    mathtool::Vector3d bk_n;
};

//an edge of the model
struct edge
{
    int vid[2];
    int fid[2];
    mathtool::Vector3d v;       //parallel vertor
    mathtool::Vector3d in_n[2]; //inface normals
    std::list<unsigned int> samples;

    //backups
    mathtool::Vector3d bk_v;       //parallel vertor
    mathtool::Vector3d bk_in_n[2]; //inface normals
};

struct model
{
    //initialization
    model(){ 
        pts=NULL; 
        vertices=NULL;
        edges=NULL;
        tris=NULL; 
        cd=NULL; 
    }

    ~model(){ destroy(); }

    void destroy(){
        delete [] pts;      pts=NULL;
        delete [] tris;     tris=NULL;
        delete [] edges;    edges=NULL;
        delete [] vertices; vertices=NULL;
        delete cd;          cd=NULL;
    }

    //build this model
    bool build(const string & name, bool fixed);

	bool build(const std::vector<std::vector<double> >& points, const std::vector<std::vector<int> >& tri_indices);

    // build collision detection model from this model
    // when reverse is true, a reverse version of the model
    // will be built
    void build_collision_detection(bool reverse);

    //sample points
    void sample(float d, int edgescale=1);
    void resample();

    //rotate points
    void rotate(const mathtool::Matrix2x2& m);
    void rotate(const mathtool::Matrix3x3& M);

    //negate point/facets ...
    void negate();

    //simply, convert points to a float array
    //caller is responsible to release the allocated memory
    //by this function
    float * pts2float(){
        float * tmp=new float[p_size*6];
        assert(tmp);
        for(unsigned int i=0;i<p_size;i++){
            int id=i*6;
            tmp[id]  =pts[i].p[0];
            tmp[id+1]=pts[i].p[1];
            tmp[id+2]=pts[i].p[2];
            tmp[id+3]=pts[i].n[0];
            tmp[id+4]=pts[i].n[1];
            tmp[id+5]=pts[i].n[2];
        }
        return tmp;
    }


    //data
    point    * pts;       //sampled points
    vertex   * vertices;  //vertices
    triangle * tris;      //triangles
    edge     * edges;     //edges
    cd_model * cd;        //collision detection model
    unsigned int p_size, v_size, e_size, t_size;
};


inline void computeRotationMatrix(float r, mathtool::Matrix2x2& m)
{
    float c_r=cos(r);
    float s_r=sin(r);
    m.set(c_r,-s_r,s_r,c_r);
}


inline void computeRotationMatrix(const mathtool::Point3d& r, mathtool::Matrix3x3& M)
{
    float c_r=cos(r[0]);
    float s_r=sin(r[0]);
    mathtool::Matrix3x3 mx(1,0,0, 0, c_r,-s_r, 0, s_r,c_r);
    c_r=cos(r[1]);
    s_r=sin(r[1]);
    mathtool::Matrix3x3 my(c_r, 0, s_r, 0, 1, 0, -s_r, 0, c_r);
    c_r=cos(r[2]);
    s_r=sin(r[2]);
    mathtool::Matrix3x3 mz(c_r,-s_r, 0, s_r,c_r,0,0,0,1);
    M=mz*my*mx;
}

}

#endif //_BF_MODEL_H_
