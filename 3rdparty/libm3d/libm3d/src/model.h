//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_MODEL_H_
#define _BF_MODEL_H_

#include <Point.h>
#include <Vector.h>
#include <Matrix.h>
using namespace mathtool;

#include <string>
#include <cassert>
using namespace std;

#include "objReader.h"
#include "cd.h"


//these two functions are defined in model.cpp
struct model; 
model& getP();
model& getQ();


//this is a sampled point from the model
struct point
{
    Point3d p;  //position
    Vector3d n; //normal
    char from;  //source, v, f, or e
    int from_id;
    
    //backups
    Point3d bk_p; 
    Vector3d bk_n;
};

//a triangle of the model
struct triangle
{
    unsigned int v[3]; //vertex id
    Vector3d n;
    list<unsigned int> samples;

    //backups
    Vector3d bk_n;
};

//a vertex of the model
struct vertex
{
    Point3d p;  //position
    Vector3d n; //normal
    list<unsigned int> m_f;
    list<unsigned int> m_e; //a list of edges
    unsigned int sample;

    //backups
    Point3d bk_p; 
    Vector3d bk_n;
};

//an edge of the model
struct edge
{
    int vid[2];
    int fid[2];
    Vector3d v;       //parallel vertor
    Vector3d in_n[2]; //inface normals
    list<unsigned int> samples;

    //backups
    Vector3d bk_v;       //parallel vertor
    Vector3d bk_in_n[2]; //inface normals
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

    // build collision detection model from this model
    // when reverse is true, a reverse version of the model
    // will be built
    void build_collision_detection(bool reverse);

    //sample points
    void sample(float d, int edgescale=1);
    void resample();

    //rotate points
    void rotate(const Matrix2x2& m);
    void rotate(const Matrix3x3& M);

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


inline void computeRotationMatrix(float r, Matrix2x2& m)
{
    float c_r=cos(r);
    float s_r=sin(r);
    m.set(c_r,-s_r,s_r,c_r);
}


inline void computeRotationMatrix(const Point3d& r, Matrix3x3& M)
{
    float c_r=cos(r[0]);
    float s_r=sin(r[0]);
    Matrix3x3 mx(1,0,0, 0, c_r,-s_r, 0, s_r,c_r);
    c_r=cos(r[1]);
    s_r=sin(r[1]);
    Matrix3x3 my(c_r, 0, s_r, 0, 1, 0, -s_r, 0, c_r);
    c_r=cos(r[2]);
    s_r=sin(r[2]);
    Matrix3x3 mz(c_r,-s_r, 0, s_r,c_r,0,0,0,1);
    M=mz*my*mx;
}

#endif //_BF_MODEL_H_
