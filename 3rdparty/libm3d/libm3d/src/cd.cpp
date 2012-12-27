//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "model.h"
#include "cd.h"
#include <Matrix.h>
using namespace mathtool;


void cd_model::build(model * m, bool reverse)
{
    if(m==NULL || m_cd!=NULL)
        return;

    m_cd = new RAPID_model;

    m_cd->BeginModel();
    for( unsigned int t=0;t<m->t_size;t++ ){
        const triangle & tri=m->tris[t];
        const Point3d& pt1=m->vertices[tri.v[0]].p;
        const Point3d& pt2=m->vertices[tri.v[1]].p;
        const Point3d& pt3=m->vertices[tri.v[2]].p;
        if(reverse){
            double p1[3]={-pt1[0],-pt1[1],-pt1[2]};
            double p2[3]={-pt2[0],-pt2[1],-pt2[2]};
            double p3[3]={-pt3[0],-pt3[1],-pt3[2]};
            m_cd->AddTri(p1,p2,p3,t);
        }
        else{
            double p1[3]={pt1[0],pt1[1],pt1[2]};
            double p2[3]={pt2[0],pt2[1],pt2[2]};
            double p3[3]={pt3[0],pt3[1],pt3[2]};
            m_cd->AddTri(p1,p2,p3,t);
        }
    }
    //end RAPID model
    m_cd->EndModel();
}

bool is_in_collision
(model* P, model * Q, double p[3], float radian)
{
    return is_in_collision(P->cd,Q->cd,p,radian);
}

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], float radian)
{
    double c_r=cos(radian);
    double s_r=sin(radian);
    double R[3][3]={ {c_r,-s_r,0} , {s_r,c_r,0} , {0,0,1} };    
    return is_in_collision(P,Q,p,R);
}

bool is_in_collision
(model* P, model * Q, double p[3], const float * r)
{
    return is_in_collision(P->cd,Q->cd,p,r);
}

bool is_in_collision
(model* P, model * Q, double p[3], double R[3][3])
{
    return is_in_collision(P->cd,Q->cd,p,R);
}

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], const float * r)
{
    double R[3][3]={ {1,0,0} , {0,1,0} , {0,0,1} };    
    
    Matrix3x3 M; //build matrix
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
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            R[i][j]=M[i][j];

    return is_in_collision(P,Q,p,R);
}

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], double R[3][3])
{
    static double I[3][3]={ {1,0,0} , {0,1,0} , {0,0,1} };
    static double O[3]={0,0,0};
    RAPIDres cres;
    RAPID_Collide(cres,R, p, P->m_cd, I, O, Q->m_cd, RAPID_FIRST_CONTACT);
    return cres.RAPID_num_contacts>0;
}



