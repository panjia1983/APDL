//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_CD_H_
#define _BF_CD_H_

#include <stdlib.h>
#include <RAPID.H>


struct model;

struct cd_model
{
    cd_model(){
        m_cd=NULL;
    }

    ~cd_model(){
        delete m_cd;
        m_cd=NULL;
    }

    void build(model * m, bool reverse);
    RAPID_model * m_cd;
};

//check for collision

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], double R[3][3]);

bool is_in_collision
(model* P, model * Q, double p[3], float radian=0);

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], float radian=0);

bool is_in_collision
(model* P, model * Q, double p[3], const float * r);


bool is_in_collision
(model* P, model * Q, double p[3], double R[3][3]);

bool is_in_collision
(cd_model* P, cd_model * Q, double p[3], const float * r);


#endif //_BF_CD_H_


