//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _BF_OCT_H_
#define _BF_OCT_H_

#include "minkowski.h"


//
//
// This file contains the definition of octree class
//
//

namespace libm3d 
{

template<int KS=8>
struct _octree_cell
{
    _octree_cell(mksum * mk)
    {
        this->mk=mk;
        parent=NULL;
        for(int i=0;i<KS;i++) kids[i]=NULL;
        visited=false;
    }

    void destroy()
    {
        if(kids[0]==NULL) return;
        for(int i=0;i<KS;i++){
            kids[i]->destroy();
            delete kids[i];
        }
    }

    //this is used for initial contruction of the octree
    void build(float bbox[6])
    {   
        //check if boundary point list is empty
        //or if not touching the external any more
        if(bd_points.empty() && inside(bbox))
            return;

        //check total size
        int total=points.size();
        if(total<KS*4)
            return;

        subdivide(true);
        if(kids[0]==NULL) return; //done!
        for(int i=0;i<KS;i++) 
            kids[i]->build(bbox);
    }

    //this is used during neighboring cell query
    void refine(float query[6])
    {   
        int total=points.size();
        if(total<KS*2)
            return;

        //not touching the external any more
        if( !intersect(query) )
            return;

        subdivide();
        if(kids[0]==NULL) return; //done!
        for(int i=0;i<KS;i++) 
            kids[i]->refine(query);
    }

    int index(float * pos, const mathtool::Point3d& mid) const 
    {
        int x=(pos[0]>mid[0])?1:0;
        int y=(pos[1]>mid[1])?1:0;
        int z=(pos[2]>mid[2])?1:0;
        return x+y*2+z*4;
    }

    void subdivide(bool split_bd_points=false){

        MKPTS subs[KS];
        MKPTS sub_pts[KS];

        //compute mid
        mathtool::Point3d mid((box[0]+box[1])/2, (box[2]+box[3])/2, (box[4]+box[5])/2);

        //split accoding to mid
        float pos[3];
        for(MKPTIT i=points.begin();i!=points.end();i++){
            mk->getMPos(*i,pos);
            subs[index(pos,mid)].push_back(*i);
        }
        points.clear(); //not used, clear to save memory

        //split bd points
        if(split_bd_points){
            for(MKPTIT i=bd_points.begin();i!=bd_points.end();i++){
                mk->getMPos(*i,pos);
                sub_pts[index(pos,mid)].push_back(*i);
            }
            bd_points.clear();
        }

        //create kids!
        for(int i=0;i<KS;i++)
        {
            kids[i]=new _octree_cell<KS>(mk);
            assert(kids[i]);
            int x=(i%2<1)?0:1;
            int y=(i%4<2)?0:1;
            int z=(i<4)?0:1;
            kids[i]->box[0]=(x==1)?mid[0]:box[0];
            kids[i]->box[1]=(x==1)?box[1]:mid[0];
            kids[i]->box[2]=(y==1)?mid[1]:box[2];
            kids[i]->box[3]=(y==1)?box[3]:mid[1];
            kids[i]->box[4]=(z==1)?mid[2]:box[4];
            kids[i]->box[5]=(z==1)?box[5]:mid[2];
            kids[i]->parent=this;
            kids[i]->points.swap(subs[i]);
            if(split_bd_points)
                kids[i]->bd_points.swap(sub_pts[i]);
        }//end for i
    }

    //find intersecting cells with the given box
    void getIntersection
    (float query[6], std::list<_octree_cell<KS> *>& result, _octree_cell<KS>* from)
    {
        for(int i=0;i<KS;i++){
            if(kids[i]==from || kids[i]->visited) continue;
            if(kids[i]->intersect(query)){
                //check if refinement is needed
                //this will create kids when the cell is too big
                if(kids[i]->kids[0]==NULL) 
                    kids[i]->refine(query);

                //add or go deep
                if(kids[i]->kids[0]==NULL){
                    kids[i]->visited=true;
                    result.push_back(kids[i]); //kids is a leaf
                }
                else
                    kids[i]->getIntersection(query,result,this);
            }
        }
        if(from!=parent && parent!=NULL)
            parent->getIntersection(query,result,this);
    }

    bool intersect( float query[6] )
    {
        if(fabs(query[0]+query[1]-box[0]-box[1])>(query[1]-query[0]+box[1]-box[0]))
            return false;
        if(fabs(query[2]+query[3]-box[2]-box[3])>(query[3]-query[2]+box[3]-box[2]))
            return false;
        if(fabs(query[4]+query[5]-box[4]-box[5])>(query[5]-query[4]+box[5]-box[4]))
            return false;
        return true;
    }

    bool inside( float query[6] )
    {
        return (box[0]>query[0]&&box[1]<query[1]&&
                box[2]>query[2]&&box[3]<query[3]&&
                box[4]>query[4]&&box[5]<query[5]);
    }

    void getNeighbors(std::list<_octree_cell<KS>*>& nei)
    {
        if(kids[0]!=NULL) return; //only leaf can call this
        float query[6]={box[0]-1e-5f, box[1]+1e-5f,
                        box[2]-1e-5f, box[3]+1e-5f,
                        box[4]-1e-5f, box[5]+1e-5f};
        this->parent->getIntersection(query,nei,this);
    }

    void getBoundaryCells(std::list<_octree_cell<KS>*>& bd, float query[6])
    {
        if( inside(query) && bd_points.empty() ) return;
        
        if(kids[0]==NULL)
            bd.push_back(this);
        else
            for(int i=0;i<KS;i++)
                kids[i]->getBoundaryCells(bd,query);
    }

    void getBoundaryCells(std::list<_octree_cell<KS>*>& bd)
    {
        if(kids[0]!=NULL && parent!=NULL) return; //only root can do this
        for(int i=0;i<KS;i++)
            kids[i]->getBoundaryCells(bd,box);
    }

    bool is_leaf(){
        return (kids[0]==NULL);
    }

    mksum * mk;
    _octree_cell<KS> * parent;
    _octree_cell<KS> * kids[KS];

    MKPTS points;    //m+ points with unknown boundary
    MKPTS bd_points; //known m+ boundary points

    float box[6];
    bool visited;
};

typedef _octree_cell<8> octree_cell; //for 3D

struct octree
{
    octree(mksum * mk):root(mk){ this->mk=mk; }

    //build the octree
    void build()
    {
        //compute bbox of pts
        float pt[3];
        root.box[0]=root.box[2]=root.box[4]=1e10;
        root.box[1]=root.box[3]=root.box[5]=-1e10;
        for(MKPTIT i=mk->pts.begin();i!=mk->pts.end();i++){
            mk->getMPos(*i,pt);
            if(pt[0]<root.box[0]) root.box[0]=pt[0];
            if(pt[0]>root.box[1]) root.box[1]=pt[0];
            if(pt[1]<root.box[2]) root.box[2]=pt[1];
            if(pt[1]>root.box[3]) root.box[3]=pt[1];
            if(pt[2]<root.box[4]) root.box[4]=pt[2];
            if(pt[2]>root.box[5]) root.box[5]=pt[2];
        }
        //store points in the octree now
        root.points.swap(mk->pts);
        root.build(root.box);
    }

    void destroy()
    {
        root.destroy();
    }

    mksum * mk;
    octree_cell root;
};


}

#endif //_BF_OCT_H_

