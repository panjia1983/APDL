#ifndef _PB_MINKOWSKI_POINTS_
#define _PB_MINKOWSKI_POINTS_

#include <Vector.h>
#include <Point.h>

#include <string>
#include <cassert>

#include "objReader.h"

namespace libm3d 
{

struct point
{
    mathtool::Point3d p;  //position
    mathtool::Vector3d n; //normal
};

struct points
{
    points(){
        pts=NULL;
        size=ref=0;
    }

    ~points(){ delete [] pts; }

    void build(const string& filename){

        if(pts!=NULL) return; //has been built

        objReader reader(filename);
        reader.Read();
        const std::vector<Vpt>& v=reader.getModel().pts;
        size=v.size();
        pts=new point[size];
        assert(pts); //make sure enough memory

        //copy points, we also
        //find a ref point at +x exterme
        //(! may want to change so not both P and Q compute ref)
        float x_max=-1e10;
        ref=0;
        for(unsigned int i=0;i<size;i++){
            pts[i].p.set(&v[i].x);
            pts[i].n.set(&v[i].nx);

            if(pts[i].p[0]>x_max){
                x_max=pts[i].p[0];
                ref=i;
            }
        }//end i
    }

    const point& getRef() const { return pts[ref]; }

    point * pts;
    unsigned int size;
    unsigned int ref; //id to ref point
};

// compute the bounding box of the minkowski sum of P and Q
inline 
void bbox_point_minkowski_sum(const points& P, const points& Q, float box[6])
{
    box[0]=box[2]=box[4]= 1e10;
    box[1]=box[3]=box[5]=-1e10;

    const point& ref=P.getRef();
    for(unsigned int i=0;i<Q.size;i++){
        mathtool::Vector3d offset=Q.pts[i].p-ref.p;

        for(unsigned int j=0;j<P.size;j++){
            mathtool::Point3d np=P.pts[j].p+offset; //
            
            if(box[0]>np[0]) box[0]=np[0];
            if(box[1]<np[0]) box[1]=np[0];
            if(box[2]>np[1]) box[2]=np[1];
            if(box[3]<np[1]) box[3]=np[1];
            if(box[4]>np[2]) box[4]=np[2];
            if(box[5]<np[2]) box[5]=np[2];
        }//end j
    }//end i
}

}

#endif //_PB_MINKOWSKI_POINTS_