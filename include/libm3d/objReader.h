//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <ctype.h>
#include <math.h>

namespace libm3d 
{

class Vpt
{
public:
    Vpt(){ x=y=z=nx=ny=nz=0; }

    void normalize()
    {
        float norm=(float)sqrt(nx*nx+ny*ny+nz*nz);
        nx/=norm;
        ny/=norm;
        nz/=norm;
    }

    float x, y, z, nx, ny, nz;
};

class V
{
public:
    V(){ x=y=z=0; }
    float x, y, z;
};

class polygon
{
public:
    std::list<int> pts;
    std::list<int> normals;
};

class objModel
{
public:
    objModel() { }
    
    void compute_v_normal()
    {
        for(std::list<polygon>::iterator i=polys.begin();i!=polys.end();i++)
        {
            std::list<int>::iterator ni=i->normals.begin();
            std::list<int>::iterator pi=i->pts.begin();

            for( ;pi!=i->pts.end();pi++,ni++){
                V& n=normals[*ni];
                Vpt& pt=pts[*pi];
                pt.nx+=n.x;
                pt.ny+=n.y;
                pt.nz+=n.z;
            }
        }//end for

        //normalize
        for(std::vector<Vpt>::iterator i=pts.begin();i!=pts.end();i++)
        {
            i->normalize();
        }
    }

    std::vector<Vpt> pts;
    std::vector<V> normals;
    std::list<polygon> polys;
};

class objReader
{
public:
    objReader(const string& name){ m_filename=name; }

    bool Read(){
        ifstream in(m_filename.c_str());
        if( !in.good() ){
            cerr<<"Can't open file "<<m_filename<<endl;
            return false;
        }
        bool r=Read(in);
        in.close();
        return r;
    }

    const objModel& getModel() const { return data; }
    objModel& getModel() { return data; }

private:

    bool Read(istream& in){

        string tmp;
        
        //read pts
        while(true){
            in>>tmp;
            if( tmp=="f" ) break;
            if( tmp=="v" ){
                Vpt pt;
                in>>pt.x>>pt.y>>pt.z;
                data.pts.push_back(pt);
            }
            else if( tmp=="vn" ){
                V pt;
                in>>pt.x>>pt.y>>pt.z;
                data.normals.push_back(pt);
            }
            getline(in,tmp);
        }

        //read faces
        polygon poly;
        do{

            in>>tmp;
            if( in.eof() ) break;
            if( isdigit(tmp[0]) ){ //this defines a vetex
                
                int pos1=tmp.find('/');
                int pos2=tmp.rfind('/');

                int id_v=atoi(tmp.substr(0,pos1).c_str())-1;
                int id_n=atoi(tmp.substr(pos2+1).c_str())-1;

                poly.pts.push_back(id_v);
                poly.normals.push_back(id_n);
            }
            else if( tmp=="f" )
            {
                data.polys.push_back(poly);
                poly.pts.clear();
                poly.normals.clear();
            }
            else 
                getline(in,tmp);

        }
        while( !in.eof() );

        data.polys.push_back(poly);
        data.compute_v_normal();

        return true;
    }

    string m_filename;
    objModel data;
};

}

#endif //_OBJ_READER_H_



