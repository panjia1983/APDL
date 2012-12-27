//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _MODEL_GRAPH_EDGE_H_
#define _MODEL_GRAPH_EDGE_H_

#include "model.h"


class CModelEdge  
{
public:
	//////////////////////////////////////////////////////////////////////
	// Constructor/Destructors
	CModelEdge(int start, int end);

	//////////////////////////////////////////////////////////////////////
	// Access Function
	bool isEndPt( int key ) const{
		for( int i=0;i<2;i++)
			if( key==m_Key[i] ) return true;
		return false;
	}
	int getStartPt()  const { return m_Key[0]; }
	int getEndPt()    const { return m_Key[1]; }

	//List Access
	CModelEdge * getNext() const { return m_Next; }
	void setNext(CModelEdge * e) { m_Next=e; }

	void setFacet0(int id) {m_Fid[0]=id;}
	void setFacet1(int id) {m_Fid[1]=id;}

	int getFacet0() const {return m_Fid[0];}
	int getFacet1() const {return m_Fid[1];}

	Vector3d& getInNormal(int id){ return m_InNormal[id]; }
	Vector3d& getV(){ return m_V; }
	//////////////////////////////////////////////////////////////////////
	// Private Stuff
private:

	int m_Key[2]; ///< 0->start, 1->end
	int m_Fid[2]; //facet id

	Vector3d m_V;  //points from start to end
	Vector3d m_InNormal[2]; //inface normal

	//list link
	CModelEdge * m_Next;
};


#endif // _MODEL_GRAPH_EDGE_H_
