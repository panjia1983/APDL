//------------------------------------------------------------------------------
//  Copyright 2007 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "ModelEdge.h"

namespace libm3d
{

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CModelEdge::CModelEdge(int start, int end)
{
    m_Key[0]=start;
    m_Key[1]=end;
    m_Next=NULL; //m_Kids[0]=m_Kids[1]=NULL;
    //m_Mid=-1;
	m_Fid[0]=m_Fid[1]=-1;
}

}