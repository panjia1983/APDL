#ifndef MESH_IO_H
#define MESH_IO_H

#include <C2A/C2A.h>
#include <string>

namespace APDL
{
	void readObjFile(C2A_Model* model, const std::string& obj_file);
	
	void readOffFile(C2A_Model* model, const std::string& off_file);
	
	void readTriFile(C2A_Model* model, const std::string& tris_file);
	
}

#endif