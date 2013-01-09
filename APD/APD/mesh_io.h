#ifndef MESH_IO_H
#define MESH_IO_H

#include <C2A/C2A.h>
#include <string>
#include "math_utility.h"

namespace APDL
{
	void readObjFile(C2A_Model*& model, const std::string& obj_file);

	void readObjFiles(std::vector<C2A_Model*>& models, const std::string& obj_file);

	void readSE3Model(std::vector<std::pair<C2A_Model*, Quaternion> >& cspace, const std::string& model_file);
	
	void readOffFile(C2A_Model*& model, const std::string& off_file);
	
	void readTriFile(C2A_Model*& model, const std::string& tris_file);
	
}

#endif