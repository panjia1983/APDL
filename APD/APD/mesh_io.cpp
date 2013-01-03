#include "mesh_io.h"
#include <iostream>
#include <fstream>
namespace APDL
{
	void readObjFile(C2A_Model* model, const std::string& obj_file)
	{
		std::ifstream in(obj_file.c_str());
		if(!in.is_open())
		{
			std::cout << "Failed to open the file " << obj_file << std::endl;
			return;
		}
		
		std::string type;
		std::vector<double> v;
		std::vector<int> id;
		
		while(!in.eof())
		{
			in >> type;
			if(type == "v")
			{
				double a, b, c;
				in >> a >> b >> c;
				v.push_back(a);
				v.push_back(b);
				v.push_back(c);
			}
			else if(type == "f")
			{
				int i1, i2, i3;
				in >> i1 >> i2 >> i3;
				id.push_back(i1 - 1);
				id.push_back(i2 - 1);
				id.push_back(i3 - 1);
			}

			getline(in, type);
		}
		
		int num_v = v.size() / 3;
		int num_f = id.size() / 3;
		
		model->BeginModel();
		
		PQP_REAL p1[3], p2[3], p3[3];
		int i1, i2, i3;
		for(int i = 0; i < num_f; ++i)
		{
			i1 = id[3 * i + 0];
			i2 = id[3 * i + 1];
			i3 = id[3 * i + 2];
			p1[0] = v[i1 * 3];
			p1[1] = v[i1 * 3 + 1];
			p1[2] = v[i1 * 3 + 2];
			
			p2[0] = v[i2 * 3];
			p2[1] = v[i2 * 3 + 1];
			p2[2] = v[i2 * 3 + 2];
			
			p3[0] = v[i3 * 3];
			p3[1] = v[i3 * 3 + 1];
			p3[2] = v[i3 * 3 + 2];
			
			model->AddTri(p1, p2, p3, i, i1, i2, i3);
		}
		
		model->EndModel();
	}
	
	void readOffFile(C2A_Model* model, const std::string& off_file)
	{
		std::ifstream in(off_file.c_str());
		if(!in.is_open())
		{
			std::cout << "Failed to open the file " << off_file << std::endl;
			return;
		}
		
		std::string str;
		int num_v, num_f;
		int dummy;
		in >> str >> num_v >> num_f >> dummy;
		
		PQP_REAL* p = new PQP_REAL[3 * num_v];
		
		model->BeginModel();
		
		PQP_REAL a, b, c;
		int i1, i2, i3;
		
		for(int i = 0; i < num_v; ++i)
		{
			in >> a >> b >> c;
			p[3 * i + 0] = a;
			p[3 * i + 1] = b;
			p[3 * i + 2] = c;
		}
		
		PQP_REAL p1[3], p2[3], p3[3];
		for(int i = 0; i < num_f; ++i)
		{
			in >> dummy >> i1 >> i2 >> i3; // dummy should be 3
			p1[0] = p[i1 * 3];
			p1[1] = p[i1 * 3 + 1];
			p1[2] = p[i1 * 3 + 2];
			
			p2[0] = p[i2 * 3];
			p2[1] = p[i2 * 3 + 1];
			p2[2] = p[i2 * 3 + 2];
			
			p3[0] = p[i3 * 3];
			p3[1] = p[i3 * 3 + 1];
			p3[2] = p[i3 * 3 + 2];
			
			model->AddTri(p1, p2, p3, i, i1, i2, i3);
		}
		
		delete [] p;
		
		model->EndModel();
	}
	
	void readTriFile(C2A_Model* model, const std::string& tris_file)
	{
		std::ifstream in(tris_file.c_str());
		if(!in.is_open())
		{
			std::cout << "Failed to open the file " << tris_file << std::endl;
			return;
		}
		
		std::string str;
		int num_v, num_f;
		in >> str >> num_v >> num_f;
		
		PQP_REAL* p = new PQP_REAL[3 * num_v];
		
		PQP_REAL a, b, c;
		int i1, i2, i3;
		
		for(int i = 0; i < num_v; ++i)
		{
			in >> a >> b >> c;
			p[3 * i + 0] = a;
			p[3 * i + 1] = b;
			p[3 * i + 2] = c;
		}
		
		model->BeginModel();
		
		PQP_REAL p1[3], p2[3], p3[3];
		for(int i = 0; i < num_f; ++i)
		{
			in >> i1 >> i2 >> i3;
			p1[0] = p[i1 * 3];
			p1[1] = p[i1 * 3 + 1];
			p1[2] = p[i1 * 3 + 2];
			
			p2[0] = p[i2 * 3];
			p2[1] = p[i2 * 3 + 1];
			p2[2] = p[i2 * 3 + 2];
			
			p3[0] = p[i3 * 3];
			p3[1] = p[i3 * 3 + 1];
			p3[2] = p[i3 * 3 + 2];
			
			model->AddTri(p1, p2, p3, i, i1, i2, i3);
		}
		
		delete [] p;
		
		model->EndModel();
	}
}