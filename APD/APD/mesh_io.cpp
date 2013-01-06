#include "mesh_io.h"
#include <iostream>
#include <fstream>
#include <sstream>

namespace APDL
{
	void readObjFiles(std::vector<C2A_Model*>& models, const std::string& obj_file)
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

		int first_o = true;
		int start_id = 1;
		double center[3] = {0, 0, 0};

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
				center[0] += a;
				center[1] += b;
				center[2] += c;
			}
			else if(type == "f")
			{
				int i1, i2, i3;
				in >> i1 >> i2 >> i3;
				id.push_back(i1 - start_id);
				id.push_back(i2 - start_id);
				id.push_back(i3 - start_id);
			}
			else if(type == "o")
			{
				if(!first_o)
				{
					int num_v = v.size() / 3;
					int num_f = id.size() / 3;

					for(int i = 0; i < id.size(); ++i)
					{
						if(id[i] < 0 || id[i] >= num_v)
							std::cout << "error " << id[i] << std::endl;
					}

					center[0] /= num_v;
					center[1] /= num_v;
					center[2] /= num_v;

					double radius = 0;
					for(std::size_t i = 0; i < num_v; ++i)
					{
						double r = 0;
						r += (v[3 * i] - center[0]) * (v[3 * i] - center[0]);
						r += (v[3 * i + 1] - center[1]) * (v[3 * i + 1] - center[1]);
						r += (v[3 * i + 2] - center[2]) * (v[3 * i + 2] - center[2]);
						if(r > radius) radius = r;
					}
					radius = sqrt(radius);

					C2A_Model* model = new C2A_Model;

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
					model->com[0] = center[0];
					model->com[1] = center[1];
					model->com[2] = center[2];
					model->radius = radius;

					center[0] = 0; center[1] = 0; center[2] = 0;

					models.push_back(model);

					start_id += num_v;
					v.clear();
					id.clear();
				}
				else
					first_o = false;

			}

			getline(in, type);
		}

		int num_v = v.size() / 3;
		int num_f = id.size() / 3;
		center[0] /= num_v;
		center[1] /= num_v;
		center[2] /= num_v;

		double radius = 0;
		for(std::size_t i = 0; i < num_v; ++i)
		{
			double r = 0;
			r += (v[3 * i] - center[0]) * (v[3 * i] - center[0]);
			r += (v[3 * i + 1] - center[1]) * (v[3 * i + 1] - center[1]);
			r += (v[3 * i + 2] - center[2]) * (v[3 * i + 2] - center[2]);
			if(r > radius) radius = r;
		}
		radius = sqrt(radius);

		C2A_Model* model = new C2A_Model;
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
		model->com[0] = center[0];
		model->com[1] = center[1];
		model->com[2] = center[2];
		model->radius = radius;

		models.push_back(model);
	}


	void readObjFile(C2A_Model*& model, const std::string& obj_file)
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
		
		double center[3] = {0, 0, 0};

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
				center[0] += a;
				center[1] += b;
				center[2] += c;
			}
			else if(type == "f")
			{
				for(int i = 0; i < 3; ++i)
				{
					std::string face;
					std::string face2;
					in >> face;
					for(int k = 0; k < face.size(); ++k)
					{
						if(face[k] == '/') break;
						face2.push_back(face[k]);
					}
					int fid;
					std::istringstream convert(face2);
					convert >> fid;
					// std::cout << face << " " << fid << std::endl;
					id.push_back(fid - 1);
				}


				//int i1, i2, i3;
				//in >> i1 >> i2 >> i3;
				//id.push_back(i1 - 1);
				//id.push_back(i2 - 1);
				//id.push_back(i3 - 1);
			}

			getline(in, type);
		}
		
		int num_v = v.size() / 3;
		int num_f = id.size() / 3;
		center[0] /= num_v;
		center[1] /= num_v;
		center[2] /= num_v;

		double radius = 0;
		for(std::size_t i = 0; i < num_v; ++i)
		{
			double r = 0;
			r += (v[3 * i] - center[0]) * (v[3 * i] - center[0]);
			r += (v[3 * i + 1] - center[1]) * (v[3 * i + 1] - center[1]);
			r += (v[3 * i + 2] - center[2]) * (v[3 * i + 2] - center[2]);
			if(r > radius) radius = r;
		}
		radius = sqrt(radius);

		model = new C2A_Model;

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
		model->com[0] = center[0];
		model->com[1] = center[1];
		model->com[2] = center[2];
		model->radius = radius;
	}
	
	void readOffFile(C2A_Model*& model, const std::string& off_file)
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

		model = new C2A_Model;
		model->BeginModel();
		
		PQP_REAL a, b, c;
		int i1, i2, i3;

		double center[3] = {0, 0, 0};
		
		for(int i = 0; i < num_v; ++i)
		{
			in >> a >> b >> c;
			p[3 * i + 0] = a;
			p[3 * i + 1] = b;
			p[3 * i + 2] = c;
			center[0] += a;
			center[1] += b;
			center[2] += c;
		}
	
		center[0] /= num_v;
		center[1] /= num_v;
		center[2] /= num_v;

		double radius = 0;
		for(std::size_t i = 0; i < num_v; ++i)
		{
			double r = 0;
			r += (p[3 * i] - center[0]) * (p[3 * i] - center[0]);
			r += (p[3 * i + 1] - center[1]) * (p[3 * i + 1] - center[1]);
			r += (p[3 * i + 2] - center[2]) * (p[3 * i + 2] - center[2]);
			if(r > radius) radius = r;
		}
		radius = sqrt(radius);

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
		model->com[0] = center[0];
		model->com[1] = center[1];
		model->com[2] = center[2];
		model->radius = radius;

	}
	
	void readTriFile(C2A_Model*& model, const std::string& tris_file)
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

		double center[3] = {0, 0, 0};

		for(int i = 0; i < num_v; ++i)
		{
			in >> a >> b >> c;
			p[3 * i + 0] = a;
			p[3 * i + 1] = b;
			p[3 * i + 2] = c;
			center[0] += a;
			center[1] += b;
			center[2] += c;
		}

		center[0] /= num_v;
		center[1] /= num_v;
		center[2] /= num_v;

		double radius = 0;
		for(std::size_t i = 0; i < num_v; ++i)
		{
			double r = 0;
			r += (p[3 * i] - center[0]) * (p[3 * i] - center[0]);
			r += (p[3 * i + 1] - center[1]) * (p[3 * i + 1] - center[1]);
			r += (p[3 * i + 2] - center[2]) * (p[3 * i + 2] - center[2]);
			if(r > radius) radius = r;
		}
		radius = sqrt(radius);

		
		model = new C2A_Model;
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
		model->com[0] = center[0];
		model->com[1] = center[1];
		model->com[2] = center[2];
		model->radius = radius;

	}
}