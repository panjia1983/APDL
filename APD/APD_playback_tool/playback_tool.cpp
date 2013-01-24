#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

void main()
{
	int i = 0;

	std::string rotation_file_name("cupspoon_rotations.txt");
	std::string obj_file_name("cupspoon_");
	std::string file_name;


	std::ofstream rot_file(rotation_file_name.c_str());

	std::vector<std::string> dir_names;

	dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result1/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result2/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result3/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result3/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result4/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result5/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result6/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon/result7/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon2/result0/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon2/result1/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon2/result2/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon3/result0/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon3/result1/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon3/result2/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon4/result0/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon4/result1/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon4/result2/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon4/result3/"));
	//dir_names.push_back(std::string("../../../mesh_m.win32/cupspoon4/result4/"));

	for(int k = 0; k < dir_names.size(); ++k)
	{
		std::string dir_name = dir_names[k];
		file_name = dir_name + rotation_file_name;
		std::ifstream in_rot_file(file_name.c_str());
		if(!in_rot_file.is_open())
		{
			std::cout << "Failed to open the file " << file_name << std::endl;
			return;
		}

		std::string line;
		int local_id = 0;
		while(!in_rot_file.eof())
		{
			std::getline(in_rot_file, line, '\n');
			std::cout << line << std::endl;
			if(line.size() == 0) break;

			std::istringstream reader(line);

			double a, b, c, d;
			reader >> a >> b >> c >> d;

			std::stringstream ss;
			ss << local_id;
			std::string src_id;
			ss >> src_id;

			std::stringstream ss2 ;
			ss2 << i;
			std::string dest_id;
			ss2 >> dest_id;

			std::string src_file = dir_name + obj_file_name + src_id + ".obj";
			std::string dest_file = obj_file_name + dest_id + ".obj";

			std::cout << src_file << std::endl;
			std::cout << dest_file << std::endl;

			std::ifstream src(src_file.c_str());
			std::ofstream dest(dest_file.c_str());

			if(!src.is_open())
			{
				std::cout << "Failed to open the file " << src_file << std::endl;
				break;
			}

			if(!dest.is_open())
			{
				std::cout << "Failed to open the file " << dest_file << std::endl;
				break;
			}

			rot_file << a << " " << b << " " << c << " " << d << std::endl;

			dest << src.rdbuf();

			local_id++;
			i++;
		}

	}

}




void main1()
{
	int i = 0;

	std::string rotation_file_name("teeth_rotations.txt");
	std::string obj_file_name("teeth_");
	std::string file_name;


	std::ofstream rot_file(rotation_file_name.c_str());

	std::vector<std::string> dir_names;

	dir_names.push_back(std::string("../../../mesh_m.win32/teeth/result0/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth/result1/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth/result2/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth/result3/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth2/result0/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth2/result1/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth2/result2/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth3/result0/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth3/result1/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth3/result2/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth4/result0/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth4/result1/"));
	dir_names.push_back(std::string("../../../mesh_m.win32/teeth4/result2/"));

	for(int k = 0; k < dir_names.size(); ++k)
	{
		std::string dir_name = dir_names[k];
		file_name = dir_name + rotation_file_name;
		std::ifstream in_rot_file(file_name.c_str());
		if(!in_rot_file.is_open())
		{
			std::cout << "Failed to open the file " << file_name << std::endl;
			return;
		}

		std::string line;
		int local_id = 0;
		while(!in_rot_file.eof())
		{
			std::getline(in_rot_file, line, '\n');
			std::cout << line << std::endl;
			if(line.size() == 0) break;

			std::istringstream reader(line);

			double a, b, c, d;
			reader >> a >> b >> c >> d;

			std::stringstream ss;
			ss << local_id;
			std::string src_id;
			ss >> src_id;

			std::stringstream ss2 ;
			ss2 << i;
			std::string dest_id;
			ss2 >> dest_id;

			std::string src_file = dir_name + obj_file_name + src_id + ".obj";
			std::string dest_file = obj_file_name + dest_id + ".obj";

			std::cout << src_file << std::endl;
			std::cout << dest_file << std::endl;

			std::ifstream src(src_file.c_str());
			std::ofstream dest(dest_file.c_str());

			if(!src.is_open())
			{
				std::cout << "Failed to open the file " << src_file << std::endl;
				break;
			}

			if(!dest.is_open())
			{
				std::cout << "Failed to open the file " << dest_file << std::endl;
				break;
			}

			rot_file << a << " " << b << " " << c << " " << d << std::endl;

			dest << src.rdbuf();

			local_id++;
			i++;
		}

	}

}