#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/math_utility.h>

#include <APD/profile.h>

namespace APDL
{
	void playback()
	{
		std::vector<std::vector<DataVector> > frames;

		std::string frame_file_name = "../data/models/Teeth/teeth.ani";

		bool use_euler = true;
		frames = readAnimationFile(frame_file_name, use_euler);

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/Teeth/lo_03de_new.obj");
		readObjFile(Q, "../data/models/Teeth/up_03de_new.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));

		for(std::size_t i = 0; i < frames.size(); ++i)
		{
			DataVector q = relative3D(frames[i][0], frames[i][1]);

			//for(int j = 0; j < q.dim(); ++j)
			//	std::cout << q[j] << " ";
			//std::cout << std::endl;
		}

	}
}

void main()
{
	APDL::playback();
}