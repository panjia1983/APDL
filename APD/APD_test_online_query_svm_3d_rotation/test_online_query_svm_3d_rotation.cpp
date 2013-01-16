#include <APD/online_query.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>

#include <APD/profile.h>

namespace APDL
{
	extern double distance_weight[7];

	void test_online_query_svm_3d_rotation()
	{
		bool use_euler = true;

		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readObjFile(P, "../data/models/CupSpoon/Cup.obj");
		readObjFile(Q, "../data/models/CupSpoon/Spoon.obj");

		P->ComputeRadius();
		Q->ComputeRadius();

		{
			double Ix, Iy, Iz;
			inertia_weight(Q, Ix, Iy, Iz);
			distance_weight[0] = 1; distance_weight[1] = 1; distance_weight[2] = 1;
			std::cout << Ix << " " << Iy << " " << Iz << std::endl;
			if(use_euler)
			{
				distance_weight[3] = Ix; distance_weight[4] = Iy; distance_weight[5] = Iz;
			}
			else
			{
				distance_weight[3] = 1; distance_weight[4] = Ix; distance_weight[5] = Iy; distance_weight[6] = Iz;
			}
		}

		ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::ofstream scaler_file("scaler_3d_rotation_cupspoon.txt");
		scaler_file << contactspace.getScaler() << std::endl;

		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		SVMLearner learner;
		learner.setDim(contactspace.active_data_dim());
		learner.setC(20);
		learner.setScaler(contactspace.getScaler());
		learner.setUseScaler(true);
		learner.setGamma(50); 

		learner.learn(contactspace_samples, contactspace.active_data_dim());
		learner.save("model.txt");

		std::cout << contactspace_samples.size() << " " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << std::endl;


		std::vector<ContactSpaceSampleData> query_samples = contactspace.uniform_sample(40);

		flann::HierarchicalClusteringIndex<ContactSpaceSE3Euler::DistanceType>* query_index = learner.constructIndexOfSupportVectorsForQuery<ContactSpaceSE3Euler, flann::HierarchicalClusteringIndex, flann::HierarchicalClusteringIndexParams>();

		std::ofstream query_file("query_results.txt");

		for(std::size_t i = 0; i < query_samples.size(); ++i)
		{
			QueryResult apprx_PD = PD_query(learner, contactspace, query_index, query_samples[i].v);
			QueryResult apprx_PD2 = PD_query2(learner, contactspace, query_index, query_samples[i].v);

			query_file << apprx_PD.v[0] << " " << apprx_PD.v[1] << " " << apprx_PD.v[2] << " " << apprx_PD.PD << std::endl;
			query_file << apprx_PD2.v[0] << " " << apprx_PD2.v[1] << " " << apprx_PD2.v[2] << " " << apprx_PD2.PD << std::endl;
			query_file << std::endl;
		}

		delete query_index;
	}
}

void main()
{
	APDL::tools::Profiler::Start();
	APDL::test_online_query_svm_3d_rotation();
	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}