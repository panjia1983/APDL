#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>

#include <APD/profile.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{
	void test_conlitron_learner_3d_rotation()
	{
		C2A_Model* P = new C2A_Model;
		C2A_Model* Q = new C2A_Model;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler2 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

		std::ofstream out("space_test_3d_rotation.txt");
		asciiWriter(out, contactspace_samples);

		DataVector w(6);
		w[0] = 1; w[1] = 1; w[2] = 1; w[3] = 1; w[4] = 1; w[5] = 1;
		MulticonlitronLearner learner(w, 0.01);
		///////////////////////////////////////
		learner.use_approximate_dist = false;
		///////////////////////////////////////

		tools::Profiler::Begin("learn without approximate knn");
		learner.learn(contactspace_samples, 6);
		tools::Profiler::End("learn without approximate knn");

		std::cout << learner.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << std::endl;
	}

	void test_conlitron_learner_3d_rotation_approximate_knn()
	{
		C2A_Model* P = new C2A_Model;
		C2A_Model* Q = new C2A_Model;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceSE3Euler2 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_3d_rotation.txt");
		asciiWriter(out, contactspace_samples);

		DataVector w(6);
		w[0] = 1; w[1] = 1; w[2] = 1; w[3] = 1; w[4] = 1; w[5] = 1;
		MulticonlitronLearner learner(w, 0.01);
		///////////////////////////////////////
		learner.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn");
		learner.learn(contactspace_samples, 6);
		tools::Profiler::End("learn with approximate knn");

		std::cout << learner.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << std::endl;





		MulticonlitronLearner learner2(w, 0.01, 1, 1);
		///////////////////////////////////////
		learner2.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn");
		learner2.learn(contactspace_samples, 6);
		tools::Profiler::End("learn with approximate knn");

		std::cout << learner2.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner2) << " " << errorRatioOnGrid(contactspace, learner2, 5) << std::endl;




		MulticonlitronLearner learner3(w, 0.01, 5, 5);
		///////////////////////////////////////
		learner3.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn");
		learner3.learn(contactspace_samples, 6);
		tools::Profiler::End("learn with approximate knn");

		std::cout << learner3.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner3) << " " << errorRatioOnGrid(contactspace, learner3, 5) << std::endl;


		MulticonlitronLearner learner4(w, 0.01, 10, 10);
		///////////////////////////////////////
		learner4.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn");
		learner4.learn(contactspace_samples, 6);
		tools::Profiler::End("learn with approximate knn");

		std::cout << learner4.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner4) << " " << errorRatioOnGrid(contactspace, learner4, 5) << std::endl;


	}
}

void main()
{
	APDL::tools::Profiler::Start();

	// APDL::test_conlitron_learner_3d_rotation();
	APDL::test_conlitron_learner_3d_rotation_approximate_knn();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}