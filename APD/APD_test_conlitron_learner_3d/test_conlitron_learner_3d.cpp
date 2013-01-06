#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/active_learning.h>

#include <APD/profile.h>


namespace APDL
{
	void test_conlitron_learner_3d()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

		std::ofstream out("space_test_3d.txt");
		asciiWriter(out, contactspace_samples);

		DataVector w(3);
		w[0] = 1; w[1] = 1; w[2] = 1;
		MulticonlitronLearner learner(w, 0.01);
		///////////////////////////////////////
		learner.use_approximate_dist = false;
		///////////////////////////////////////

		tools::Profiler::Begin("learn without approximate knn");
		learner.learn(contactspace_samples, 3);
		tools::Profiler::End("learn without approximate knn");

		std::cout << learner.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 50) << std::endl;

		delete P;
		delete Q;
	}

	void test_conlitron_learner_3d_approximate_knn()
	{
		C2A_Model* P = NULL;
		C2A_Model* Q = NULL;
		readOffFile(P, "../data/cup.off");
		readOffFile(Q, "../data/spoon.off");

		P->ComputeRadius();
		Q->ComputeRadius();

		ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
		std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

		std::ofstream out("space_test_3d.txt");
		asciiWriter(out, contactspace_samples);

		DataVector w(3);
		w[0] = 1; w[1] = 1; w[2] = 1;
		MulticonlitronLearner learner(w, 0.01, 0, 0);
		///////////////////////////////////////
		learner.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn");
		learner.learn(contactspace_samples, 3);
		tools::Profiler::End("learn with approximate knn");

		std::cout << learner.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 50) << std::endl;



		MulticonlitronLearner learner2(w, 0.01, 1, 0);
		///////////////////////////////////////
		learner2.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn2");
		learner2.learn(contactspace_samples, 3);
		tools::Profiler::End("learn with approximate knn2");

		std::cout << learner2.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner2) << " " << errorRatioOnGrid(contactspace, learner2, 50) << std::endl;




		MulticonlitronLearner learner3(w, 0.01, 20, 20);
		///////////////////////////////////////
		learner3.use_approximate_dist = true;
		///////////////////////////////////////

		tools::Profiler::Begin("learn with approximate knn3");
		learner3.learn(contactspace_samples, 3);
		tools::Profiler::End("learn with approximate knn3");

		std::cout << learner3.model.numOfHyperPlanes() << std::endl;

		std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner3) << " " << errorRatioOnGrid(contactspace, learner3, 50) << std::endl;

		delete P;
		delete Q;

	}
}

void main()
{
	APDL::tools::Profiler::Start();

	// APDL::test_conlitron_learner_3d();
	APDL::test_conlitron_learner_3d_approximate_knn();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}