#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/mesh_io.h>
#include <APD/active_learning.h>

void* user_conlitron_model;
double* user_conlitron_data;


namespace APDL
{
	void test_svm_learner_3d()
	{

		{
			C2A_Model* P = new C2A_Model;
			C2A_Model* Q = new C2A_Model;
			readOffFile(P, "../data/cup.off");
			readOffFile(Q, "../data/spoon.off");

			P->ComputeRadius();
			Q->ComputeRadius();

			ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);
				
			std::ofstream out("space_test_3d.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 
			//learner.setGamma(0.1 * contactspace.getScaler().getScale());
			//std::cout << 0.1 * contactspace.getScaler().getScale() << std::endl;



			std::ofstream scaler_file("scaler_3d.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_3d.txt");

			std::cout << "model saved" << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << std::endl;

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;
		}

	}

	//void test()
	//{
	//	std::vector<DataVector> points = readPM3dFile("../data/grate12.m+");
	//	for(std::size_t i = points.size() - 10; i < points.size(); ++i)
	//	{
	//		for(std::size_t j = 0; j < points[i].dim(); ++j)
	//			std::cout << points[i][j] << " ";
	//		std::cout << std::endl;
	//	}
	//	std::cout << points.size() << std::endl;
	//}

}

void main()
{
	// APDL::test();
	APDL::test_svm_learner_3d();
}