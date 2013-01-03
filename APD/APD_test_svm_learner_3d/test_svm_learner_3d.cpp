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

	void test_svm_learner_3d_rings()
	{

		{
			C2A_Model* P = new C2A_Model;
			C2A_Model* Q = new C2A_Model;
			readObjFile(P, "../data/models/Bullet/ringz.obj");
			readObjFile(Q, "../data/models/Bullet/ringz.obj");

			P->ComputeRadius();
			Q->ComputeRadius();

			ContactSpaceR3 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);

			std::ofstream out("space_test_3d_ringz.txt");
			asciiWriter(out, contactspace_samples);

			std::ifstream in("space_test_3d_ringz.txt");
			std::vector<ContactSpaceSampleData> contactspace_samples2;
			asciiReader(in, contactspace_samples2);

			std::ofstream out2("space_test_3d_ringz2.txt");
			asciiWriter(out2, contactspace_samples2);		

			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 
			//learner.setGamma(0.1 * contactspace.getScaler().getScale());
			//std::cout << 0.1 * contactspace.getScaler().getScale() << std::endl;



			std::ofstream scaler_file("scaler_3d_ringz.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_3d_ringz.txt");

			std::cout << "model saved" << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << std::endl;

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;
		}
	}
}

void main()
{
	// APDL::test_svm_learner_3d();
	APDL::test_svm_learner_3d_rings();
}