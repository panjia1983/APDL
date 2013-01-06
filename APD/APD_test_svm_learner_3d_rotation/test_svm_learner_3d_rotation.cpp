#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/mesh_io.h>
#include <APD/active_learning.h>


namespace APDL
{
	void test_svm_learner_3d_rotation()
	{

		{
			C2A_Model* P = NULL;
			C2A_Model* Q = NULL;
			readOffFile(P, "../data/cup.off");
			readOffFile(Q, "../data/spoon.off");

			P->ComputeRadius();
			Q->ComputeRadius();

			ContactSpaceSE3Euler2 contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(100000);
				
			std::ofstream out("space_test_3d_rotation.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 


			std::ofstream scaler_file("scaler_3d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_3d_rotation.txt");

			std::cout << "model saved" << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 5) << std::endl;


			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			delete P;
			delete Q;
		}

	}

}

void main()
{
	APDL::test_svm_learner_3d_rotation();
}