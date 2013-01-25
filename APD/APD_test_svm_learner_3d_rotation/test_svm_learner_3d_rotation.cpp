#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/mesh_io.h>
#include <APD/active_learning.h>

#include <APD/profile.h>

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

			tools::Profiler::Begin("learner");
			ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);
			tools::Profiler::End("learner");
				
			std::ofstream out("space_test_3d_rotation.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 


			std::ofstream scaler_file("scaler_3d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			
			tools::Profiler::Begin("learner");
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			tools::Profiler::End("learner");
			learner.save("model_3d_rotation.txt");

			std::cout << "model saved" << std::endl;


			std::vector<ContactSpaceSampleData> test_samples = contactspace.uniform_sample(10000);
			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << empiricalErrorRatio(test_samples, learner) << std::endl;

		    //std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " <<  errorRatioOnGrid(contactspace, learner, 5) << " ";

			//for(std::size_t i = 0; i < contactspace_samples.size(); ++i)
			//	std::cout << "(" << results[i].label << "," << contactspace_samples[i].col << ")";
			//std::cout << std::endl;

			delete P;
			delete Q;
		}

	}


	void test_svm_learner_3d_rotation_quat()
	{

		{
			C2A_Model* P = NULL;
			C2A_Model* Q = NULL;
			readOffFile(P, "../data/cup.off");
			readOffFile(Q, "../data/spoon.off");

			P->ComputeRadius();
			Q->ComputeRadius();

			tools::Profiler::Begin("learner");
			ContactSpaceSE3Quat contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(10000);
			tools::Profiler::End("learner");
				
			std::ofstream out("space_test_3d_rotation_quat.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 


			std::ofstream scaler_file("scaler_3d_rotation_quat.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			
			tools::Profiler::Begin("learner");
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			tools::Profiler::End("learner");
			learner.save("model_3d_rotation_quat.txt");

			std::cout << "model saved" << std::endl;

			// std::vector<ContactSpaceSampleData> test_samples = contactspace.uniform_sample(10000);
			// std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << empiricalErrorRatio(test_samples, learner) << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " <<  errorRatioOnGrid(contactspace, learner, 5) << " ";


			delete P;
			delete Q;
		}

	}



	void test_svm_learner_3d_rotation_carseat()
	{

		{
			C2A_Model* P = NULL;
			C2A_Model* Q = NULL;
			readObjFile(P, "../data/car_seat/car.obj");
			readObjFile(Q, "../data/car_seat/seat.obj");

			P->ComputeRadius();
			Q->ComputeRadius();

			tools::Profiler::Begin("learner");
			ContactSpaceSE3Euler contactspace(P, Q, 0.05 * (P->radius + Q->radius));
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(500000);
			tools::Profiler::End("learner");
				
			std::ofstream out("space_test_3d_rotation.txt");
			asciiWriter(out, contactspace_samples);
			
			SVMLearner learner;
			learner.setC(50);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 


			std::ofstream scaler_file("scaler_3d_rotation.txt");
			scaler_file << contactspace.getScaler() << std::endl;
			
			tools::Profiler::Begin("learner");
			learner.learn(contactspace_samples, contactspace.active_data_dim());
			tools::Profiler::End("learner");
			learner.save("model_3d_rotation.txt");

			std::cout << "model saved" << std::endl;


			std::vector<ContactSpaceSampleData> test_samples = contactspace.uniform_sample(10000);
			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << empiricalErrorRatio(test_samples, learner) << std::endl;

		    //std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " <<  errorRatioOnGrid(contactspace, learner, 5) << " ";

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
	APDL::tools::Profiler::Start();

	// APDL::test_svm_learner_3d_rotation();

	// APDL::test_svm_learner_3d_rotation_quat();

	APDL::test_svm_learner_3d_rotation_carseat();

	APDL::tools::Profiler::Stop();

	APDL::tools::Profiler::Status();
}