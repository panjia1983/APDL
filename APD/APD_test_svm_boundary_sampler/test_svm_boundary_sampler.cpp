#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/active_learning.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{	
	void test_svm_boundary_sampler()
	{
		{
			std::ifstream room_file("../data/rooms_star.dat");

			if(!room_file.is_open())
			{
				std::cerr << "Failed to open the input file." << std::endl;
				return;
			}

			Minkowski_Cspace_2D::Polygon_2 P, Q;

			room_file >> P >> Q;
			room_file.close();

			Polygon p1 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(P);
			Polygon p2 = toPolygon<Minkowski_Cspace_2D::Polygon_2, Minkowski_Cspace_2D::Kernel>(Q);

			ContactSpaceR2 contactspace(p1, p2, 2);
			std::vector<ContactSpaceSampleData> contactspace_samples = contactspace.uniform_sample(1000);

			std::ofstream out("space_test_2d.txt");
			asciiWriter(out, contactspace_samples);

			SVMLearner learner;
			learner.setC(10);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			SpatialTreeParam param0;
			SpatialTreeEParam param;

			std::ofstream scaler_file("scaler_2d.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			if(1) // not scaled
			{
				param.max_depth = 10;
				param.initial_depth = 4;
				param.stop_abs_diff = 0.2;
				param.stop_related_diff = 0.1;
				param.epsilon = 0;
				param.result_eps = 0;

				param0.max_depth = 6;
				param0.initial_depth = 4;
				param0.result_eps = 0.5;
			}
			else // scaled
			{
				learner.setUseScaler(true);
				learner.setGamma(0.5 * contactspace.getScaler().getScale() * 10);

				param.max_depth = 10;
				param.initial_depth = 5;
				param.stop_abs_diff = 1e-4;
				param.stop_related_diff = 0.02;
				param.epsilon = 0.002;
				param.result_eps = 1e-3;


				param0.max_depth = 6;
				param0.initial_depth = 4;
				param0.result_eps = 0.01;
			}
		

			learner.learn(contactspace_samples, contactspace.active_data_dim());

			SVMEvaluator evaluator(learner);

			learner.save("model.txt");

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;


			std::vector<DataVector> samples1, samples2, samples3;
			std::vector<DataVector> samples1_kcentriods, samples2_kcentriods, samples3_kcentriods;
			
			
			sample_decision_boundary_interpolation(learner, samples1, 50);
			samples1 = filter(evaluator, samples1, 1);
			samples1_kcentriods = sampleSelectionKCentroids(samples1, 100, 10);

			//std::cout << samples1_kcentriods.size() << std::endl;
			//
			//for(std::size_t i = 0; i < samples1_kcentriods.size(); ++i)
			//{
			//	double value = evaluator.evaluate(samples1_kcentriods[i]);
			//	std::cout << value << " ";
			//}
			//std::cout << std::endl;


			//std::ofstream boundary_sample_file("boundary_sample.txt");
			//for(std::size_t i = 0; i < samples1_kcentriods.size(); ++i)
			//{
			//	boundary_sample_file << checkStatus(contactspace.collider, contactspace.data_dim(), samples1_kcentriods[i]) << " ";
			//	for(std::size_t j = 0; j < samples1_kcentriods[i].dim(); ++j)
			//	{
			//		boundary_sample_file << samples1_kcentriods[i][j] << " ";
			//	}
			//	boundary_sample_file << std::endl;
			//}
			//boundary_sample_file.close();

			
			sample_decision_boundary_hierarchial_tree<SVMDistanceToDecisionBoundary_OptimizationGradient>(param0, learner, samples2);
			samples2_kcentriods = sampleSelectionKCentroids(samples2, 100, 10);

			for(std::size_t i = 0; i < samples2_kcentriods.size(); ++i)
			{
				double value = evaluator.evaluate(samples2_kcentriods[i]);
				std::cout << value << " ";
			}
			std::cout << std::endl;

			std::ofstream boundary_sample_file("boundary_sample.txt");
			for(std::size_t i = 0; i < samples2_kcentriods.size(); ++i)
			{
				boundary_sample_file << checkStatus(contactspace.collider, contactspace.data_dim(), samples2_kcentriods[i]) << " ";
				for(std::size_t j = 0; j < samples2_kcentriods[i].dim(); ++j)
				{
					boundary_sample_file << samples2_kcentriods[i][j] << " ";
				}
				boundary_sample_file << std::endl;
			}
			boundary_sample_file.close();
			




			sample_decision_boundary_hierarchial_tree_E<SVMEvaluator>(param, learner, samples3);
			samples3_kcentriods = sampleSelectionKCentroids(samples3, 100, 10);

			//std::ofstream boundary_sample_file("boundary_sample.txt");
			//for(std::size_t i = 0; i < samples3_kcentriods.size(); ++i)
			//{
			//	boundary_sample_file << checkStatus(contactspace.collider, contactspace.data_dim(), samples3_kcentriods[i]) << " ";
			//	for(std::size_t j = 0; j < samples3_kcentriods[i].dim(); ++j)
			//	{
			//		boundary_sample_file << samples3_kcentriods[i][j] << " ";
			//	}
			//	boundary_sample_file << std::endl;
			//}
			//boundary_sample_file.close();
		}

	}
}

void main()
{
	APDL::test_svm_boundary_sampler();
}


