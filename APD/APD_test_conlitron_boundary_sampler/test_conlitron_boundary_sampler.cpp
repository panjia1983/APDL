#include <APD/minkowski_cspace.h>
#include <APD/contact_space_learning.h>
#include <APD/decision_boundary_sampler.h>
#include <APD/active_learning.h>

void* user_conlitron_model;
double* user_conlitron_data;

namespace APDL
{	
	void test_conlitron_boundary_sampler()
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

			DataVector w(2);
			w[0] = 1; w[1] = 1;
			MulticonlitronLearner learner(w, 0.01);
			learner.setScaler(contactspace.getScaler());
			learner.setDim(contactspace.active_data_dim());
			learner.learn(contactspace_samples, 2);
			SpatialTreeEParam param;

			param.max_depth = 10;
			param.initial_depth = 3;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;	
			param.result_eps = 0;


			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 100) << std::endl;

			learner.saveVisualizeData("conlitron_2d_vis.txt", contactspace.getScaler(), 100);
			MulticonlitronEvaluator evaluator(learner);

			std::vector<DataVector> samples1;
			sample_decision_boundary_hierarchial_tree_E<SVMLearner, MulticonlitronEvaluator>(param, learner, samples1);
			std::vector<DataVector> samples1_kcentriods = sampleSelectionKCentroids(samples1, 100, 10);
			std::cout << samples1.size() << std::endl;
			
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





			//std::vector<DataVector> samples1_kmeans = sampleSelectionKMeans(samples1, 100);
			//for(std::size_t i = 0; i < samples1_kmeans.size(); ++i)
			//{
			//	double value = evaluator.evaluate(samples1_kmeans[i]);
			//	std::cout << value << " ";
			//}
			//std::cout << std::endl;

			std::vector<DataVector> samples2;
			sample_decision_boundary_interpolation(learner, samples2, 50);
			samples2 = filter(evaluator, samples2, 0.01);
			std::cout << samples2.size() << std::endl;

			std::vector<DataVector> samples2_kcentriods = samples2; // sampleSelectionKCentroids(samples2, 100, 10);
			//for(std::size_t i = 0; i < samples2_kcentriods.size(); ++i)
			//{
			//	double value = evaluator.evaluate(samples2_kcentriods[i]);
			//	std::cout << value << " ";
			//}
			//std::cout << std::endl;


			//std::ofstream boundary_sample_file("boundary_sample.txt");
			//for(std::size_t i = 0; i < samples2_kcentriods.size(); ++i)
			//{
			//	boundary_sample_file << checkStatus(contactspace.collider, contactspace.data_dim(), samples2_kcentriods[i]) << " ";
			//	for(std::size_t j = 0; j < samples2_kcentriods[i].dim(); ++j)
			//	{
			//		boundary_sample_file << samples2_kcentriods[i][j] << " ";
			//	}
			//	boundary_sample_file << std::endl;
			//}
			//boundary_sample_file.close();
		}
	}
}

void main()
{
	APDL::test_conlitron_boundary_sampler();
}