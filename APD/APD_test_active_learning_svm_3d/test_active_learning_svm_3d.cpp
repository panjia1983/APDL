#include <APD/active_learning.h>
#include <APD/contact_space_learning.h>
#include <APD/minkowski_cspace.h>
#include <APD/decision_boundary_sampler.h>


namespace APDL
{
	void test_svm_learner_3d_cupspoon()
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


		// original learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 


			std::ofstream scaler_file("scaler_3d_cupspoon.txt");
			scaler_file << contactspace.getScaler() << std::endl;

			learner.learn(contactspace_samples, contactspace.active_data_dim());
			learner.save("model_3d_cupspoon.txt");

			std::cout << "model saved" << std::endl;

			std::cout << contactspace_samples.size() << ": " << empiricalErrorRatio(contactspace_samples, learner) << " " << errorRatioOnGrid(contactspace, learner, 20) << std::endl;
		}

		// active learner, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 9);
			aparam.debug = true;
			aparam.num_grid = 20;
			active_learning(contactspace, learner, decision_boundary_sampler, aparam);
		}

		// active learner2, scaled
		{
			SVMLearner learner;
			learner.setC(20);
			learner.setProbability(true);
			learner.setScaler(contactspace.getScaler());
			learner.setUseScaler(true);
			learner.setGamma(50); 

			SpatialTreeEParam param;
			param.max_depth = 4;
			param.initial_depth = 2;
			param.stop_abs_diff = 0.2;
			param.stop_related_diff = 0.1;
			param.epsilon = 0;
			param.result_eps = 0;

			SVMEvaluator evaluator(learner);
			FilterParam fparam;
			fparam.filter_threshold = 1;
			DecisionBoundaryHierarchialTreeESampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(param, fparam, learner);
			// DecisionBoundaryInterpolationSampler<SVMLearner, SVMEvaluator> decision_boundary_sampler(10, fparam, evaluator, learner);


			ActiveLearningParam aparam(1000, 500, 500, 18);
			aparam.debug = true;
			aparam.num_grid = 20;
			active_learning2(contactspace, learner, decision_boundary_sampler, aparam);
		}

		delete P;
		delete Q;
	}
}

void main()
{
	APDL::test_svm_learner_3d_cupspoon();
}