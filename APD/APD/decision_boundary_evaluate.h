#ifndef DECISION_BOUNDARY_EVALUATE_H
#define DECISION_BOUNDARY_EVALUATE_H

#include "contact_space_learning.h"
#include "learning/svm.h"

namespace APDL
{
	struct SVMEvaluator
	{
		SVMEvaluator(const SVMLearner& learner_) : learner(learner_)
		{
			node = new svm_node[learner.feature_dim + 1];
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].index = i + 1;
			node[learner.feature_dim].index = -1;
		}

		~SVMEvaluator()
		{
			delete [] node;
		}

		double evaluate(const DataVector& v) const
		{
			for(std::size_t i = 0; i < learner.feature_dim; ++i)
				node[i].value = v[i];

			return svm_predict_values_twoclass(learner.model, node);
		}

		const SVMLearner& learner;

		mutable svm_node* node;
	};

	struct MulticonlitronEvaluator
	{
		MulticonlitronEvaluator(const MulticonlitronLearner& learner_) : learner(learner_)
		{

		}

		double evaluate(const DataVector& v) const
		{
			return learner.model.evaluate2(v);
		}



		const MulticonlitronLearner& learner;
	};
}

#endif