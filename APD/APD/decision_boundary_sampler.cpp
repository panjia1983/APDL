#include "decision_boundary_sampler.h"


namespace APDL
{
	void sample_decision_boundary_interpolation(const SVMLearner& learner,
												std::vector<DataVector>& samples)
	{
		const struct svm_model* model = learner.model;
		std::size_t dim = learner.feature_dim;
		double w_norm = sqrt(learner.hyperw_normsqr);

		std::size_t start[2];
		std::size_t number[2];
		start[0] = 0;
		number[0] = model->nSV[0];
		start[1] = start[0] + number[0];
		number[1] = model->nSV[1];



		svm_node* midpoint = new svm_node[dim + 1];
		for(std::size_t i = 0; i < dim; ++i) midpoint[i].index = i;
		midpoint[dim].index = -1;

		DataVector p(dim);

		for(std::size_t i = start[0]; i < start[0] + number[0]; ++i)
		{
			for(std::size_t j = start[1]; j < start[1] + number[1]; ++j)
			{
				svm_node* x = model->SV[i];
				svm_node* y = model->SV[j];

				feature_space_midpoint(model, x, y, midpoint);

				double f_value = svm_predict_values_twoclass(model, midpoint);

				// std::cout << f_value / w_norm << std::endl;

				if(std::abs(f_value) / w_norm > 0.015) continue;

				for(std::size_t k = 0; k < dim; ++k)
					p[k] = midpoint[k].value;

				samples.push_back(p);
			}
		}

		delete [] midpoint;
	}
}