#ifndef SPATIAL_TREE_H
#define SPATIAL_TREE_H

#include "data_vector.h"

namespace APDL
{
	struct SpatialTreeParam
	{
		std::size_t max_depth;
		std::size_t initial_depth;
		double stop_abs_diff;
		double stop_related_diff;
		double epsilon;
		double result_eps;

		SpatialTreeParam()
		{
			max_depth = 10;
			initial_depth = 3;
			stop_abs_diff = 0;
			stop_related_diff = 0.1;
			epsilon = 0;
			result_eps = 0;
		}
	};

	struct SpatialTreeNode
	{
		SpatialTreeNode(std::size_t dim_) : dim(dim_), center(dim_), scale(dim_)
		{
			parent = NULL;
			children = NULL;
			num_of_children = 0;
		}

		std::size_t getNumOfChildren() const
		{
			return num_of_children;
		}

		double getRadius() const 
		{
			double r = 0;
			for(std::size_t i = 0; i < dim; ++i)
				r += scale[i] * scale[i];
			return std::sqrt(r);
		}

		void initChildren() 
		{
			best_split_dir = getBestSplitDirection();
			bool is_binary = isBinarySplit(best_split_dir);
			if(is_binary)
			{
				num_of_children = 2;
				children = new SpatialTreeNode* [num_of_children];
			}
			else
			{
				num_of_children = 1 << dim;
				children = new SpatialTreeNode* [num_of_children];
			}
		}

		std::size_t dim;
		DataVector center;
		DataVector scale; // half length of each dimension
		SpatialTreeNode* parent;
		SpatialTreeNode** children;
		std::size_t num_of_children;

		std::size_t best_split_dir;

		double dist_to_decision_boundary;

	private:

		std::size_t getBestSplitDirection() const
		{
			std::size_t best_dir = 0;
			for(std::size_t i = 1; i < dim; ++i)
			{
				if(scale[i] > scale[best_dir]) best_dir = i;
			}

			return best_dir;
		}

		bool isBinarySplit(std::size_t best_dir) const
		{
			bool is_binary = false;
			for(std::size_t i = 0; i < dim; ++i)
			{
				if(scale[i] <= 0.5 * scale[best_dir]) 
				{
					is_binary = true;
					break;
				}
			}

			return is_binary;
		}
	};


	// general kdtree
	
	template<typename Learner, typename DistanceOracle>
	class SpatialTree
	{
	public:
		SpatialTree(const DataVector& lower, const DataVector& upper, const SpatialTreeParam& param, const Learner& learner_) : oracle(learner_)
		{
			dim = learner_.feature_dim;
			max_depth = param.max_depth;
			initial_depth = param.initial_depth;
			stop_abs_diff_percentage = param.stop_abs_diff;
			stop_related_diff_percentage = param.stop_related_diff;
			epsilon = param.epsilon;
			result_eps = param.result_eps;

			root = new SpatialTreeNode(dim);
			root->center = 0.5 * (upper + lower);
			root->scale = 0.5 * (upper - lower);

			root->parent = NULL;

			buildRecurse(root, 0);
		}

		~SpatialTree()
		{
			deleteRecurse(root);
			delete root;
		}

		void deleteRecurse(SpatialTreeNode* root)
		{
			if(root->num_of_children == 0)
			{
				return;
			}

			for(std::size_t i = 0; i < root->num_of_children; ++i)
			{
				deleteRecurse(root->children[i]);
				delete root->children[i];
			}

			delete [] root->children;
		}

		void collectBoundarySamples(std::vector<DataVector>& samples) const
		{
			collectRecurse(root, samples);
		}

		void collectRecurse(SpatialTreeNode* root, std::vector<DataVector>& samples) const
		{
			if(root->num_of_children > 0)
			{
				for(std::size_t i = 0; i < root->num_of_children; ++i)
					collectRecurse(root->children[i], samples);
			}
			else
			{
				if(root->dist_to_decision_boundary < result_eps)
					samples.push_back(root->center);
			}
		}

		void buildRecurse(SpatialTreeNode* root, std::size_t depth)
		{
			root->dist_to_decision_boundary = oracle.distance(root->center);

			if((root->dist_to_decision_boundary > root->getRadius() && depth > initial_depth) || depth > max_depth)
			{
				return;
			}
			else
			{
				root->initChildren();

				if(root->num_of_children == 2) // binary
				{
					for(std::size_t i = 0; i < 2; ++i)
					{
						root->children[i] = new SpatialTreeNode(dim);
						root->children[i]->parent = root;
					}

					DataVector center1(dim), center2(dim), scale(dim);
					for(std::size_t i = 0; i < dim; ++i)
					{
						if(i != root->best_split_dir)
						{
							center1[i] = root->center[i];
							center2[i] = root->center[i];
							scale[i] = root->scale[i];
						}
						else
						{
							center1[i] = root->center[i] + 0.5 * root->scale[i];
							center2[i] = root->center[i] - 0.5 * root->scale[i];
							scale[i] = 0.5 * root->scale[i];
						}
					}

					root->children[0]->center = center1; 
					root->children[1]->center = center2;
					root->children[0]->scale = scale;
					root->children[1]->scale = scale;

				}
				else
				{
					bool* bin = new bool[dim];
					for(std::size_t i = 0; i < root->num_of_children; ++i)
					{
						root->children[i] = new SpatialTreeNode(dim);
						root->children[i]->parent = root;
						intToBin(i, bin);
						DataVector center(root->center);
						for(std::size_t j = 0; j < dim; ++j)
						{
							if(bin[j]) center[j] += 0.5 * root->scale[j];
							else center[j] -= 0.5 * root->scale[j];
						}

						root->children[i]->center = center;
						root->children[i]->scale = 0.5 * root->scale;
					}

					delete [] bin;
				}

			
				for(std::size_t i = 0; i < root->num_of_children; ++i)
					buildRecurse(root->children[i], depth + 1);
			}
		}

		DistanceOracle oracle;

		SpatialTreeNode* root;

	private:

		std::size_t dim;

		void intToBin(std::size_t id, bool* bin) const
		{
			std::size_t v = id;
			for(std::size_t i = 0; i < dim; ++i)
			{
				if(v % 2) bin[i] = 1;
				else bin[i] = 0;
				v >>= 1;
			}
		}

		std::size_t max_depth;
		std::size_t initial_depth;
		double stop_abs_diff_percentage;
		double stop_related_diff_percentage;
		double epsilon;
		double result_eps;
	};



	struct SpatialTreeNodeE
	{
		SpatialTreeNodeE(std::size_t dim_) : dim(dim_), center(dim_), scale(dim_)
		{
			parent = NULL;
			children = NULL;
			num_of_children = 0;
			homogeneous = true;
		}

		std::size_t getNumOfChildren() const
		{
			return num_of_children;
		}

		void initChildren() 
		{
			best_split_dir = getBestSplitDirection();
			bool is_binary = isBinarySplit(best_split_dir);
			if(is_binary)
			{
				num_of_children = 2;
				children = new SpatialTreeNodeE* [num_of_children];
			}
			else
			{
				num_of_children = (1 << dim);
				children = new SpatialTreeNodeE* [num_of_children];
			}
		}

		std::size_t dim;
		DataVector center;
		DataVector scale; // half length of each dimension
		SpatialTreeNodeE* parent;
		SpatialTreeNodeE** children;
		std::size_t num_of_children;

		std::size_t best_split_dir;

		bool homogeneous;
		double value;

	private:

		std::size_t getBestSplitDirection() const
		{
			std::size_t best_dir = 0;
			for(std::size_t i = 1; i < dim; ++i)
			{
				if(scale[i] > scale[best_dir]) best_dir = i;
			}

			return best_dir;
		}

		bool isBinarySplit(std::size_t best_dir) const
		{
			bool is_binary = false;
			for(std::size_t i = 0; i < dim; ++i)
			{
				if(scale[i] <= 0.5 * scale[best_dir]) 
				{
					is_binary = true;
					break;
				}
			}

			return is_binary;
		}
	};

	struct SpatialTreeEParam
	{
		std::size_t max_depth;
		std::size_t initial_depth;
		double stop_abs_diff;
		double stop_related_diff;
		double epsilon;
		double result_eps;

		SpatialTreeEParam()
		{
			max_depth = 10;
			initial_depth = 3;
			stop_abs_diff = 0;
			stop_related_diff = 0.1;
			epsilon = 0;
			result_eps = 0;
		}
	};

	template<typename Learner, typename EvaluateOracle>
	class SpatialTreeE
	{
	public:
		SpatialTreeE(const DataVector& lower, const DataVector& upper, const SpatialTreeEParam& param, const Learner& learner_) : oracle(learner_)
		{
			dim = learner_.feature_dim;
			max_depth = param.max_depth;
			initial_depth = param.initial_depth;
			stop_abs_diff_percentage = param.stop_abs_diff;
			stop_related_diff_percentage = param.stop_related_diff;
			epsilon = param.epsilon;
			result_eps = param.result_eps;

			root = new SpatialTreeNodeE(dim);
			root->center = 0.5 * (upper + lower);
			root->scale = 0.5 * (upper - lower);

			root->parent = NULL;

			buildRecurse(root, 0);
		}                                                     

		~SpatialTreeE()
		{
			deleteRecurse(root);
			delete root;
		}

		void deleteRecurse(SpatialTreeNodeE* root)
		{
			if(root->num_of_children == 0)
			{
				return;
			}

			for(std::size_t i = 0; i < root->num_of_children; ++i)
			{
				deleteRecurse(root->children[i]);
				delete root->children[i];
			}

			delete [] root->children;
		}

		void collectBoundarySamples(std::vector<DataVector>& samples) const
		{
			collectRecurse(root, 0, samples);
		}

		void collectRecurse(SpatialTreeNodeE* root, std::size_t depth, std::vector<DataVector>& samples) const
		{
			if(root->num_of_children > 0) 
			{
				for(std::size_t i = 0; i < root->num_of_children; ++i)
					collectRecurse(root->children[i], depth + 1, samples);
			}
			else
			{
				if(!root->homogeneous || abs(root->value) < result_eps)
					samples.push_back(root->center);
			}
		}

		void buildRecurse(SpatialTreeNodeE* root, std::size_t depth)
		{
			bool homogeneous = true;
			double center_v = oracle.evaluate(root->center);
			root->value = center_v;

			bool sign  = (center_v > 0);
			
			bool* bin = new bool[dim];
			std::size_t num_of_corners = (1 << dim);
			for(std::size_t i = 0; i < num_of_corners; ++i)
			{
				intToBin(i, bin);
				DataVector corner(root->center);
				for(std::size_t j = 0; j < dim; ++j)
				{
					if(bin[j]) corner[j] += root->scale[j];
					else corner[j] -= root->scale[j];
				}

				double corner_v = oracle.evaluate(corner);

				bool corner_sign;
				if(sign)
				{
					if(corner_v < epsilon) corner_sign = false; 
					else corner_sign = true;
				}
				else
				{
					if(corner_v > -epsilon) corner_sign = true;
					else corner_sign = false;
				}

				if(corner_sign != sign) { homogeneous = false; break; }
			}


			delete [] bin;

			root->homogeneous = homogeneous;

			if((homogeneous && depth > initial_depth) || depth > max_depth)
			{
				return;
			}
			else
			{
				root->initChildren();

				if(root->num_of_children == 2) // binary
				{
					for(std::size_t i = 0; i < 2; ++i)
					{
						root->children[i] = new SpatialTreeNodeE(dim);
						root->children[i]->parent = root;
					}

					DataVector center1(dim), center2(dim), scale(dim);
					for(std::size_t i = 0; i < dim; ++i)
					{
						if(i != root->best_split_dir)
						{
							center1[i] = root->center[i];
							center2[i] = root->center[i];
							scale[i] = root->scale[i];
						}
						else
						{
							center1[i] = root->center[i] + 0.5 * root->scale[i];
							center2[i] = root->center[i] - 0.5 * root->scale[i];
							scale[i] = 0.5 * root->scale[i];
						}
					}

					root->children[0]->center = center1; 
					root->children[1]->center = center2;
					root->children[0]->scale = scale;
					root->children[1]->scale = scale;

				}
				else
				{
					bool* bin = new bool[dim];
					for(std::size_t i = 0; i < root->num_of_children; ++i)
					{
						root->children[i] = new SpatialTreeNodeE(dim);
						root->children[i]->parent = root;
						intToBin(i, bin);
						DataVector center(root->center);
						for(std::size_t j = 0; j < dim; ++j)
						{
							if(bin[j]) center[j] += 0.5 * root->scale[j];
							else center[j] -= 0.5 * root->scale[j];
						}

						root->children[i]->center = center;
						root->children[i]->scale = 0.5 * root->scale;
					}

					delete [] bin;
				}

				for(std::size_t i = 0; i < root->num_of_children; ++i)
					buildRecurse(root->children[i], depth + 1);
			}
		}

		EvaluateOracle oracle;

		SpatialTreeNodeE* root;

		std::size_t max_depth;
		std::size_t initial_depth;
		double stop_abs_diff_percentage;
		double stop_related_diff_percentage;
		double epsilon;
		double result_eps;


	private:

		std::size_t dim;

		void intToBin(std::size_t id, bool* bin) const
		{
			std::size_t v = id;
			for(std::size_t i = 0; i < dim; ++i)
			{
				if(v % 2) bin[i] = 1;
				else bin[i] = 0;
				v >>= 1;
			}
		}
	};
}

#endif