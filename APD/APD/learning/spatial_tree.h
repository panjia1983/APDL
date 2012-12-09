#ifndef SPATIAL_TREE_H
#define SPATIAL_TREE_H

#include "data_vector.h"

namespace APDL
{
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
		SpatialTree(const DataVector& lower, const DataVector& upper, const Learner& learner_) : oracle(learner_)
		{
			dim = learner_.feature_dim;

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

		void collectBoundarySamples(double on_decision_threshold, std::vector<DataVector>& samples) const
		{
			collectRecurse(root, on_decision_threshold, samples);
		}

		void collectRecurse(SpatialTreeNode* root, double on_decision_threshold, std::vector<DataVector>& samples) const
		{
			if(root->num_of_children > 0)
			{
				for(std::size_t i = 0; i < root->num_of_children; ++i)
					collectRecurse(root->children[i], on_decision_threshold, samples);
			}
			else
			{
				if(root->dist_to_decision_boundary < on_decision_threshold)
				{
					samples.push_back(root->center);
				}
			}
		}

		void buildRecurse(SpatialTreeNode* root, std::size_t depth)
		{
			root->dist_to_decision_boundary = oracle.distance(root->center);
			//std::cout << depth << " ";
			//for(std::size_t i = 0; i < dim; ++i)
			//	std::cout << root->center[i] << " ";
			//for(std::size_t i = 0; i < dim; ++i)
			//	std::cout << root->scale[i] << " ";
			//std::cout << root->dist_to_decision_boundary << std::endl;
			if(root->dist_to_decision_boundary > root->getRadius() || depth > 15)
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
	};



	struct SpatialTreeNodeE
	{
		SpatialTreeNodeE(std::size_t dim_) : dim(dim_), center(dim_), scale(dim_)
		{
			parent = NULL;
			children = NULL;
			num_of_children = 0;
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


	template<typename Learner, typename EvaluateOracle>
	class SpatialTreeE
	{
	public:
		SpatialTreeE(const DataVector& lower, const DataVector& upper, const Learner& learner_) : oracle(learner_)
		{
			dim = learner_.feature_dim;

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

		void collectBoundarySamples(double on_decision_threshold, std::vector<DataVector>& samples) const
		{
			collectRecurse(root, on_decision_threshold, samples);
		}

		void collectRecurse(SpatialTreeNodeE* root, double on_decision_threshold, std::vector<DataVector>& samples) const
		{
			if(root->num_of_children > 0)
			{
				for(std::size_t i = 0; i < root->num_of_children; ++i)
					collectRecurse(root->children[i], on_decision_threshold, samples);
			}
			else
			{
				if(root->dist_to_decision_boundary < on_decision_threshold)
				{
					samples.push_back(root->center);
				}
			}
		}

		void buildRecurse(SpatialTreeNode* root, std::size_t depth)
		{
			bool homogeneous = true;
			bool sign  = oracle.evaluate(root->center);
			
			bool* bin = new bool[dim];
			std::size_t num_of_corners = 1 << dim;
			for(std::size_t i = 0; i < root->num_of_corners; ++i)
			{
				intToBin(i, bin);
				DataVector center(root->center);
				for(std::size_t j = 0; j < dim; ++j)
				{
					if(bin[j]) center[j] += root->scale[j];
					else center[j] -= root->scale[j];
				}

				bool corner_sign = oracle.evaluate(center);
				if(corner_sign != sign) homogeneous = false;
			}

			delete [] bin;


			if(homogeneous || depth > 15)
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

		EvaluateOracle oracle;

		SpatialTreeNodeE* root;

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