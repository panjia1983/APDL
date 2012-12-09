#include <iostream>
#include "incsvm.h"
#include <vector>
#include <fstream>

using namespace std;

vector<vector<double> > kernel_values;
vector<int> labels;
double C = 100;

void load_kernel_from_file(const char * filename)
{
  ifstream ifs(filename, ios_base::in);
  if(!ifs.is_open()) 
  {
    cout << "Unable to read file " << filename << endl;
    exit(1);
  }

  int m = -1, n = -1;
  ifs >> m >> n;
  if(!ifs.good() || m < 0 || n < 0) 
  {
    cout << "Incorrect row or column size: " << m <<", "<< n <<endl;
    exit(1);
  }

  kernel_values.resize(m);
  for(int i=0;i<m;i++)
  {
    kernel_values[i].resize(n);
    int label;
    ifs >> label;
    labels.push_back(label);

    double val;
    for(int j=0;j<n;j++)
	{
      ifs >> val;
      if(!ifs.good())
	  {
		cout << "Error reading row "<< m+1<<", column "<< j+2 << endl;
		exit(1);	
      }
      kernel_values[i][j] = val;
    }
  } 
  ifs.close();
}

// implement kernel function that returns kernel values between examples with
// indices i, j
double kernel(int i, int j, void *kparam)
{
  return kernel_values[i][j];
}

void process_instance(IncSVM& svm, int idx, int label)
{
  svm.setY(idx,label);
  svm.learn(idx,1);
}

void unlearn_instance(IncSVM& svm, int idx, int label)
{
  svm.unlearn(idx);
  svm.setY(idx,label);
}

double predict_value(IncSVM& svm, int idx)
{
    return svm.svmeval(idx,NULL);
}

void test_incsvm(double C, const std::string& kernel_file)
{
	load_kernel_from_file(kernel_file.c_str());

	IncSVM svm;
  
	/* intialize(m, n) where m is the upper limit on number of training examples and n is the upper limit on number of support vectors that you expect
	* (never got around to removing this restriction, please see L122 in libincsvm.h for details)
	*/
	svm.initialize(kernel_values.size(),1000);

	for(int i=0;i<kernel_values.size();i++)
	{
		// register index i, any dummy label, svm C parameter 
		svm.adddata(i,labels[i],C);
	}

	//always call adddata whenever you're adding new examples (train, unlabeled or test) into the classifier
	//make sure each example has a unique index 
	//the IncSVM class works directly through the example index values and you need to keep track of loading the feature
	//vector and computing the kernel values etc

	//training loop
	for(int i = 0; i< kernel_values.size(); i++)
	{
		process_instance(svm, i, labels[i]);
	}

	int num_correct = 0;
	//unlearning loop
	for(int i = 0; i< kernel_values.size(); i++)
	{
		//leave one out
		unlearn_instance(svm, i, labels[i]);

		//predict label
		double out = predict_value(svm, i);
		if(out * labels[i] > 0)
			num_correct++;

		//add example back
		process_instance(svm, i, labels[i]);
	}
}
