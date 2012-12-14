#include <APD/learning/incsvm.h>

namespace APDL
{
	void test_incremental_svm()
	{
		test_incsvm(1.0, "../data/data_kernel.txt");
		test_incsvm(1.0, "../data/diabetes_kernel.txt");
	}

}

void main()
{
	APDL::test_incremental_svm();
}