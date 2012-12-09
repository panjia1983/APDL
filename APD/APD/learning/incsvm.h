#ifndef INCSVM_H
#define INCSVM_H 1

#include <vector>
#include <time.h>
#include <map>
#include <cstdio>
#include <cstdlib>

void my_assert(int assertion, char * msg);

using namespace std;

class stopwatch
{
public:
	stopwatch() {}
	void reset() {}
	~stopwatch() {}
	double get_time() {}
};

class Matrix
{
public:
	double **array;
	double *p;
	int size;
	int* row_size;
	int m, n;
	
	Matrix()
	{
	}
	
	Matrix(int _m, int _n)
	{
		initialize(_m, _n);
	}
	
	void initialize(int _m, int _n)
	{
		size = 0;
		m = _m;
		n = _n;
		p = new double[m * n];
		row_size = new int[m];
		array = new double*[m];
		for(int i = 0; i < m; i++)
		{
			array[i] = &p[i * n];
			row_size[i] = 0;
			for(int j = 0; j < n; j++)
				array[i][j] = 0;
		}
	}
	
	~Matrix()
	{
		//printf("Deleting matrix\n");
		delete[] row_size;
		delete[] p;
		delete[] array;
	}
	
	void del(int i, int j)
	{
		for(int k = j; k < row_size[i] - 1; k++)
			array[i][k] = array[i][k + 1];
		array[i][row_size[i] - 1] = 0;
		row_size[i]--;
	}
	
	void del_col(int j)
	{
		for(int k = 0; k < size; k++)
			if(j < row_size[k])
				del(k, j);
	}
	
	void del_row(int row)
	{
		for(int i = row; i < size - 1; i++)
		{
			for(int j = 0; j < row_size[i]; j++)
				array[i][j] = array[i + 1][j];
			row_size[i] = row_size[i + 1];
		}
		size--;
	}
	
	void print_row(int ind)
	{
		printf("Row %d: ", ind);
		for(int i = 0; i < row_size[ind]; i++)
			printf("%lf ", array[ind][i]);
		printf("\n");
	}
	
	void push_back(int i, double val)
	{
		array[i][row_size[i]] = val;
		row_size[i]++;
	}
	
	void push_back(double val)
	{
		push_back(size - 1, val);
	}
	
	void clear(int i)
	{
		row_size[i] = 0;
	}
	
	void addrow()
	{
		row_size[size] = 0;
		size++;
		if(size >= m)
		{
			//TODO: add code to automatically increase matrix size when this happens
			fprintf(stderr, "Matrix (%d, %d) out of space %d\n", m, n, size);
			fflush(stderr);
			exit(-3);
		}
	}
	
	double* get_row(int row)
	{
		return array[row];
	}
	
	int get_row_size(int row)
	{
		return row_size[row];
	}
	
};



class SparseMatrix
{
	double **array;
public:
	double *p;
	int size;
	int * row_size;
	int m, n;
	map<int, int> index_to_row;
	map<int, int> row_to_index;
	
	SparseMatrix()
	{
	}
	
	SparseMatrix(int _m, int _n)
	{
		initialize(_m, _n);
	}
	
	void initialize(int _m, int _n)
	{
		size = 0;
		m = _m;
		n = _n;
		p = new double[m * n];
		row_size = new int[m];
		array = new double*[m];
		for(int i = 0; i < m; i++)
		{
			array[i] = &p[i * n];
			row_size[i] = 0;
			for(int j = 0; j < n; j++)
				array[i][j] = 0;
		}
	}
	
	~SparseMatrix()
	{
		//printf("Deleting matrix\n");
		delete[] row_size;
		delete[] p;
		delete[] array;
	}
	
	void del(int i, int j)
	{
		//my_assert(i<size && j<row_size[i],"Line 67: libincsvm");
		for(int k = j; k < row_size[i] - 1; k++)
			array[i][k] = array[i][k + 1];
		array[i][row_size[i] - 1] = 0;
		row_size[i]--;
	}
	
	void del_col(int j)
	{
		for(int k = 0; k < size; k++)
			if(j < row_size[k])
				del(k, j);
	}
	
	void del_row(int ind)
	{
		int row = index_to_row[ind];
		for(int i = row; i < size - 1; i++)
		{
			for(int j = 0; j < row_size[i]; j++)
			{
				array[i][j] = array[i + 1][j];
			}
			row_to_index[i] = row_to_index[i + 1];
			index_to_row[row_to_index[i]] = i;
			row_size[i] = row_size[i + 1];
		}
		row_to_index[size] = -1;
		index_to_row[ind] = -1;
		size--;
	}
	
	void print_row(int ind)
	{
		int row = index_to_row[ind];
		printf("Row %d: ", row);
		for(int i = 0; i < row_size[row]; i++)
			printf("%lf ", array[row][i]);
		printf("\n");
	}
	
	void push_back(int ind, double val)
	{
		int i = index_to_row[ind];
		array[i][row_size[i]] = val;
		row_size[i]++;
	}
	
	void push_back(double val)
	{
		push_back(size - 1, val);
	}
	
	void clear(int i)
	{
		del_row(i);
		//row_size[i] = 0;
	}
	
	void addrow(int ind)
	{
		row_size[size] = 0;
		index_to_row[ind] = size;
		row_to_index[size] = ind;
		size++;
		if(size >= m)
		{
			fprintf(stderr, "Matrix (%d, %d) out of space %d\n", m, n, size);
			fflush(stderr);
			exit(-3);
		}
	}
	
	double * get_row(int ind)
	{
		int row = index_to_row[ind];
		return array[row];
	}
	
	int get_row_size(int ind)
	{
		if(index_to_row.find(ind) != index_to_row.end() && index_to_row[ind] != -1)
			return row_size[index_to_row[ind]];
		return 0;
	}
};



class IncSVM
{

public:
	//alpha values
	vector<double> a;
	
	//classifier bias
	double b;
	
	vector<int> ind[4];
	
	virtual double svmeval(int indc, std::vector<double>* Qc);
	
	virtual int learn(int indc, int rflag);
	
	virtual int unlearn(int indc);
	
	virtual void initialize(int m, int n);
	
	virtual void adddata(int idx, int y_new, double C_new);
	
	virtual void setY(int idx, int y_new);
	
	std::vector<double>& get_model(double& b0);
	
	virtual ~IncSVM() {};
	
	double abs(double a);
	
	double get_w2();
	
private:
	Matrix _Rs;
	SparseMatrix _Q;
	vector<double> c;
	double deps;
	vector<double> g;
	double scale;
	double type;
	vector<double> y;
	int perturbations;
	int max_reserve_vectors;
	int kernel_evals;
	stopwatch watch;
	
	
	int find(int indc, vector<int>& indices);
	
	void move_ind(int cstatus, int nstatus, int index);
	
	int move_indr(int cstatus, int indc);
	
	void bookkeeping(int indss, int cstatus, int nstatus, int& indco, int& idx);
	
	void select_indices(vector<double>& array, vector<int>& indices, vector<double>& res);
	
	void updateRQ(int indc);
	
	void updateRQ(vector<double>& beta, vector<double>& tmp1, int indc);
	
	void minof(double * a, int size, double& res, int& ind);
	
	void min_delta(vector<int>& flags, vector<double> psi_initial, vector<double> psi_final, vector<double> psi_sens, double& min_d, int& k);
	
	void min_delta_acb(int indc, vector<double>& gamma, vector<double>& beta, int polc, int rflag, double& min_dacb, int& indss_, int& cstatus_, int& nstatus_);
	
};

void test_incsvm(double C, const std::string& kernel_file);





#endif
