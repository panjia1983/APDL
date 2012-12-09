#include "incsvm.h"
#include <vector>
#include <limits>
#include <cmath>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


extern "C" {
#include <cblas.h>
}

using namespace std;

#define MARGIN 0
#define ERROR  1
#define RESERVE 2
#define UNLEARNED 3


double kernel(int i, int j, void * params);

void print_vector(vector<double>& a)
{
	for(unsigned int i = 0; i < a.size(); i++)
		//if(a[i]!=0)
		printf("%g ", a[i]);
	printf("\n");
}


void print_vector(vector<int>& a)
{
	for(unsigned int i = 0; i < a.size(); i++)
		printf("%d ", a[i]);
	printf("\n");
}

void print_matrix(vector<vector<double>*>& mat)
{
	for(unsigned int i = 0; i < mat.size(); i++)
		print_vector(*mat[i]);
}

int IncSVM::find(int indc, vector<int>& indices)
{
	for(unsigned int i = 0; i < indices.size(); i++)
		if(indices[i] == indc)
			return i;
	return -1;
}

double IncSVM::svmeval(int indc, vector<double>* Qc)
{
	double fx = b;
	int sets[] = {MARGIN, ERROR, RESERVE}; //,UNLEARNED};
	for(int j = 0; j < 3; j++)
	{
		int set = sets[j];
		for(unsigned int i = 0; i < ind[set].size(); i++)
		{
			int index = ind[set][i];
			double k = kernel(indc, index, NULL);
			if(Qc != NULL)
			{
				Qc[set].push_back(k * y[index]*y[indc]);
			}
			if(set == MARGIN || set == ERROR || a[index] > 0)
				fx = fx + a[index] * y[index] * k;
		}
	}
	return fx;
}


void IncSVM::move_ind(int cstatus, int nstatus, int index)
{
	int indc = ind[cstatus][index];
	ind[cstatus].erase(ind[cstatus].begin() + index);
	ind[nstatus].push_back(indc);
}

int IncSVM::move_indr(int cstatus, int indc)
{
	int removed_i = -1;
	move_ind(cstatus, RESERVE, indc);
	// implement max reserve vectors thing
	// (optional)
	return removed_i;
}

void IncSVM::bookkeeping(int indss, int cstatus, int nstatus, int& indco, int& idx)
{
	indco = -1;
	idx = -1;
	if(cstatus == nstatus)
		return;
	switch(nstatus)
	{
	case RESERVE:
		a[indss] = 0;
		break;
	case ERROR:
		a[indss] = c[indss];
		break;
	}
	
	int index = find(indss, ind[cstatus]);
	if(cstatus == MARGIN)
		indco = index + 1;
		
	if(nstatus == RESERVE)
		idx = move_indr(cstatus, index);
	else
		move_ind(cstatus, nstatus, index);
}

void vector_plus_vector(int start1, int c1, vector<double>& v1, int start2, int c2, vector<double>& v2, int start_res, vector<double>& result, int size)
{
	for(int i = 0; i < size; i++)
		result[start_res + i] = c1 * v1[start1 + i] + c2 * v2[start2 + i];
}

void vector_plus_vector(int start1, int c1, vector<double>& v1, int start2, int c2, vector<double>& v2, vector<int>& indices, vector<double>& result)
{
	for(unsigned int i = 0; i < indices.size(); i++)
		result[indices[i]] = c1 * v1[start1 + i] + c2 * v2[start2 + i];
}


void vector_plus_vector(vector<int>& ind1, int c1, vector<double>& v1, int start2, int c2, vector<double>& v2, vector<int>& indres, vector<double>& result)
{
	for(unsigned int i = 0; i < ind1.size(); i++)
		result[indres[i]] = c1 * v1[ind1[i]] + c2 * v2[start2 + i];
}

double vector_star_vector(vector<double>& v1, vector<double>& v2)
{
	double res = 0;
	for(unsigned int i = 0; i < v1.size(); i++)
		res += v1[i] * v2[i];
	return res;
}

double vector_star_vector(double v1[], vector<double>& v2)
{
	double res = 0;
	for(unsigned int i = 0; i < v2.size(); i++)
		res += v1[i] * v2[i];
	return res;
}


double vector_star_vector(double v1[], double v2[], int size)
{
	double res = 0;
	for(int i = 0; i < size; i++)
		res += v1[i] * v2[i];
	return res;
}


void matrix_star_vector(int c, vector<vector<double>*>& A, vector<double>& x, vector<double>& res)
{
	res.clear();
	for(unsigned int i = 0; i < A.size(); i++)
		res.push_back(c * vector_star_vector(*(A[i]), x));
}

void matrix_star_vector(int c, vector<vector<double>*>& A, double x[], vector<double>& res)
{
	res.clear();
	for(unsigned int i = 0; i < A.size(); i++)
		res.push_back(c * vector_star_vector(x, *(A[i])));
}

void matrix_star_vector(int c, Matrix& A, vector<double>& x, vector<double>& res)
{
	res.clear();
	for(int i = 0; i < A.size; i++)
		res.push_back(c * vector_star_vector(A.array[i], x));
}

void matrix_star_vector(int c, Matrix& A, double x[], vector<double>& res)
{
	res.clear();
	for(int i = 0; i < A.size; i++)
		res.push_back(c * vector_star_vector(x, A.array[i], A.row_size[i]));
}



void matrix_star_vector(int c, vector<double>** A, int size, vector<double>& x, vector<double>& res)
{
	res.clear();
	for(int i = 0; i < size; i++)
		res.push_back(c * vector_star_vector(*(A[i]), x));
}




void IncSVM::updateRQ(int indc)
{
	int rows = _Rs.size;
	stopwatch sw;
	if(rows > 2)
	{
		double val = _Rs.array[indc][indc];
		for(int i = 0; i < rows; i++)
			for(int j = 0; j < rows; j++)
				if(i != indc && j != indc)
					_Rs.array[i][j] = _Rs.array[i][j] -  _Rs.array[i][indc] * _Rs.array[indc][j] / val;
		_Rs.del_col(indc);
		_Rs.del_row(indc);
	}
	else
	{
		_Rs.row_size[0] = 0;
		_Rs.row_size[1] = 0;
		_Rs.size = 0;
	}
	_Q.del_col(indc);
}

void IncSVM::updateRQ(vector<double>& beta, vector<double>& tmp1, int indc)
{
	stopwatch watch;
	double gamma = tmp1[0];
	int rows = _Rs.size;
	if(gamma < deps)
		gamma = deps;
	if(rows > 1)
	{
		for(int i = 0; i < _Rs.size; i++)
		{
			for(int j = 0; j < _Rs.row_size[i]; j++)
				_Rs.array[i][j] = _Rs.array[i][j] + beta[i] * beta[j] / gamma;
			_Rs.push_back(i, beta[i] / gamma);
		}
		_Rs.addrow();
		for(unsigned int i = 0; i < beta.size(); i++)
			_Rs.push_back(beta[i] / gamma);
		_Rs.push_back(1 / gamma);
	}
	else
	{
		_Rs.addrow();
		_Rs.push_back(-(kernel(indc, indc, NULL) + deps));
		_Rs.push_back(y[indc]);
		_Rs.addrow();
		_Rs.push_back(y[indc]);
		_Rs.push_back(0);
	}
	
	for(unsigned int i = 0; i < y.size(); i++)
	{
		double val = 0;
		if(i == (unsigned int)indc)
			val = deps;
		if(_Q.get_row_size(i) > 0)
			_Q.push_back(i, y[indc]*y[i]*kernel(indc, i, NULL) + val);
	}
}

void IncSVM::minof(double * a, int size, double& res, int& ind)
{
	res = std::numeric_limits<double>::infinity();
	ind = 0;
	for(int i = 0; i < size; i++)
		if(a[i] < res)
		{
			res = a[i];
			ind = i;
		}
}

double IncSVM::abs(double val)
{
	if(val < 0)
		val = -val;
	return val;
}


void IncSVM::min_delta(vector<int>& flags, vector<double> psi_initial, vector<double> psi_final, vector<double> psi_sens, double& min_d, int& k)
{
	vector<int> ind;
	int idx;
	for(unsigned int i = 0; i < flags.size(); i++)
		if(flags[i])
			ind.push_back(i);
	if(ind.size() > 0)
	{
		double* deltas = new double[ind.size()];
		for(unsigned int i = 0; i < ind.size(); i++)
		{
			deltas[i] = (psi_final[ind[i]] - psi_initial[ind[i]]) / psi_sens[ind[i]];
		}
		minof(deltas, ind.size(), min_d, idx);
		k = ind[idx];
		double max_sens = -std::numeric_limits<double>::infinity();
		for(unsigned int i = 0; i < ind.size(); i++)
		{
			if(deltas[i] == min_d && abs((double)psi_sens[ind[i]]) > max_sens)
			{
				max_sens = abs((double)psi_sens[ind[i]]);
				k = ind[i];
			}
		}
		delete [] deltas;
	}
	else
	{
		min_d = std::numeric_limits<double>::infinity();
		k = -1;
	}
}




void IncSVM::select_indices(vector<double>& array, vector<int>& indices, vector<double>& res)
{
	for(unsigned int i = 0; i < indices.size(); i++)
	{
		res.push_back(array[indices[i]]);
	}
}

void IncSVM::min_delta_acb(int indc, vector<double>& gamma, vector<double>& beta, int polc, int rflag, double& min_dacb, int& indss_, int& cstatus_, int& nstatus_)
{
	vector<int> indss(6, 0);
	indss[0] = indss[1] = indss[2] = indc;
	vector<int> cstatus(6, 0);
	vector<int> nstatus(6, 0);
	double delta_m, delta_e, delta_r;
	if(polc == 1)
	{
		delta_m = -g[indc] / gamma[indc];
		nstatus[0] = MARGIN;
		if(beta.size() > 1)
		{
			delta_e = c[indc] - a[indc];
			nstatus[1] = ERROR;
		}
		else
			delta_e = std::numeric_limits<double>::infinity();
		delta_r = std::numeric_limits<double>::infinity();
	}
	else
	{
		if(g[indc] > 0)
		{
			delta_m = g[indc] / gamma[indc];
			nstatus[0] = MARGIN;
		}
		else
			delta_m = std::numeric_limits<double>::infinity();
			
		if(a[indc] <= c[indc])
		{
			delta_e = std::numeric_limits<double>::infinity();
			delta_r = a[indc];
			if(g[indc] > 0)
				nstatus[2] = RESERVE;
			else
				nstatus[2] = UNLEARNED;
		}
		else
		{
			delta_e  = a[indc] - c[indc];
			delta_r = std::numeric_limits<double>::infinity();
			nstatus[1] = ERROR;
		}
	}
	
	
	// change in a_c or b that causes a margin vector to an error or reserve vector
	vector<double> beta_s;
	vector<int> flags;
	double delta_mer;
	int idx;
	if(beta.size() > 1)
	{
		for(unsigned int i = 1; i < beta.size(); i++)
		{
			beta_s.push_back(polc * beta[i]);
			if(abs(polc * beta[i]) > 0)
				flags.push_back(1);
			else
				flags.push_back(0);
		}
		
		vector<double> a_tmp;
		vector<double> c_tmp;
		
		select_indices(a, ind[MARGIN], a_tmp);
		select_indices(c, ind[MARGIN], c_tmp);
		for(unsigned int i = 0; i < c_tmp.size(); i++)
			c_tmp[i] *= (beta_s[i] > 0);
		min_delta(flags, a_tmp, c_tmp, beta_s, delta_mer, idx);
		if(delta_mer < std::numeric_limits<double>::infinity())
		{
			indss[3] = ind[MARGIN][idx];
			cstatus[3] = MARGIN;
			nstatus[3] = ERROR * (beta_s[idx] > 0) + RESERVE * (beta_s[idx] < 0);
		}
	}
	else
		delta_mer = std::numeric_limits<double>::infinity();
		
	// change in a_c or b that causes an error vector to change to a margin vector
	vector<double> gamma_e;
	flags.clear();
	double delta_em;
	select_indices(gamma, ind[ERROR], gamma_e);
	for(unsigned int i = 0; i < gamma_e.size(); i++)
	{
		gamma_e[i] *= polc;
		flags.push_back(gamma_e[i] > 0);
	}
	vector<double> g_tmp;
	vector<double> z_tmp(ind[ERROR].size(), 0);
	select_indices(g, ind[ERROR], g_tmp);
	min_delta(flags, g_tmp, z_tmp, gamma_e, delta_em, idx);
	if(delta_em < std::numeric_limits<double>::infinity())
	{
		indss[4] = ind[ERROR][idx];
		cstatus[4] = ERROR;
		nstatus[4] = MARGIN;
	}
	
	
	// change in a_c or b that causes a reserve vector to change to a margin vector
	double delta_rm;
	if(rflag)
	{
		vector<double> gamma_r;
		vector<int> flags;
		vector<double> g_tmp;
		vector<double> z_tmp(ind[RESERVE].size(), 0);
		select_indices(g, ind[RESERVE], g_tmp);
		select_indices(gamma, ind[RESERVE], gamma_r);
		for(unsigned int i = 0; i < gamma_r.size(); i++)
		{
			gamma_r[i] *= polc;
			flags.push_back(g_tmp[i] >= 0 && gamma_r[i] < 0);
		}
		min_delta(flags, g_tmp, z_tmp, gamma_r, delta_rm, idx);
		if(delta_rm < std::numeric_limits<double>::infinity())
		{
			indss[5] = ind[RESERVE][idx];
			cstatus[5] = RESERVE;
			nstatus[5] = MARGIN;
		}
	}
	else
		delta_rm = std::numeric_limits<double>::infinity();
		
	double tmp[] = {delta_m, delta_e, delta_r, delta_mer, delta_em, delta_rm};
	int min_ind;
	minof(tmp, 6, min_dacb, min_ind);
	indss_ = indss[min_ind];
	cstatus_ = cstatus[min_ind];
	nstatus_ = nstatus[min_ind];
	min_dacb = polc * min_dacb;
}


void IncSVM::initialize(int _m, int _n)
{
	_Rs.initialize(_m, _m);
	_Q.initialize(_n, _m);
	deps = 1e-3;
	max_reserve_vectors = 3000;
	b = 0;
	kernel_evals = 0;
	perturbations = 0;
}



void IncSVM::adddata(int indc, int y_new, double c_new)
{
	a.push_back(0);
	c.push_back(c_new);
	g.push_back(0);
	//_Q.addrow();
	y.push_back(y_new);
	ind[UNLEARNED].push_back(indc);
}


int IncSVM::learn(int indc, int rflag)
{
	_Q.addrow(indc);
	_Q.push_back(indc, y[indc]);
	for(unsigned int i = 0; i < ind[MARGIN].size(); i++)
	{
		_Q.push_back(indc, y[ind[MARGIN][i]]*y[indc]*kernel(indc, ind[MARGIN][i], NULL));
	}
	stopwatch sw;
	sw.reset();
	vector<double> Qc[4];
	double fx = svmeval(indc, Qc);
	g[indc] = y[indc] * fx - 1;
	if(g[indc] >= 0)
	{
		int indco, idx;
		bookkeeping(indc, UNLEARNED, RESERVE, indco, idx);
		return RESERVE;
	}
	double Qcc = kernel(indc, indc, NULL) + deps;
	int converged = 0;
	int cstatus, nstatus;
	double min_delta_param;
	int num_MVs = ind[MARGIN].size();
	stopwatch s1;
	
	while(!converged)
	{
		s1.reset();
		perturbations++;
		vector<double> beta;
		vector<double> gamma(y.size(), 0);
		if(num_MVs > 0)
		{
			//calculate beta
			for(int i = 0; i < _Rs.size; i++)
			{
				double val = _Rs.array[i][0] * y[indc] + vector_star_vector(_Rs.array[i] + 1, Qc[MARGIN]);
				beta.push_back(-val);
			}
			
			//calculate gamma
			int sets[] = {ERROR, RESERVE};
			for(int i = 0; i < 2; i++)
				for(unsigned int j = 0; j < ind[sets[i]].size(); j++)
				{
					double gamma_tmp = vector_star_vector(_Q.get_row(ind[sets[i]][j]), beta);
					gamma[ind[sets[i]][j]] = gamma_tmp + Qc[sets[i]][j];
				}
			double gamma_tmp = vector_star_vector(_Q.get_row(indc), beta);
			gamma[indc] = gamma_tmp + Qcc;
			if(gamma[indc] < 0)
			{
				print_vector(beta);
				_Q.print_row(indc);
				printf("indc = %d gamma_tmp = %lf Qcc = %lf\n", indc, gamma_tmp, Qcc);
				printf("Kernel val is %lf with ind %d\n", kernel(indc, ind[MARGIN][1], NULL), ind[MARGIN][1]);
				fprintf(stderr, "Kernel is not positive\n");
				exit(21);
			}
		}
		else
		{
			beta.push_back(y[indc]);
			for(unsigned int i = 0; i < y.size(); i++)
				gamma[i] = y[indc] * y[i];
		}
		
		int indss;
		watch.reset();
		min_delta_acb(indc, gamma, beta, 1, rflag, min_delta_param, indss, cstatus, nstatus);
		watch.reset();
		if(num_MVs > 0)
		{
			a[indc] += min_delta_param;
			for(unsigned int i = 0; i < ind[MARGIN].size(); i++)
				a[ind[MARGIN][i]] += beta[i + 1] * min_delta_param;
		}
		b += beta[0] * min_delta_param;
		for(unsigned int i = 0; i < g.size(); i++)
			g[i] += min_delta_param * gamma[i];
		converged = (indss == indc);
		if(converged)
		{
			cstatus = UNLEARNED;
			Qc[nstatus].push_back(Qcc);
		}
		else
		{
			int ind_temp = find(indss, ind[cstatus]);
			Qc[nstatus].push_back(Qc[cstatus][ind_temp]);
			Qc[cstatus].erase(Qc[cstatus].begin() + ind_temp);
		}
		int indco, idx;
		bookkeeping(indss, cstatus, nstatus, indco, idx);
		//remove reserve vectors
		for(unsigned int i = 0; i < ind[MARGIN].size(); i++)
			g[ind[MARGIN][i]] = 0;
		watch.reset();
		if(nstatus == MARGIN)
		{
			num_MVs++;
			if(num_MVs > 1)
			{
				double gamma1;
				if(converged)
					gamma1 = gamma[indss];
				else
				{
					matrix_star_vector(-1, _Rs, _Q.get_row(indss), beta);
					gamma1 = kernel(indss, indss, NULL) + deps + vector_star_vector(_Q.get_row(indss), beta);
				}
				gamma.clear();
				gamma.push_back(gamma1);
			}
			updateRQ(beta, gamma, indss);
		}
		else if(cstatus == MARGIN)
		{
			num_MVs--;
			updateRQ(indco);
		}
	}
	//if(indc >=130 && indc <=150){
	//for(int i = 130;i<=150;i++)
	//  printf("After learning %d, %d is in MARGIN:%d, ERROR: %d, RESERVE: %d, g[indc] = %lf\n",indc,i,find(i,ind[MARGIN]),find(i,ind[ERROR]),find(i,ind[RESERVE]),g[i]);
	//}
	return nstatus;
}


int IncSVM::unlearn(int indc)
{
	int num_MVs = ind[MARGIN].size();
	if(g[indc] < 0)
	{
		int idx = find(indc, ind[ERROR]);
		if(idx == -1)
		{
			int a = find(indc, ind[MARGIN]);
			int b = find(indc, ind[RESERVE]);
			int c = find(indc, ind[UNLEARNED]);
			printf("Index %d not found in ERROR, MARGIN:%d, RESERVE:%d, UNLEARNED: %d\n", indc, a, b, c);
		}
		ind[ERROR].erase(ind[ERROR].begin() + idx);
	}
	else if(g[indc] == 0)
	{
		int idx = find(indc, ind[MARGIN]);
		if(idx != -1)
		{
			ind[MARGIN].erase(ind[MARGIN].begin() + idx);
			num_MVs --;
			updateRQ(idx + 1);
		}
		else
		{
			idx = find(indc, ind[RESERVE]);
			ind[RESERVE].erase(ind[RESERVE].begin() + idx);
		}
	}
	else
	{
		int idx = find(indc, ind[RESERVE]);
		ind[RESERVE].erase(ind[RESERVE].begin() + idx);
		ind[UNLEARNED].push_back(indc);
		_Q.clear(indc);
		return RESERVE;
	}
	
	vector<double> Qc[4];
	svmeval(indc, Qc);
	double Qcc = kernel(indc, indc, NULL) + deps;
	int cstatus, nstatus;
	double min_delta_param;
	int converged = 0;
	while(!converged)
	{
		vector<double> beta;
		vector<double> gamma(y.size(), 0);
		if(num_MVs > 0)
		{
			vector<vector<double>*> Q_tmp;
			//calculate beta
			for(int i = 0; i < _Rs.size; i++)
			{
				double val = _Rs.array[i][0] * y[indc] + vector_star_vector(_Rs.array[i] + 1, Qc[MARGIN]);
				beta.push_back(-val);
			}
			
			//calculate gamma
			int sets[] = {ERROR, RESERVE};//,UNLEARNED};
			for(int i = 0; i < 2; i++)
			{
				for(unsigned int j = 0; j < ind[sets[i]].size(); j++)
				{
					double gamma_tmp = vector_star_vector(_Q.get_row(ind[sets[i]][j]), beta);
					gamma[ind[sets[i]][j]] = gamma_tmp + Qc[sets[i]][j];
				}
			}
			double gamma_tmp = vector_star_vector(_Q.get_row(indc), beta);
			gamma[indc] = gamma_tmp + Qcc;
		}
		else
		{
			beta.push_back(y[indc]);
			for(unsigned int i = 0; i < y.size(); i++)
				gamma[i] = y[indc] * y[i];
		}
		int indss;
		min_delta_acb(indc, gamma, beta, -1, 1, min_delta_param, indss, cstatus, nstatus);
		converged = (indss == indc);
		if(num_MVs > 0)
		{
			a[indc] += min_delta_param;
			for(unsigned int i = 0; i < ind[MARGIN].size(); i++)
				a[ind[MARGIN][i]] += beta[i + 1] * min_delta_param;
		}
		b += beta[0] * min_delta_param;
		for(unsigned int i = 0; i < g.size(); i++)
			g[i] += min_delta_param * gamma[i];
		if(!converged)
		{
			int ind_temp = find(indss, ind[cstatus]);
			Qc[nstatus].push_back(Qc[cstatus][ind_temp]);
			Qc[cstatus].erase(Qc[cstatus].begin() + ind_temp);
			int indco, idx;
			bookkeeping(indss, cstatus, nstatus, indco, idx);
			//reserve vectors thingy
			if(nstatus == MARGIN)
			{
				num_MVs++;
				if(num_MVs > 1)
				{
					matrix_star_vector(-1, _Rs, _Q.get_row(indss), beta);
					double gamma1 = kernel(indss, indss, NULL) + deps + vector_star_vector(_Q.get_row(indss), beta);
					gamma.clear();
					gamma.push_back(gamma1);
				}
				updateRQ(beta, gamma, indss);
			}
			else if(cstatus == MARGIN)
			{
				num_MVs --;
				updateRQ(indco);
			}
		}
		else
		{
			ind[nstatus].push_back(indc);
			if(nstatus == MARGIN)
			{
				num_MVs++;
				if(num_MVs > 1)
				{
					matrix_star_vector(-1, _Rs, _Q.get_row(indss), beta);
					double gamma1 = kernel(indss, indss, NULL) + deps + vector_star_vector(_Q.get_row(indss), beta);
					gamma.clear();
					gamma.push_back(gamma1);
				}
				updateRQ(beta, gamma, indss);
			}
		}
		for(unsigned int i = 0; i < ind[MARGIN].size(); i++)
			g[ind[MARGIN][i]] = 0;
	}
	_Q.clear(indc);
	return MARGIN;
	//if(indc >=130 && indc <=150){
	//for(int i = 130;i<=150;i++)
	//  printf("After unlearning %d, %d is in MARGIN:%d, ERROR: %d, RESERVE: %d, g[indc] = %lf\n",indc,i,find(i,ind[MARGIN]),find(i,ind[ERROR]),find(i,ind[RESERVE]),g[i]);
	//}
	
}


double IncSVM::get_w2()
{
	unsigned int i;
	double s = 0;
	for(i = 0; i < a.size(); i++)
	{
		if(a[i] <= 0)
			continue;
		if(y[i] > 0)
			s += a[i] * (g[i] + 1.0);
		else
			s += y[i] * a[i] * (g[i] - 1.0);
	}
	return s / 2.0;
}


vector<double>& IncSVM::get_model(double& b_)
{
	b_ = b;
	return a;
}


void IncSVM::setY(int idx, int y_new)
{
	int y_old = (int)y[idx];
	y[idx] = y_new;
	for(int i = 0; i < _Q.get_row_size(idx); i++)
		_Q.get_row(idx)[i] *= y_old * y_new;
}


