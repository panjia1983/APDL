#ifndef _PB_MINKOWSKI_CELL_
#define _PB_MINKOWSKI_CELL_

#include <list>
using namespace std;

//a registration of a minkowski point to a cell
struct regi
{
	regi(unsigned int p, unsigned int q){ pid=p; qid=q; }
	unsigned int pid, qid; 
};

struct ms_cell
{
	ms_cell(){ visited=false; }

	void add(unsigned int pid,unsigned int qid)
	{ enclosed.push_back(regi(pid,qid)); }

	bool isempty() const { return enclosed.empty(); }

	// check if the cell is an incident cell
	bool is_incident() const;

	// check if the cell is an incident cell
	bool has_gap() const;

	list<regi> enclosed;  // a list of enclosed minkowski points

	bool visited; //1 byte....
};

#endif //_PB_MINKOWSKI_CELL_





