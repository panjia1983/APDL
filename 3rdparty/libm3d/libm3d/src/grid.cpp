#include "grid.h"



//register points
void ms_grid::registering(const points& P,const points& Q)
{
	unsigned int id;
	this->P=&P;
	this->Q=&Q;
	const point& ref=P.getRef();

	for(unsigned int i=0;i<Q.size;i++){
		const point& q=Q.pts[i];
		Vector3d offset=q.p-ref.p;

		for(unsigned int j=0;j<P.size;j++){
			const point& p=P.pts[j];

			Point3d np=p.p+offset;
			Vector3d v=q.p-np;
			if(v*q.n<0) continue; //opposite....

			id=pt2id(np);
			cells[id].add(j,i);

		}//end j
	}//end i
}

//-----------------------------------------------------------------------------

inline void classify
(unsigned int id, ms_cell& c, list<int>& ext, list<int>& bd)
{
	if(c.visited) return;
	if( c.isempty() )
		ext.push_back(id); //external cell
	else
		bd.push_back(id); //boundary
	c.visited=true;
}

// first wave
// compute bd_cells from external cells
void ms_grid::first_wave()
{
	list<int> ext;
	list<int> nei; //neighbors

	classify_grid_boundary(ext);
	while(!ext.empty()){
		int cid=ext.front();
		ext.pop_front();
		getNeighbors(cid,nei);
		//check the neighboring cells
		for(list<int>::iterator i=nei.begin();i!=nei.end();i++)
			classify(*i,cells[*i],ext,bd_cells);
	}//end while
}

//-----------------------------------------------------------------------------

// second wave
// compute more bd_cells by identifying gaps
void ms_grid::second_wave()
{

}

//-----------------------------------------------------------------------------

void ms_grid::getBoundaryPts(list<Point3d>& M)
{
	const point& ref=P->getRef();

	for(list<int>::iterator i=bd_cells.begin();i!=bd_cells.end();i++){
		ms_cell& c=cells[*i];
		list<regi>& en=c.enclosed;
		for(list<regi>::iterator i=en.begin();i!=en.end();i++){
			Point3d p=P->pts[i->pid].p+(Q->pts[i->qid].p-ref.p);
			M.push_back(p);
		}
	}//end i
}

//-----------------------------------------------------------------------------

//classify the most external grids to external or boundary
void ms_grid::classify_grid_boundary(list<int>& ext)
{
	unsigned int x_id, y_id, z_id, id1, id2;

	//the the first wave front
	y_id=0;
	x_id=(cell_count[0]-1)*cell_yz;
	for(unsigned int y=0; y<cell_count[1]; y++, y_id+=cell_count[2]){
		for(unsigned int z=0;z<cell_count[2];z++){
			id1=y_id+z;
			id2=x_id+id1;
			classify(id1, cells[id1], ext, bd_cells);
			classify(id2, cells[id2], ext, bd_cells);
		}//end z
	}//end y

	x_id=0;
	y_id=(cell_count[1]-1)*cell_count[2];
	for(unsigned int x=0; x<cell_count[0]; x++, x_id+=cell_yz){
		for(unsigned int z=0;z<cell_count[2];z++){
			id1=x_id+z;
			id2=y_id+id1;
			classify(id1, cells[id1], ext, bd_cells);
			classify(id2, cells[id2], ext, bd_cells);
		}//end z
	}//end y

	x_id=0;
	y_id=0;
	z_id=cell_count[2]-1;
	for(unsigned int x=0; x<cell_count[0]; x++, x_id+=cell_yz){
		for(unsigned int y=0;y<cell_count[1];y++, y_id+=cell_count[2]){
			id1=x_id+y_id;
			id2=z_id+id1;
			classify(id1, cells[id1], ext, bd_cells);
			classify(id2, cells[id2], ext, bd_cells);
		}//end z
	}
}




