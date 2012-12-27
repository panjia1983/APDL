#ifndef _PB_MINKOWSKI_GRID_
#define _PB_MINKOWSKI_GRID_

#ifdef WIN32
#pragma warning(disable:4018)
#endif

#include "cell.h"
#include "points.h"

class ms_grid
{
public:

	ms_grid(){ cells=NULL; }

	void build(float box[6], float d){

		if(cells!=NULL) return; //has been build

		cell_total=1;
		cell_size=d;

		for(int i=0;i<3;i++){
			int id1=i*2;
			int id2=id1+1;
			bbox[id1]=box[id1];
			bbox[id2]=box[id2];
			bbox_size[i]=box[id2]-box[id1];
			cell_count[i]=(unsigned int) ceil( ((float)bbox_size[i])/d );
			cell_total*=cell_count[i];
		}
		
		cell_yz=cell_count[1]*cell_count[2];

		//allocate
		cells=new ms_cell[cell_total];
		assert(cells); //make sure we have enough memory
	}

	~ms_grid(){ destroy(); }

	void destroy(){ delete [] cells; cells=NULL; }

	//register points
	void registering(const points& P,const points& Q);

	// first wave
	// compute bd_cells from external cells
	void first_wave();

	// second wave
	// compute more bd_cells by identifying gaps
	void second_wave();


	//get points on the boundary
	void getBoundaryPts(list<Point3d>& M);

private:

	//classify the most external grids to external or boundary
	void classify_grid_boundary(list<int>& ext);

	//given a point, return its id
	unsigned int pt2id(const Point3d& p)
	{
		int idx=(int)((p[0]-bbox[0])/cell_size);
		int idy=(int)((p[1]-bbox[2])/cell_size);
		int idz=(int)((p[2]-bbox[4])/cell_size);
		return idx*cell_yz+idy*cell_count[2]+idz;
	}

	void getNeighbors(int cid, list<int>& nei)
	{
		nei.clear();
		unsigned int cx=cid/cell_yz;
		unsigned int r=cid-cx*cell_yz;
		unsigned int cy=r/cell_count[2];
		unsigned int cz=r-cy*cell_count[2];

		//get the neighboring cells
		unsigned int nid=0;
		for(int x=-1;x<2;x++){
			int _x=cx+x;
			if(_x<0 || _x>= cell_count[0]) continue;
			for(int y=-1;y<2;y++){
				int _y=cy+y;
				if(_y<0 || _y>= cell_count[1]) continue;
				for(int z=-1;z<2;z++){
					int _z=cz+z;
					if(_z<0||_z>=cell_count[2]||(x==0&&y==0&&z==0)) continue;
					nid=_x*cell_yz+_y*cell_count[2]+_z;
					nei.push_back(nid);
				}//end z
			}//end y
		}//end x
	}

	//bbox
	float bbox[6];
	float bbox_size[3]; //computed from bbox

	//cells
	ms_cell *cells;
	unsigned int cell_count[3];
	unsigned int cell_total;
	unsigned int cell_yz;    //size of cell in yz plane
	float cell_size;

	list<int> bd_cells;

	//associated points
	const points* P;
	const points* Q;
};

#endif //_PB_MINKOWSKI_GRID_