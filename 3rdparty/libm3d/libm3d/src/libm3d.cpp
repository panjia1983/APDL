#include "libm3d.h"
#include "octree.h"

#include <pthread.h>
pthread_mutex_t mutexBD = PTHREAD_MUTEX_INITIALIZER;

double current_rot[3][3]={{1,0,0},{0,1,0},{0,0,1}};

void collectBoundaryPts_tryall(const model& P, const model& Q,
							   mksum* M, const LibM3dParam& param,
							   MKPTS& allBds)
{
	double pos[3];
	MKPTS& pts=M->pts;

	int count=0;
	for(MKPTIT i= pts.begin(); i!= pts.end(); i++, count++)
	{
		M->getMPos_offset(*i, pos);
		if( !is_in_collision(P.cd,Q.cd,pos,current_rot) )
		{
			if(param.threadsize > 1) pthread_mutex_lock( &mutexBD );
			allBds.push_back(*i);
			if(param.threadsize > 1) pthread_mutex_unlock( &mutexBD );
		}
		if(count%10000==0)
			cout<< "." << flush;
	}//for i

	if(param.threadsize > 1) pthread_mutex_lock( &mutexBD );
	if(param.threadsize > 1) pthread_mutex_unlock( &mutexBD );
}

//-----------------------------------------------------------------------------
// validate points using an octree

bool isBoundaryCell(const model& P, const model& Q,
						   mksum* M, const LibM3dParam& param,
						   octree_cell* cell, list<mksum_pt>& tmp)
{
	MKPTS& pts=cell->points;
	//check if empty cell
	if(pts.empty()) return false;

	bool collision_free=false;
	double pos[3];

	//record total cd calls
	if(param.threadsize>1) pthread_mutex_lock( &mutexBD );
	if(param.threadsize>1) pthread_mutex_unlock( &mutexBD );

	for(MKPTIT i=pts.begin();i!=pts.end();i++)
	{
		M->getMPos_offset(*i,pos);

		if( !is_in_collision(P.cd,Q.cd,pos,current_rot) )
		{
			//lock here
			if(param.threadsize>1) pthread_mutex_lock( &mutexBD );
			tmp.push_back(*i);
			if(param.threadsize>1) pthread_mutex_unlock( &mutexBD );
			//lock here
			collision_free=true;
		}
	}//end i

	return collision_free;
}

void collectBoundaryPts_using_octree(const model& P, const model& Q,
											mksum* M, const LibM3dParam& param, 
											octree& O,
											MKPTS& allBds)
{
	list<octree_cell*> open;
	O.root.getBoundaryCells(open);

	for(list<octree_cell*>::iterator i=open.begin();i!=open.end();i++)
		(*i)->visited=true;

	int counter=0;

	while(!open.empty())
	{
		octree_cell * c=open.front();
		open.pop_front();
		bool bd=isBoundaryCell(P, Q, M, param, c, allBds);
		if(bd){
			list<octree_cell*> tmp;
			c->getNeighbors(tmp);
			for(list<octree_cell*>::iterator i=tmp.begin();i!=tmp.end();i++)
				open.push_back(*i);
		}

		counter++;
		if(counter%100==0)
			cout<<"."<<flush;

	}//end while
	cout<<endl;
}



struct MKSData
{
	model* P;
	model* Q;
	mksum M;
	LibM3dParam param;
	MKPTS* allBds;

	MKSData() : M(mksum())
	{
		P = NULL;
		Q = NULL;
		allBds = NULL;
	}
};

void* MKS(void* _data)
{
	MKSData* data = (MKSData*)_data;
	mksum * M = &(data->M);
	if(!data->param.bOctfilter)
	{
		collectBoundaryPts_tryall(*(data->P), *(data->Q), M, data->param, *(data->allBds));
	}
	else    //build octree and collect more boundary points
	{
		//build octree
		octree O(M);
		O.build();
		//collect points
		collectBoundaryPts_using_octree(*(data->P), *(data->Q), M, data->param, O, *(data->allBds));
		O.destroy();
	}
	//---------------------------------------------------------------

	M->destroy();
	return _data;
}


void ComputePointMKSum(model& P, model& Q, 
					   double R[3][3],
					   const LibM3dParam& param, MKPTS& allBds)
{
	mksum M;

	if(param.bNormalfilter)
		M.build_with_filter(&P, &Q);
	else
		M.build_without_filter(&P, &Q);

	MKSData mksdata;

	mksdata.M = M;
	mksdata.P = &P;
	mksdata.Q = &Q;
	mksdata.allBds = &allBds;
	mksdata.param = param;

	if(!param.bCDfilter)
		allBds.swap(M.pts);
	else
	{
		if(param.threadsize <= 1)
		{
			MKS((void*)&mksdata);
		}
		else
		{
			//create threads
			//split M
			vector<MKSData> mksdatav(param.threadsize);
			for(int i=0; i<param.threadsize;i++)
			{
				mksdatav[i].M.P = M.P;
				mksdatav[i].M.Q = M.Q;
				mksdatav[i].M.pt_size = 0;
				mksdatav[i].P = &P;
				mksdatav[i].Q = &Q;
				mksdatav[i].allBds = &allBds;
				mksdatav[i].param = param;
			}

			int threadindex=0;
			for(MKPTIT ptr=M.pts.begin();ptr!=M.pts.end();ptr++){
				mksdatav[threadindex].M.pts.push_back(*ptr);
				mksdatav[threadindex].M.pt_size++;
				threadindex=(threadindex+1)%param.threadsize;
			}//end for

			M.destroy();

			//create threads
			cout<<"\t- Create "<<param.threadsize<<" threads"<<endl;
			vector<pthread_t> threads(param.threadsize,pthread_t());
			for(int i=0;i<param.threadsize;i++)
				pthread_create(&threads[i], NULL, MKS, (void*)&mksdatav[i]);

			//wait for threads to finish
			for(int i=0;i<param.threadsize;i++)
				pthread_join(threads[i], NULL);			
		}
	}

	P.destroy();
	Q.destroy();
}
