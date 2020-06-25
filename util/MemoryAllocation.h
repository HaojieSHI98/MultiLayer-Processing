/*
 * Memory_Allocation.h
 *
 *  Created on: 2017年1月16日
 *      Author: lostrong
 */

#ifndef SRC_MEMORYALLOCATION_H_
#define SRC_MEMORYALLOCATION_H_

struct mem_struct{
	long* dist;
	long* dist1;
	int* visited;
	int* visited1;
	DijkstraQueue* q;
	DijkstraQueue* q1;
	int* marked_first;
	int* colors;
	int* pre;
};

//This is kind of a package allocation for Dijkstra's algorithm
void allocate_mem(mem_struct& mems, int size){
	mems.dist=new long[size];
	mems.dist1 = new long[size];
	mems.visited=new int[size];
	mems.visited1=new int[size];
	mems.q=new DijkstraQueue(size);
	mems.q1=new DijkstraQueue(size);
	mems.marked_first=new int[size];
	mems.colors=new int[size];
	mems.pre=new int[size];
	memset(mems.visited, 0, sizeof(int)*size);
	memset(mems.visited1, 0, sizeof(int)*size);
	memset(mems.marked_first, 0, sizeof(int)*size);
	memset(mems.colors, 0, sizeof(int)*size);
}

void delete_mems(mem_struct& mems){
	delete[] mems.dist;
	delete[] mems.dist1;
	delete[] mems.visited;
	delete[] mems.visited1;
    mems.q->deleteDijkstraQueue();
    mems.q1->deleteDijkstraQueue();
	delete[] mems.marked_first;
	delete[] mems.colors;
	delete[] mems.pre;
}


#endif /* SRC_MEMORYALLOCATION_H_ */
