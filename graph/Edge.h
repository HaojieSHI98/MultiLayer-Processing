/*
 * Edge.h
 *
 *  Created on: 2017年2月14日
 *      Author: lostrong
 */

#ifndef SRC_EDGE_H_
#define SRC_EDGE_H_



struct Edge
{
	int to;
	int dis;
	Edge() :to(0), dis(0){}
	Edge(int to, int dis) :to(to), dis(dis){}
};



#endif /* SRC_EDGE_H_ */
