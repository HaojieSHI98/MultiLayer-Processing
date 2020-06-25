/*
 * CarPos.h
 *
 *  Created on: 2016年9月18日
 *      Author: lostrong
 */

#ifndef SRC_CARPOS_H_
#define SRC_CARPOS_H_

class CarPos {
public:
	int id;
	int node1;
	int node2;
	int dis1;
	int dis2;
	time_t time;
	CarPos() {

	}
	CarPos(int id, int node1, int node2, int dis1, int dis2, time_t time) {
		this->id = id;
		this->node1 = node1;
		this->node2 = node2;
		this->dis1 = dis1;
		this->dis2 = dis2;
		this->time = time;
	}

	static bool compare_time(const CarPos&a, const CarPos&b) {
		if(a.time<b.time) return true;
		else if(a.time==b.time) return a.id<b.id;
		else return false;
//		return a.time < b.time;
	}
};

class PointCar {
public:
	int carid;
	int dis;
	PointCar(int carid, int dis) {
		this->carid = carid;
		this->dis = dis;
	}
	static bool compare_dis(const PointCar&a, const PointCar&b) {
			return a.dis < b.dis;
		}
};


class NodeRank {
public:
	int node_order;
	int in_node_order;
	int id;
	NodeRank(int node_order, int in_node_order, int id) {
		this->node_order = node_order;
		this->in_node_order = in_node_order;
		this->id=id;
	}
	static bool compare(const NodeRank&a, const NodeRank&b) {
			if(a.node_order<b.node_order) return true;
			else if(a.node_order==b.node_order && a.in_node_order<b.in_node_order) return true;
			else if(a.node_order==b.node_order && a.in_node_order==b.in_node_order && a.id<b.id) return true;
			else return false;

		}
	static bool compare_level(const NodeRank&a, const NodeRank&b) {
				if(a.node_order<b.node_order) return true;
				else return false;

			}
};

#endif /* SRC_CARPOS_H_ */
