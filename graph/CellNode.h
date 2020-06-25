/*
 * TreeNode.h
 *
 *  Created on: 2016年8月28日
 *      Author: lostrong
 */

#ifndef SRC_TREENODE_H_
#define SRC_TREENODE_H_
#include<vector>
using namespace std;
class CellNode {
	// set access as public for convenience
public:
	int coor_x_min;
	int coor_y_min;
	int coor_x_max;
	int coor_y_max;
	vector<int> inside_nodes;

	CellNode(int coor_x_min_p, int coor_x_max_p, int coor_y_min_p,
			int coor_y_max_p, vector<int> inside_nodes_p) {
		coor_x_min = coor_x_min_p;
		coor_x_max = coor_x_max_p;
		coor_y_min = coor_y_min_p;
		coor_y_max = coor_y_max_p;
		for (int i = 0; i < inside_nodes_p.size(); i++)
			inside_nodes.push_back(inside_nodes_p[i]);
	}
};


class Shell {
	// set access as public for convenience
public:
	int cell_id;
	int coor_x_min;
	int coor_y_min;
	int coor_x_max;
	int coor_y_max;

	Shell(int cell_id_p, int coor_x_min_p, int coor_x_max_p, int coor_y_min_p,
			int coor_y_max_p) {
		cell_id = cell_id_p;
		coor_x_min = coor_x_min_p;
		coor_x_max = coor_x_max_p;
		coor_y_min = coor_y_min_p;
		coor_y_max = coor_y_max_p;

	}
};

struct KNode {
	int id;
	int dis;
	KNode() :
			id(0), dis(0) {
	}
	KNode(int a, int b) :
			id(a), dis(b) {
	}
	static bool compare_dis(const KNode&a, const KNode&b) {
		return a.dis < b.dis;
	}
	KNode(const KNode& C) {
		id = C.id;
		dis = C.dis;
	}
	KNode& operator=(const KNode& s) {
		id = s.id;
		dis = s.dis;
		return *this;
	}
//	bool operator<( KNode a, KNode b){
//	    if(a.dis==b.dis) return a.id<b.id;
//	    return a.dis<b.dis;
//	}

};

#endif /* SRC_TREENODE_H_ */
