/*
 * NaiveKNN.h
 *
 *  Created on: 2016年8月29日
 *      Author: lostrong
 */

#ifndef SRC_NAIVEKNN_H_
#define SRC_NAIVEKNN_H_

#include<vector>
#include<unordered_map>
#include<set>
#include<limits.h>
#include "../object/CarPos.h"
#include "../util/GeneralUtility.h"
using namespace std;

// in practice, vector dis and visited can be reused for different queries,
//so set them as parameters for simulating practical scenarios
vector<KNode> naiveKNN(int k, vector<int>& cars, int query, long* dis,
		int* visited, DijkstraQueue* q, int* car_nodes) {

    // the parameter cars is not used
//	if (cars.size() < k) {
//		cout << "The number of cars is less than k!" << endl;
//		stop_here();
//	}
	vector<int> visitedlist;
	vector<KNode> result;
	dis[query] = 0;
	visited[query] = 1;
	visitedlist.push_back(query);
	q->push(query, 0);

	while (!q->empty()) {
		int u = q->pop();
		if (car_nodes[u]) {
			KNode* c = new KNode(u, dis[u]);
			result.push_back(*c);
			if (result.size() >= k) {
				break;
			}
		}
		for (KNode knode : level_scs_arr[0][u]) {
			int v = knode.id;
			if (visited[v]) {
				if (knode.dis < dis[v] - dis[u]) {
					dis[v] = dis[u] + knode.dis;
					q->decrease_priority(v, dis[v]);
				}
			} else {
				dis[v] = dis[u] + knode.dis;
				q->push(v, dis[v]);
				visited[v] = 1;
				visitedlist.push_back(v);
			}

		}
	}
	for (int item : visitedlist)
		visited[item] = 0;
	q->clear();
//    cout<<"visit "<<visitedlist.size()<<" nodes"<<endl;
	return result;
}

vector<KNode> naiveKNN_carpos(int k, vector<CarPos>& carposes,unordered_map<int, vector<PointCar> >& car_hash, int query,
		long* dis, int* visited, DijkstraQueue* q) {
	q->clear();
	vector<int> visitedlist;
	vector<KNode> result;



	dis[query] = 0;
	visited[query] = 1;
	visitedlist.push_back(query);
	q->push(query, 0);

	long start = clock();
	while (!q->empty()) {
		int u = q->pop();
		if (u > test_n) {
			KNode knode(u, dis[u]);
			result.push_back(knode);
//			if (result.size() >= k) {
//				break;
//			}
			continue;
		}

		if (result.size() >= k && dis[u] > result[k - 1].dis)
			break;
		if (car_hash.find(u) != car_hash.end()) {
			for (PointCar pc : car_hash[u]) {
				int carid = pc.carid + test_n;
				if (visited[carid]) {
					if (dis[u] + pc.dis < dis[carid]) {
						dis[carid] = dis[u] + pc.dis;
						q->decrease_priority(carid, dis[carid]);
					}
				} else {
					dis[carid] = dis[u] + pc.dis;
					q->push(carid, dis[carid]);
					visited[carid] = 1;
					visitedlist.push_back(carid);
				}
			}
		}

		for (int i = vStart[u]; i != 0; i = eNext[i]) {
			int v = eNode[i];

			if (visited[v]) {
				if (eDis[i] < dis[v] - dis[u]) {
					dis[v] = dis[u] + eDis[i];
					q->decrease_priority(v, dis[v]);
				}
			} else {
				dis[v] = dis[u] + eDis[i];
				q->push(v, dis[v]);
				visited[v] = 1;
				visitedlist.push_back(v);
			}

		}

	}
	for (int i = 0; i < visitedlist.size(); i++)
		visited[visitedlist[i]] = 0;
	q->clear();
	vector<KNode> tmp;
	for (int i = 0; i < k && i < result.size(); i++)
		tmp.push_back(result[i]);

	return tmp;
}



vector<KNode> naiveKNN_withHash(int k, vector<int>& cars, int query, int* dis,
		int* visited, DijkstraQueue* q, unordered_map<int, int>& car_hash) {
	if (cars.size() < k) {
		cout << "The number of cars is less than k!" << endl;
        stop_here();
	}
	vector<int> visitedlist;
	vector<KNode> result;
	dis[query] = 0;
	visited[query] = 1;
	visitedlist.push_back(query);
	q->push(query, 0);
	long start = clock();
	while (!q->empty()) {
		int u = q->pop();
		if (car_hash.find(u) != car_hash.end()) {
			KNode knode(u, dis[u]);
			result.push_back(knode);
			if (result.size() >= k) {
				break;
			}
		}
		for (int i = vStart[u]; i != 0; i = eNext[i]) {
			int v = eNode[i];
			if (visited[v]) {
				if (eDis[i] < dis[v] - dis[u]) {
					dis[v] = dis[u] + eDis[i];
					q->decrease_priority(v, dis[v]);
				}
			} else {
				dis[v] = dis[u] + eDis[i];
				q->push(v, dis[v]);
				visited[v] = 1;
				visitedlist.push_back(v);
			}

		}
	}


	for (int item : visitedlist)
		visited[item] = 0;
	q->clear();

	return result;
}

#endif /* SRC_NAIVEKNN_H_ */
