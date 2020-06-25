/*
 * DijkstraQueue.h
 *
 *  Created on: 2017年2月14日
 *      Author: lostrong
 */

#ifndef SRC_DIJKSTRAQUEUE_H_
#define SRC_DIJKSTRAQUEUE_H_
#include <limits.h>


class DijkstraQueue
{
	int *label;
	int *dis;
	int *state;
	int num;
	int n;
	int l(int x){ return x << 1; }
	int r(int x){ return (x << 1) | 1; }
	int f(int x){ return x >> 1; }

	void swap_it(int x, int y)
	{
		state[label[x]] = y;
		state[label[y]] = x;
		std::swap(label[x], label[y]);
		std::swap(dis[x], dis[y]);
	}

	void keep_up(int x)
	{
		int fx = f(x);
		while (fx > 0 && dis[fx] > dis[x]){
			swap_it(fx, x);
			x = fx;
			fx = f(x);
		}
	}

	void keep_down(int x)
	{
		int lx = l(x), rx = r(x), a = 0;
		while (rx <= num){
			a = dis[lx] < dis[rx] ? lx : rx;
			if (dis[x] <= dis[a]) break;
			swap_it(x, a);
			x = a;
			lx = l(x), rx = r(x);
		}
		if (lx == num && dis[lx] < dis[x])
			swap_it(x, lx);
	}

public:
	DijkstraQueue(int n,int s);
	void init(int s);
	DijkstraQueue(int n);
	void deleteDijkstraQueue();
	~DijkstraQueue();
	void clear(){ num = 0; }
	bool empty(){ return num == 0; }
	void push(int s, int d){
//		cout << "num " << num << endl;
		state[s] = ++num;
		dis[num] = d;
		label[num] = s;
		keep_up(num);
	}
	int pop(){
		int u = label[1];
		swap_it(1, num);
		--num;
		if (num > 0) keep_down(1);
		return u;
	}
	int top(){
		return label[1];
	}
	int size(){ return num; }
	void decrease_priority(int v, int d){
		int u = state[v];
		dis[u] = d;
		keep_up(u);
	}
};

DijkstraQueue::DijkstraQueue(int n, int s) {
	dis = new int[n + 1];
	state = new int[n + 1];
	label = new int[n + 1];
	this->n = n;
	dis[0] = 0;
	dis[1] = 0;
	state[s] = 1;
	label[1] = s;
	for (int i = 2; i < n + 1; ++i)
		dis[i] = INT_MAX;
	for (int i = 1; i < s; ++i) {
		label[i + 1] = i;
		state[i] = i + 1;
	}
	for (int i = s + 1; i < n + 1; ++i) {
		label[i] = i;
		state[i] = i;
	}
	num = n;
}




void DijkstraQueue::init(int s) {
	num = 0;
	dis[0] = 0;
	dis[1] = 0;
	state[s] = 1;
	label[1] = s;
	for (int i = 2; i < n + 1; ++i)
		dis[i] = INT_MAX;
	for (int i = 1; i < s; ++i) {
		label[i + 1] = i;
		state[i] = i + 1;
	}
	for (int i = s + 1; i < n + 1; ++i) {
		label[i] = i;
		state[i] = i;
	}
	num = n;
}



DijkstraQueue::DijkstraQueue(int n) {
	dis = new int[n + 1];
	state = new int[n + 1];
	label = new int[n + 1];
	this->n = n;
	num = 0;
	dis[0] = 0;
	state[0] = 0;
	label[0] = 0;
}



void DijkstraQueue::deleteDijkstraQueue(){
	delete[] dis;
	delete[] state;
	delete[] label;
}
DijkstraQueue::~DijkstraQueue() {
	delete[] dis;
	delete[] state;
	delete[] label;
}



#endif /* SRC_DIJKSTRAQUEUE_H_ */
