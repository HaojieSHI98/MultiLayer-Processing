//
// Created by lostrong on 18/3/10.
//

#ifndef TOAIN_DIJKSTRAKNN_H
#define TOAIN_DIJKSTRAKNN_H
#include "../graph/CellNode.h"
#include "../dijk/DijkstraQueue.h"
#include "../para/GlobalVariables.h"
#include "../util/GeneralUtility.h"


// in practice, vector dis and visited can be reused for different queries,
//so set them as parameters for simulating practical scenarios
vector<KNode> DijkstraKNNQuery(int k, int query, long* dis,
                       int* visited, DijkstraQueue* q, int*& car_nodes) {
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
//    cout<<"visit "<<visitedlist.size()<<" nodes"<<endl;
    return result;
}

vector<KNode> DijkstraKNNQueryThreshold(int k, int query, long* dis,
                               int* visited, DijkstraQueue* q, int*& car_nodes, vector<int>& threshold, int query_id) {
    vector<int> visitedlist;
    vector<KNode> result;
    dis[query] = 0;
    visited[query] = 1;
    visitedlist.push_back(query);
    q->push(query, 0);

    while (!q->empty()) {
        int u = q->pop();
        if(dis[u]> threshold[query_id]) break;
        if (car_nodes[u]) {
            KNode* c = new KNode(u, dis[u]);
            result.push_back(*c);
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
//    cout<<"visit "<<visitedlist.size()<<" nodes"<<endl;
    return result;
}

void DijkstraKNNDelete(int object_node, int*& car_nodes) {
//    cout<<"enter delete function"<<endl;
    if(car_nodes[object_node]==0){
        cout<<"Not proper delete of object "<<object_node<<"! The object does not exist"<<endl;
        stop_here();
    }
    car_nodes[object_node]=0;
}

void DijkstraKNNInsert(int object_node, int*& car_nodes) {
//    cout<<"enter insert function"<<endl;
    if(car_nodes[object_node]==1){
        cout<<"Not proper insert of object "<<object_node<<"! The object exists"<<endl;
        stop_here();
    }
//    cout<<"before"<<endl;
    car_nodes[object_node]=1;
//    cout<<"after"<<endl;
}

#endif //TOAIN_DIJKSTRAKNN_H
