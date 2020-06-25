#ifndef SRC_OBJECTFUNCTIONS_H_
#define SRC_OBJECTFUNCTIONS_H_
#include<vector>
#include "../para/GlobalVariables.h"
#include "CarPos.h"
#include "../util/GeneralUtility.h"
#include "../RoadKNNUtilities.h"
using namespace std;

string add_dataspace_to_name(const char* name) {
	string res_name = data_space;
	res_name = res_name.append(name);
	return res_name;
}
void get_car_position_ny() {
	vector<CarPos> carposes;
	char buffer1[100], buffer2[100], buffer3[100], buffer4[100];
	unordered_map<long, long> id_hash;
	FILE* fInput = fopen("green_tripdata_2015-01.csv", "r");
	if (fInput == NULL) {
		cout << "fInput==NULL" << endl;
		char c;
		cin >> c;
	}
	double x1 = 1, y1 = 1, x2 = 1, y2 = 1;
	int pcnt = 1;
	int id = 1;
	double nw_dist = 1;
	long x3, y3;
	fscanf(fInput, "%*s");
	int id_cnt = 1;
	int count = 0;
	while (fscanf(fInput, "%s %s %13s%lf,%lf,%lf,%lf,%d,%lf,%s", buffer1,
			buffer2, buffer3, &x1, &y1, &x2, &y2, &pcnt, &nw_dist, buffer4)
			== 10) {
		count++;
		//	while (fscanf(fInput, "%s %s %8s,%d,%lf,%lf,%lf,%d,%c,%lf,%lf,%s", buffer1, buffer2, buffer3, &pcnt, &nw_dist, &x1, &y1, &d2, &c1, &x2, &y2, buffer4) == 12){

		x3 = (int) (x1 * 1000000);
		y3 = (int) (y1 * 1000000);
		int min_dist = INT_MAX;
		int min_node1 = 0;
		int min_node2 = 0;
		int cell_x = (x3 - x_min) / x_len;
		int cell_y = (y3 - y_min) / y_len;
		if (cell_x < 0 || cell_y < 0 || cell_x >= cell_num
				|| cell_y >= cell_num)
			continue;
		for (int t = 0; t < cells[cell_x][cell_y].size(); t++) {
			int u = cells[cell_x][cell_y][t];
			int x1 = x[u - 1];
			int y1 = y[u - 1];
			double d1 = sqrt(
					((x1 - x3) / 1000000.0) * ((x1 - x3) / 1000000.0)
							+ ((y1 - y3) / 1000000.0)
									* ((y1 - y3) / 1000000.0));

			for (int i = vStart[u]; i != 0; i = eNext[i]) {
				int v = eNode[i];
				int x2 = x[v - 1];
				int y2 = y[v - 1];
				double d2 = sqrt(
						((x1 - x2) / 1000000.0) * ((x1 - x2) / 1000000.0)
								+ ((y1 - y2) / 1000000.0)
										* ((y1 - y2) / 1000000.0));
				double d3 = sqrt(
						((x2 - x3) / 1000000.0) * ((x2 - x3) / 1000000.0)
								+ ((y2 - y3) / 1000000.0)
										* ((y2 - y3) / 1000000.0));
				double p = (d1 + d2 + d3) / 2;
				double p1 = p - d1;
				double p2 = p - d2;
				double p3 = p - d3;
				if (d3 < 0.00001)
					continue;
				double h = sqrt(p * p1 * p2 * p3) * 2.0 / d3;
				if (h < min_dist) {
					min_dist = h;
					min_node1 = u;
					min_node2 = v;
				}
			}
		}
		if (min_dist < INT_MAX) {
			double d1 = ((x3 - x[min_node1 - 1]) / 1000000.0)
					* ((x3 - x[min_node1 - 1]) / 1000000.0)
					+ ((y3 - y[min_node1 - 1]) / 1000000.0)
							* ((y3 - y[min_node1 - 1]) / 1000000.0);
			double d2 = ((x3 - x[min_node2 - 1]) / 1000000.0)
					* ((x3 - x[min_node2 - 1]) / 1000000.0)
					+ ((y3 - y[min_node2 - 1]) / 1000000.0)
							* ((y3 - y[min_node2 - 1]) / 1000000.0);
			double dis1 = sqrt(d1 - min_dist * min_dist);
			double dis2 = sqrt(d2 - min_dist * min_dist);
			double edge_len = get_distance_lat_lng(double(x[min_node1]),
					double(y[min_node1]), double(x[min_node2]),
					double(y[min_node2]), 1000000.0);
			if (id_hash.find(id) == id_hash.end()) {
				id_hash[id] = id_cnt;
				id_cnt++;
			}
			carposes.push_back(
					CarPos(id_hash[id], min_node1, min_node2,
							(int) (edge_len * dis1 / (dis1 + dis2)),
							(int) (edge_len * dis2 / (dis1 + dis2)), 100));
			id++;
		}

	}

	sort(carposes.begin(), carposes.end(), CarPos::compare_time);
	string file_name1 = data_space;
	file_name1 = file_name1.append("car_positions.txt");
	FILE* file1 = fopen(file_name1.c_str(), "w");
	for (CarPos car : carposes) {
		fprintf(file1, "%d %d %d %d %d %ld\n", car.id, car.node1, car.node2,
				car.dis1, car.dis2, car.time);
	}
	fclose(fInput);
	fclose(file1);
}

void get_car_position_random() {

	string file_name1 = data_space;
	file_name1 = file_name1.append("car_positions.txt");
	FILE* file1 = fopen(file_name1.c_str(), "w");
	int rand_num = 1200;
	int rand_num_cnt = 1;
	for (int i = 0; i < rand_num;) {
		int node1 = get_rand_node();
		if (valid_node_hash.find(node1) != valid_node_hash.end()) {
			for (int i = vStart[node1]; i != 0; i = eNext[i]) {
				int node2 = eNode[i];
				fprintf(file1, "%d %d %d %d %d %d\n", rand_num_cnt, node1,
						node2, 0, eDis[i], 100);
				break;
			}
			i++;
			rand_num_cnt++;
		}

	}

	fclose(file1);
}

void get_car_position_beijing() {
	FILE* file = fopen("ucarGPSBeijing.csv", "r");
	if (file == NULL) {
		cout << "file is NULL in read_cars()..." << endl;
        stop_here();
	}
	long id, timestamp, x3, y3;
	double lat, lon, speed;
	int angle;
	unordered_map<long, long> id_hash;
	vector<long> times;
	vector<CarPos> carposes;
	int id_cnt = 1;
	while (fscanf(file, "%ld, %ld, %lf, %lf, %lf, %d", &id, &timestamp, &lat,
			&lon, &speed, &angle) == 6) {
		x3 = (int) (lon * 1000000);
		y3 = (int) (lat * 1000000);
//		cout<<x3<<" "<<x_min<<" "<<x_len<<endl;
//		cout<<y3<<" "<<y_min<<" "<<x_len<<endl;
		int min_dist = INT_MAX;
		int min_node1 = 0;
		int min_node2 = 0;
		int cell_x = (x3 - x_min) / x_len;
		int cell_y = (y3 - y_min) / y_len;
//		cout<<"there "<<cell_x<<" "<<cell_y<<" "<<cells[cell_x][cell_y].size()<<endl;
		if (cell_x < 0 || cell_y < 0 || cell_x >= cell_num
				|| cell_y >= cell_num)
			continue;
		for (int t = 0; t < cells[cell_x][cell_y].size(); t++) {
//			cout<<"here"<<endl;
			int u = cells[cell_x][cell_y][t];
//			cout<<cells[cell_x][cell_y][t]<<endl;
			int x1 = x[u - 1];
			int y1 = y[u - 1];

//			cout<<x1<<" "<<y1<<endl;
			double d1 = sqrt(
					((x1 - x3) / 1000000.0) * ((x1 - x3) / 1000000.0)
							+ ((y1 - y3) / 1000000.0)
									* ((y1 - y3) / 1000000.0));

			for (int i = vStart[u]; i != 0; i = eNext[i]) {
				int v = eNode[i];
				int x2 = x[v - 1];
				int y2 = y[v - 1];
				double d2 = sqrt(
						((x1 - x2) / 1000000.0) * ((x1 - x2) / 1000000.0)
								+ ((y1 - y2) / 1000000.0)
										* ((y1 - y2) / 1000000.0));
				double d3 = sqrt(
						((x2 - x3) / 1000000.0) * ((x2 - x3) / 1000000.0)
								+ ((y2 - y3) / 1000000.0)
										* ((y2 - y3) / 1000000.0));
				double p = (d1 + d2 + d3) / 2;
				double p1 = p - d1;
				double p2 = p - d2;
				double p3 = p - d3;
				if (d3 < 0.00001)
					continue;
				double h = sqrt(p * p1 * p2 * p3) * 2.0 / d3;
				if (h < min_dist) {
					min_dist = h;
					min_node1 = u;
					min_node2 = v;
				}
			}
		}
		if (min_dist < INT_MAX) {
			double d1 = ((x3 - x[min_node1 - 1]) / 1000000.0)
					* ((x3 - x[min_node1 - 1]) / 1000000.0)
					+ ((y3 - y[min_node1 - 1]) / 1000000.0)
							* ((y3 - y[min_node1 - 1]) / 1000000.0);
			double d2 = ((x3 - x[min_node2 - 1]) / 1000000.0)
					* ((x3 - x[min_node2 - 1]) / 1000000.0)
					+ ((y3 - y[min_node2 - 1]) / 1000000.0)
							* ((y3 - y[min_node2 - 1]) / 1000000.0);
			double dis1 = sqrt(d1 - min_dist * min_dist);
			double dis2 = sqrt(d2 - min_dist * min_dist);
			double edge_len = get_distance_lat_lng(double(x[min_node1]),
					double(y[min_node1]), double(x[min_node2]),
					double(y[min_node2]), 1000000.0);
			if (id_hash.find(id) == id_hash.end()) {
				id_hash[id] = id_cnt;
				id_cnt++;
			}
			carposes.push_back(
					CarPos(id_hash[id], min_node1, min_node2,
							(int) (edge_len * dis1 / (dis1 + dis2)),
							(int) (edge_len * dis2 / (dis1 + dis2)),
							timestamp));
		}

	}
//	cout<<"size: "<<carposes.size()<<endl;
	sort(carposes.begin(), carposes.end(), CarPos::compare_time);
	string file_name1 = data_space;
	file_name1 = file_name1.append("car_positions.txt");
	FILE* file1 = fopen(file_name1.c_str(), "w");
	for (CarPos car : carposes) {
		fprintf(file1, "%d %d %d %d %d %ld\n", car.id, car.node1, car.node2,
				car.dis1, car.dis2, car.time);
	}
	fclose(file);
	fclose(file1);
}

void get_init_car_position_random(int num) {
	string file_name_input = data_space;
	file_name_input = file_name_input.append("car_positions.txt");
	FILE* file = fopen(file_name_input.c_str(), "r");
	if (file == NULL) {
		cout << "file is NULL in read_cars()..." << endl;
        stop_here();
	}
	int id, node1, node2, dis1, dis2;
	long timestamp;
	unordered_map<int, int> id_hash;
	vector<CarPos> carposes;
	int count = 0;
	while (fscanf(file, "%d %d %d %d %d %ld", &id, &node1, &node2, &dis1, &dis2,
			&timestamp) == 6) {
		if (id_hash.find(id) == id_hash.end()) {
			count++;
			id_hash[id] = 1;
			carposes.push_back(
					CarPos(count, node1, node2, dis1, dis2, timestamp));
			if (count >= num)
				break;
		}
	}
	cout << "size: " << carposes.size() << endl;
	sort(carposes.begin(), carposes.end(), CarPos::compare_time);
	string file_name_out = data_space;
	file_name_out = file_name_out.append("init_car_positions.txt");
	FILE* file1 = fopen(file_name_out.c_str(), "w");
	for (CarPos car : carposes) {
		fprintf(file1, "%d %d %d %d %d %ld\n", car.id, car.node1, car.node2,
				car.dis1, car.dis2, car.time);
	}
	fclose(file);
	fclose(file1);
}

void get_init_car_position() {
	string file_name_input = data_space;
	file_name_input = file_name_input.append("car_positions.txt");
	FILE* file = fopen(file_name_input.c_str(), "r");
	if (file == NULL) {
		cout << "file is NULL in read_cars()..." << endl;
        stop_here();
	}
	int id, node1, node2, dis1, dis2;
	long timestamp;
	unordered_map<int, int> id_hash;
	vector<CarPos> carposes;
	while (fscanf(file, "%d %d %d %d %d %ld", &id, &node1, &node2, &dis1, &dis2,
			&timestamp) == 6) {
		if (id_hash.find(id) == id_hash.end()) {
			id_hash[id] = 1;
			carposes.push_back(CarPos(id, node1, node2, dis1, dis2, timestamp));
		}
	}
	cout << "size: " << carposes.size() << endl;
	sort(carposes.begin(), carposes.end(), CarPos::compare_time);
	string file_name_out = data_space;
	file_name_out = file_name_out.append("init_car_positions.txt");
	FILE* file1 = fopen(file_name_out.c_str(), "w");
	for (CarPos car : carposes) {
		fprintf(file1, "%d %d %d %d %d %ld\n", car.id, car.node1, car.node2,
				car.dis1, car.dis2, car.time);
	}
	fclose(file);
	fclose(file1);
}

void read_car_positions(const char* path, vector<CarPos>& carposes,
		unordered_map<int, CarPos>& carposes_map,
		unordered_map<int, vector<PointCar> >& pointcar_hash) {
	carposes.clear();
	carposes_map.clear();
	pointcar_hash.clear();
	FILE* read_init_carpos = fopen(add_dataspace_to_name(path).c_str(), "r");
	if (read_init_carpos == NULL) {
		cout << "read init carpos is null" << endl;
        stop_here();
	}

	int id, node1, node2, dis1, dis2;
	long time;
	int init_car_cnt = 0;
	while (fscanf(read_init_carpos, "%d %d %d %d %d %ld", &id, &node1, &node2,
			&dis1, &dis2, &time) == 6) {
		carposes.push_back(CarPos(id, node1, node2, dis1, dis2, time));
		init_car_cnt++;
	}
	for (CarPos cp : carposes) {
		CarPos newCp(cp.id, cp.node1, cp.node2, cp.dis1, cp.dis2, cp.time);
		carposes_map[cp.id] = newCp;
	}

	FILE* read_all_carpos = fopen(
			add_dataspace_to_name("car_positions.txt").c_str(), "r");
	if (read_all_carpos == NULL) {
		cout << "read all carpos is null" << endl;
        stop_here();
	}
	while (fscanf(read_all_carpos, "%d %d %d %d %d %ld", &id, &node1, &node2,
			&dis1, &dis2, &time) == 6) {
		allcarposes.push_back(CarPos(id, node1, node2, dis1, dis2, time));
	}
	fclose(read_init_carpos);
	fclose(read_all_carpos);

	cout << "initially, we have " << init_car_cnt << " cars." << endl;

	for (CarPos cp : carposes) {
		if (pointcar_hash.find(cp.node1) == pointcar_hash.end()) {
			vector<PointCar> temp;
			pointcar_hash[cp.node1] = temp;
		}
		PointCar pointcar1(cp.id, cp.dis1);
		pointcar_hash[cp.node1].push_back(pointcar1);

		if (pointcar_hash.find(cp.node2) == pointcar_hash.end()) {
			vector<PointCar> temp;
			pointcar_hash[cp.node2] = temp;
		}
		PointCar pointcar2(cp.id, cp.dis2);
		pointcar_hash[cp.node2].push_back(pointcar2);
	}
	unordered_map<int, vector<PointCar> >::iterator iter;
	for (iter = pointcar_hash.begin(); iter != pointcar_hash.end(); iter++) {
		sort((iter->second).begin(), (iter->second).end(),
				PointCar::compare_dis);
	}
}
// Node Version: random
double make_random_nodes(double num_objs, int* car_nodes, vector<int>& objs) {
	objs.clear();
	unordered_map<int, int> hash_objs;
	memset(car_nodes, 0, sizeof(int) * (test_n + 1));
	for (int i = 0; i < num_objs;) {
		int node;
		node = get_rand_node();
		if (valid_node_hash.find(node) != valid_node_hash.end()) {
			if (hash_objs.find(node) == hash_objs.end()) {
				objs.push_back(node);
				hash_objs[node] = 1;
				car_nodes[node] = 1;
				i++;
			}
		}
	}
	return num_objs;
}


// Node version: shenzhou
double make_real_nodes(const char* path, int* car_nodes, vector<int>& objs) {
	objs.clear();
	double obj_num=0.0;
	unordered_map<int, int> hash_objs;
	memset(car_nodes, 0, sizeof(int) * (test_n + 1));
	string path_str(path, 100);
	FILE* file = fopen((input_parameters.input_data_dir+path_str).c_str(),"r");
	cout<<"reading objects from "<<(input_parameters.input_data_dir+path_str).c_str()<<endl;
	if(file == NULL){
		cout<<"file is NULL in make_shenzhou_nodes(), reading "<<(input_parameters.input_data_dir+path_str).c_str()<<endl;
        stop_here();
	}
	int node;
	while(fscanf(file, "%d", &node)==1){
		if (valid_node_hash.find(node) != valid_node_hash.end()) {
					if (hash_objs.find(node) == hash_objs.end()) {
						objs.push_back(node);
						hash_objs[node] = 1;
						car_nodes[node] = 1;
						obj_num+=1.0;
					}
				}
	}
	fclose(file);
	return obj_num;
}

vector<int> get_random_nodes(int num_objs) {
	vector<int> objs;
	unordered_map<int, int> hash_objs;
	for (int i = 0; i < num_objs;) {
		int node = get_rand_node();
		if (valid_node_hash.find(node) != valid_node_hash.end()) {
			if (hash_objs.find(node) == hash_objs.end()) {
				objs.push_back(node);
				hash_objs[node] = 1;
				i++;
			}
		}
	}
	return objs;
}

// Edge Version: random
void read_rand_car_positions(int rand_num, vector<CarPos>& carposes) {
	carposes.clear();
	unordered_map<int, int> hash_objs;

	int rand_num_cnt = 1;
	for (int i = 0; i < rand_num;) {
		int node1 = get_rand_node();
		if (valid_node_hash.find(node1) != valid_node_hash.end()) {
			if (hash_objs.find(node1) == hash_objs.end()) {

				for (int j = vStart[node1]; j != 0; j = eNext[j]) {

					int node2 = eNode[j];
					if (valid_node_hash.find(node2) != valid_node_hash.end()) {
						carposes.push_back(
								CarPos(rand_num_cnt, node1, node2, 0, eDis[j],
										100));
						i++;
						rand_num_cnt++;
						hash_objs[node1]=1;
						break;
					}
				}
			}

		}

	}
}
#endif /* SRC_OBJECTFUNCTIONS_H_ */
