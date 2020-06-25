#ifndef SRC_DATAFUNCTIONS_H_
#define SRC_DATAFUNCTIONS_H_

#include <vector>
#include "../para/GlobalVariables.h"
#include "../object/CarPos.h"
#include "../object/ObjectFunctions.h"
#include "../object/UCarInfo.h"
#include "../object/CarPeriodLocation.h"

#include <strstream>
using namespace std;

vector<string> split(string str, string separator) {
	vector<string> result;
	int cutAt;
	while ((cutAt = str.find_first_of(separator)) != str.npos) {
		if (cutAt > 0) {
			result.push_back(str.substr(0, cutAt));
		}
		str = str.substr(cutAt + 1);
	}
	if (str.length() > 0) {
		result.push_back(str);
	}
	return result;
}

int get_nearest_node(double coor_x, double coor_y) {
	int returned_node = -1;
	//	double value = numeric_limits<double>::max();
	double value = 100000000000000.0;
	for (int i = 1; i < test_n + 1; i++) {
		double x1 = x[i - 1] / 1000000.0;
		double y1 = y[i - 1] / 1000000.0;
		//		cout << x1 << " " << y1 << " " << coor_x << " " << coor_y << endl;
//		double eu_dist = sqrt((coor_x - x1) * (coor_x - x1) + (coor_y - y1) * (coor_y - y1));
		double eu_dist = (coor_x - x1) * (coor_x - x1) + (coor_y - y1) * (coor_y - y1);
		if (eu_dist < value) {
			value = eu_dist;
			returned_node = i;
		}
	}
	return returned_node;
}


void gen_fla_poi(char const* outPath) {
	int FLA_avg_edge = 4838;
	FILE* fInput = fopen("fla_poi.txt", "r");
	if (fInput == NULL) {
		cout << "cannot find fla_poi.txt" << endl;
		return;
	}
	double coor_x, coor_y;
	vector<int> nNodes;
	// note reverse
	while (fscanf(fInput, "%lf, %lf", &coor_y, &coor_x) == 2) {
		int node = get_nearest_node(coor_x, coor_y);
		nNodes.push_back(node);
	}
	FILE* fOutput = fopen(outPath, "w");
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			if (j == i)
				continue;
			for (int k = 0; k < 10; k++) {
				if (k == j || k == i)
					continue;
				fprintf(fOutput, "%d %d %d %d 0\n", nNodes[i], nNodes[j], nNodes[k], FLA_avg_edge * 5);
				fprintf(fOutput, "%d %d %d %d 0\n", nNodes[i], nNodes[j], nNodes[k], FLA_avg_edge * 10);
				fprintf(fOutput, "%d %d %d %d 0\n", nNodes[i], nNodes[j], nNodes[k], FLA_avg_edge * 15);
				fprintf(fOutput, "%d %d %d %d 0\n", nNodes[i], nNodes[j], nNodes[k], FLA_avg_edge * 20);
			}
		}
	}
	fclose(fInput);
	fclose(fOutput);
}


double get_ucar_object_vector(vector<vector<CarPeriodLocation> >& car_period_location_list, int* object_nodes_hash,
                                   vector<int> &objects_vector, int period_num, vector<int>& indexes){
    double valid_car_num = 0.0;
    objects_vector.clear();
    memset(object_nodes_hash, 0, sizeof(int) * (test_n + 1));
    for(int i = 0;i < car_period_location_list.size(); i++){
        if(valid_cars[i]) {
            valid_car_num+=1.0;
            vector<CarPeriodLocation> &list_i = car_period_location_list[i];
            for (int j = indexes[i]; j < list_i.size(); j++) {
                if (list_i[j].period_num >= period_num) {
                    indexes[i] = j;
                    objects_vector.push_back(list_i[j].nearest_node_id);
                    object_nodes_hash[list_i[j].nearest_node_id]=1;
                    break;
                }
            }
        }
    }
    return valid_car_num;
}
string intToStr(int int_value) {
	stringstream ss;
	ss << int_value;
	return ss.str();
}

string doubleToStr(int double_value) {
	stringstream ss;
	ss << double_value;
	return ss.str();
}
void gen_input_data_ran(const char *network_name, const char *density_str, double density_value,
                        const char *method_type, int obj_type, int num, int test_n) {
	for (int i = 0; i < num; i++) {
		string path = network_name;
		path = path.append("_");
		path = path.append(method_type);
		path = path.append("_");
		path = path.append(density_str);
		path = path.append("_");
		path = path.append(intToStr(obj_type));
		if (i < 10)
			path = path.append("_00");
		else if (i < 100)
			path = path.append("_0");
		path = path.append(intToStr(i));
		path = path.append(".txt");
		FILE* file = fopen(path.c_str(), "w");
		vector<int> objs = get_random_nodes((int) (test_n * density_value));
		fprintf(file, "%s %s %s %d %d\n", network_name, method_type, density_str, obj_type,
				(int) (test_n * density_value));
		for (int obj : objs) {
			fprintf(file, "%d\n", obj);
		}
		fclose(file);
	}
    stop_here();

}


vector<int> read_queries(const char* network_name, int num) {
	vector<int> queries;
	string path = network_name;
	path = path.append("_");
	path = path.append(intToStr(num));
	path = path.append(".txt");
	FILE* file = fopen((input_parameters.input_data_dir+path).c_str(), "r");
	cout << "query read from " << path.c_str() << endl;
	if (file == NULL) {
		cout << "file == NULL" << endl;
        stop_here();
	}
	int query;
	char s1[20], s2[20];
	int a, b, c;
	fscanf(file, "%s %s %d %d %d", s1, s2, &a, &b, &c);
	while (fscanf(file, "%d", &query) == 1) {
		queries.push_back(query);
	}
	return queries;
}

string make_out_file_name(string method_name, int k, long double obj_ratio, string postfix){
	string mode_name;
	string obj_name;
    if(test_parameters.arrival_model==1)
		mode_name="uber";
	else
		mode_name="poke";
	if(test_parameters.data==0)
		obj_name="ran";
	else
		obj_name="real";

    stringstream ss;
    string k_str;
	if(k>0) {
		ss << k;
		ss >> k_str;
	}

    string obj_ratio_str="";
	if(obj_ratio>0) obj_ratio_str=std::to_string(obj_ratio);
	string res= mode_name + "_" + network_name + "_" + obj_name + "_" + method_name;

	if(test_parameters.data==0) {
		if(obj_ratio>0)
			res=res+ "_" + obj_ratio_str;

	}

	if(k>0)
		res=res+ "_" + k_str;
    if(ZIPF_OBJECT)
        res=res+ "_" + "zipfobj";
    if(ZIPF_QUERY)
        res = res+"_zipfq";

    string alpha_str;

    if(ZIPF_OBJECT || ZIPF_QUERY){
        if(ALPHA==1.0)
            alpha_str="1_0";
        if(ALPHA==1.5)
            alpha_str="1_5";
        if(ALPHA==2.0)
            alpha_str="2_0";
        if(ALPHA==2.5)
            alpha_str="2_5";
        if(ALPHA==3.0)
            alpha_str="3_0";
        if(ALPHA==1.5)
            alpha_str="1_5";
        res = res+"_"+alpha_str;
    }

	return res+"_" + postfix + ".txt";
}

double n_hundred(int n)
{
    double sum = 1.0;
    for(int i = 1;i <= n;++i)
        sum = sum *10;
    return sum;
}
double process(double t,int n)
{
    //如果要求保留t小数点后n位数字的话
    int ival = (int)(t * n_hundred(n));
    //看小数点后第n+1位数字是否大于5来进行四舍五入
    int temp = (int)(t * n_hundred(n+1))%10;
    if(temp >= 5)
        ival = ival + 1;
    double dval = ival/n_hundred(n);
    return dval;
}


double make_obj_data(double &obj_num, int *object_nodes_hash, vector<int> &objects_vector, int period_num = 0) {
    if (test_parameters.data == DATA_BJ_REAL) {// for UCAR data

        obj_num = make_real_nodes("beijing_car_positions_node.txt", object_nodes_hash, objects_vector);


    } else if (test_parameters.data == DATA_NW_REAL) {// for NW POI
        obj_num = make_real_nodes("NW_POI.txt", object_nodes_hash, objects_vector);
    } else {

            obj_num = make_random_nodes((obj_num), object_nodes_hash, objects_vector);

    }
    cout << (obj_num) << " objects" << endl;
    return obj_num;

}

int getNearestNode(double coor_x, double coor_y) {
    int returned_node = -1;
    //	double value = numeric_limits<double>::max();
    double value = 1000000000.0;
    for (int i = 1; i < test_n + 1; i++) {
        double x1 = x[i - 1] / 1000000.0;
        double y1 = y[i - 1] / 1000000.0;
        double eu_dist =
                (coor_x - x1) * (coor_x - x1) + (coor_y - y1) * (coor_y - y1);
        if (eu_dist < value) {
            value = eu_dist;
            returned_node = i;
        }
    }
    return returned_node;
}



#endif /* SRC_DATAFUNCTIONS_H_ */