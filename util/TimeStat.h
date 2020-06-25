/*
 * TimeStat.h
 *
 *  Created on: 2017年1月17日
 *      Author: lostrong
 */

#ifndef SRC_TIMESTAT_H_
#define SRC_TIMESTAT_H_
#include<vector>

//vector<int> accumulated_comb_ids;
struct times {
    int k;
    double obj_ratio;
    string network_name;
    string method_name;
    string obj_name;
	double query;
	double update;
	int init;
	double query_variance;
	double update_variance;
    int scob_configuration_id;
    int is_valid=1;
};

struct simulate_times{
    int k;
    double obj_ratio;
    string network_name;
    string method_name;
    string obj_name;
    vector<vector<int> > query_times;
    vector<int> update_times;
    int mode; //1: query_first 0: fifo

};
double get_mean(vector<long>& times) {
	double mean = 0.0;
	for (long val : times) {
		mean += val;
	}
	if(times.size()>0)
		return mean / times.size();
	else
		return 0.0;
}
double get_var(vector<long>& times) {
	double mean = get_mean(times);
	double var = 0.0;
	for (long val : times) {
		var += (val - mean) * (val - mean);
	}
	if(times.size()>0)
		return var / times.size();
	else
		return 0.0;
}

double get_mean_int(vector<int>& times) {
    double mean = 0.0;
    for (int val : times) {
        mean += val;
    }
    return mean / times.size();
}
double get_var_int(vector<int>& times) {
    double mean = get_mean_int(times);
    double var = 0.0;
    for (int val : times) {
        var += (val - mean) * (val - mean);
    }
    return var / times.size();
}

double get_mean_double(vector<double>& times) {
    double mean = 0.0;
    for (double val : times) {
        mean += val;
    }
    return mean / times.size();
}
double get_var_double(vector<double>& times) {
    double mean = get_mean_double(times);
    double var = 0.0;
    for (double val : times) {
        var += (val - mean) * (val - mean);
    }
    return var / times.size();
}

void estimate_hito(FILE* file, vector<long>& cgraph_each_query_time, int sep) {
//	cout<<"start estimation"<<endl;
	long max = -1;
	long min = INT_MAX;
	for (long val : cgraph_each_query_time) {
		if (val > max)
			max = val;
		if (val < min)
			min = val;
	}
//	cout<<"start estimation"<<endl;
	double gap = (max - min) * 1.0 / sep;
//	cout<<gap<<endl;
	int * cnt = new int[sep + 1];
	memset(cnt, 0, sizeof(int) * (sep + 1));
	for (long val : cgraph_each_query_time) {
//		cout<<(val-min)/gap<<endl;
		cnt[(int) ((val - min) / gap)]++;
	}
	for (int i = 0; i < sep; i++) {
		fprintf(file, "%lf %lf %d\n", min + gap * i, min + gap * (i + 1), cnt[i]);
	}
	fprintf(file, "\n\n");
	delete[] cnt;
}



#endif /* SRC_TIMESTAT_H_ */
