/*
 * GenUtility.h
 *
 *  Created on: 2017年2月14日
 *      Author: lostrong
 */

#ifndef SRC_GENUTILITY_H_
#define SRC_GENUTILITY_H_
#include<string.h>
#include "../graph/Point.h"

string get_network_name_with_dataspace(string dataspace);
int get_level_num_with_dataspace();
inline double get_eu_dist(int s1, int s2);
void set_base_cell_num();
void test_valid_file(FILE* file, string file_name, string function_name);
void stop_here();


unsigned int get_rand_node() {
//	return my_rand() % (test_n) + 1;
    return rand()%(test_n)+1;
}
string get_network_name_with_dataspace(string dataspace) {
    if (data_space.compare("1") == 0)
        return "NY";
    if (data_space.compare("0") == 0)
        return "BJ";
    if (data_space.compare("7") == 0)
        return "NW";
    cout<<"error dataset"<<endl;
    stop_here();
    return "ERROR";
}

int get_level_num_with_dataspace() {
    if (data_space.compare("1") == 0)
        return 6;
    else if (data_space.compare("0") == 0)
        return 11; // for each smallest cell contains at most 50 nodes
    else if (data_space.compare("7") == 0)
        return 11;
    else if (data_space.compare("8") == 0)
        return 6;
    else if (data_space.compare("9") == 0)
        return 7;
    else if (data_space.compare("10") == 0)
        return 7;
    else if (data_space.compare("11") == 0)
        return 9;
    else if (data_space.compare("12") == 0)
        return 9;
    else if (data_space.compare("13") == 0)
        return 9;
    else if (data_space.compare("14") == 0)
        return 11;
    else
        return -1;
}


inline double get_eu_dist(int s1, int s2) {
    double x1 = x[s1 - 1];
    double y1 = y[s1 - 1];
    double x2 = x[s2 - 1];
    double y2 = y[s2 - 1];
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}





//=========================================================================
//= Multiplicative LCG for generating uniform(0.0, 1.0) random numbers    =
//=   - x_n = 7^5*x_(n-1)mod(2^31 - 1)                                    =
//=   - With x seeded to 1 the 10000th x value should be 1043618065       =
//=   - From R. Jain, "The Art of Computer Systems Performance Analysis," =
//=     John Wiley & Sons, 1991. (Page 443, Figure 26.2)                  =
//=========================================================================
double rand_val(int seed)
{
    const long  a =      16807;  // Multiplier
    const long  m = 2147483647;  // Modulus
    const long  q =     127773;  // m div a
    const long  r =       2836;  // m mod a
    static long x;               // Random int value
    long        x_div_q;         // x divided by q
    long        x_mod_q;         // x modulo q
    long        x_new;           // New x value

    // Set the seed if argument is non-zero and then return zero
    if (seed > 0)
    {
        x = seed;
        return(0.0);
    }

    // RNG using integer arithmetic
    x_div_q = x / q;
    x_mod_q = x % q;
    x_new = (a * x_mod_q) - (r * x_div_q);
    if (x_new > 0)
        x = x_new;
    else
        x = x_new + m;

    // Return a random value between 0.0 and 1.0
    return((double) x / m);
}


inline int get_valid_random_node(){
    int query;
    do {
         query = get_rand_node() % test_n;
    } while (valid_node_hash[query] == 0);
    return query;
}
void inline set_base_cell_num() {
    if (data_space.compare("0") == 0)
        cell_num = 8192;
    if (data_space.compare("1") == 0)
        cell_num = 256;
    if (data_space.compare("7") == 0)
        cell_num = 8192;
    else
        cell_num = 1024;
//		cell_num = 1024;				//leaf_num=256;
    cout << "cell_num: " << cell_num << endl;
}
string inline int_to_string(int int_value){
    stringstream ss;
    string str_value;
    ss<<int_value;
    ss>>str_value;
    return str_value;
}
void test_valid_file(FILE* file, string file_name, string function_name){
    if (file == NULL) {
        cout<<file_name << " is invalid in function "<<function_name<<endl;
        stop_here();
    }

}

void stop_here() {
    cout << "input a char to go on ..." << endl;
    char c;
    cin >> c;
}



double rad(double d) {
    return d * PI / 180.0;
}

double round(double d) {
    return floor(d + 0.5);
}

//m
double get_distance_lat_lng(double lat1, double lng1, double lat2, double lng2, double scale) {
    lat1 /= scale;
    lat2 /= scale;
    lng1 /= scale;
    lng2 /= scale;
    double radLat1 = rad(lat1);
    double radLat2 = rad(lat2);
    double a = radLat1 - radLat2;
    double b = rad(lng1) - rad(lng2);
    double s = 2 * asin(sqrt(pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2)));
    s = s * EARTH_RADIUS;
    s = round(s * 10000) / 10000;
    return s * 10000;
}



//After get_input
void init_cells(int n) {
    node_to_cell = new int[test_n + 1];
    cells = new vector<int>*[n];
    circle_cells = new vector<int>*[n];
    for (auto i = 0; i < n; ++i) {
        cells[i] = new vector<int> [n];
        circle_cells[i] = new vector<int> [n];
    }
}


void dispose_cells(int n) {

    for (auto i = 0; i < n; ++i) {
        delete[] cells[i];
        delete[] circle_cells[i];

    }
    delete[] cells;
    delete[] circle_cells;
    delete[] node_to_cell;
}

inline int get_location(int i, int j, int u) {
    int xz = cell_loc[u] / cell_num, yz = cell_loc[u] % cell_num;
    return circle_cell_loc[(xz - i + 2) * 5 + (yz - j + 2)][u];
}

int set_cells(int cell_num) {
    x_min = INT_MAX, y_min = INT_MAX, x_max = INT_MIN, y_max = INT_MIN, x_len =
                                                                                0, y_len = 0;
    Point *a = new Point[test_n + 1];
    vector<Point>* b = new vector<Point> [cell_num];
    cout << test_n << endl;
    for (int i = 0; i < test_n; ++i) {
        x_min = min(x_min, x[i]);
        x_max = max(x_max, x[i]);
        y_min = min(y_min, y[i]);
        y_max = max(y_max, y[i]);
    }
    int gap = (y_max - y_min) - (x_max - x_min);
    if (gap > 0) {
        //extend x_len
        x_min -= gap / 2;
        x_max += gap / 2;
    } else {
        gap = -gap;
        y_min -= gap / 2;
        y_max += gap / 2;
    }
    x_len = (x_max - x_min) / cell_num + 1;
    y_len = (y_max - y_min) / cell_num + 1;

    for (int i = 1; i < test_n + 1; ++i) {
        a[i].x = x[i - 1], a[i].y = y[i - 1];
        a[i].flag = i;
    }
    sort(a + 1, a + test_n + 1, Point::compare_x);
    for (int i = 1; i < test_n + 1; ++i) {
        int j = (a[i].x - x_min) / x_len;
        b[j].push_back(a[i]);
    }
    delete a;
    for (int i = 0; i < cell_num; ++i) {

        sort(b[i].begin(), b[i].end(), Point::compare_y);
        for (int j = 0; j < b[i].size(); ++j) {
            int k = (b[i][j].y - y_min) / y_len;
            cells[i][k].push_back(b[i][j].flag);
            cell_loc[b[i][j].flag] = i * cell_num + k;

            if (k >= cell_num) {
                cout << "k>=cell_num" << endl;
                while (true)
                    ;
            }
        }
    }
    const int xc[25] = { -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2,
                         -1, 0, 1, 2, -2, -1, 0, 1, 2 };
    const int yc[25] = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0,
                         1, 1, 1, 1, 1, 2, 2, 2, 2, 2 };

    for (int i = 0; i < cell_num; ++i) {

        for (int j = 0; j < cell_num; ++j) {
            for (int k = 0; k < 25; ++k) {
                int ux = i + xc[k], uy = j + yc[k];
                if (ux >= 0 && ux < cell_num && uy >= 0 && uy < cell_num) {
                    size_t len = circle_cells[i][j].size();
                    for (int l = 0; l < cells[ux][uy].size(); ++l) {
                        // the rank of the node in a specified 5*5 grid, the first dimension is in [0,24]
                        //	cout << (xc[k] + 2) * 5 + yc[k] + 2 << " " << cells[ux][uy][l] << " " << len + 1 << endl;
                        circle_cell_loc[(xc[k] + 2) * 5 + yc[k] + 2][cells[ux][uy][l]] =
                                len + l;

                        int loc = cell_loc[cells[ux][uy][l]];
                        if (loc / cell_num != ux || loc % cell_num != uy) {
                            cout << "asign false" << endl;
                            cout << ux << " " << uy << " " << l << endl;
                            while (true)
                                ;
                        }
                        //else
                        //	cout << "assign true" << endl;
                        int ts = get_location(i, j, cells[ux][uy][l]);
                        if (ts != len + l) {
                            cout << "inconsistent!" << endl;
                            while (true)
                                ;
                        }

                        //			get_location_array[i][j][cells[ux][uy][l]] = len + l;
                    }
                    //			circle_cells[i][j].insert(circle_cells[i][j].end(),
                    //					cells[ux][uy].begin(), cells[ux][uy].end());
                    for (auto it = cells[ux][uy].begin();
                         it != cells[ux][uy].end(); it++)
                        circle_cells[i][j].insert(circle_cells[i][j].end(),
                                                  *it);

                }
            }
        }
    }
    cell_loc[0] = -1;
    return 0;
}


#endif /* SRC_GENUTILITY_H_ */
