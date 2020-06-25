//
// Created by lostrong on 18/3/14.
//

#ifndef TOAIN_COMPUTEKSTAR_H
#define TOAIN_COMPUTEKSTAR_H

double calculate_C(int k,int h, int num_threads_update, double alpha){
    //calculate C(k,h)
    double res = 1.0;
    for(int i = 0;i<h;i++){
        res*=(k-i)*1.0/(h-i);

    }
    for(int i =0;i<h;i++){
        res*=alpha/num_threads_update;
    }
    for(int i = 0;i<k-h;i++){
        res*=(1.0-alpha/num_threads_update);
    }
    return res;
}

long double prob_at_least(int k, int h, int num_threads_update, double alpha){
    long double prob_at_least_h= 0.0;
    for(int i = h; i<=k;i++){
        prob_at_least_h+=calculate_C(k, i, num_threads_update, alpha);
    }
    return prob_at_least_h * num_threads_update;

}
int compute_k_star(int k, int num_threads_update, double alpha, double fail_p){
    for(int k_star=1;k_star<=k;k_star++){
        long double prob = prob_at_least(k, k_star, num_threads_update, alpha);
        cout<<k_star<<" "<<prob<<endl;
        if(prob<fail_p) return k_star;
    }
    return k;

}

#endif //TOAIN_COMPUTEKSTAR_H
