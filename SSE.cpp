#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <random>


using namespace std;

const int N = 10;
const int L = 10;
const double h = 4.0;
const int inf = 10000000;
mt19937 eng(11);
uniform_real_distribution<double> dist{0.0, 1.0};
uniform_int_distribution<> int_dist(0, N-1);
uniform_int_distribution<> int_L_dist(0, L-1);



void print(int F, vector<bool> gamma, vector<vector<bool>> &p, vector<vector<bool>> &q){
  for(int i = 0; i < N; i++){
    cout<<i<<" : ";
    if(gamma[i]==1) cout<<"- ";
    else cout<<"+ ";
    vector<string> outp;
    stringstream ss;
    for(int j = 0; j < N;j++){
      if(p[i][j]==1) {
        ss<<"X_"<<j;
        outp.push_back(ss.str());
      }
    }

    for(int j = 0; j < N;j++){
      stringstream ss;
      if(q[i][j]==1) {
        ss<<"Z_"<<j;
        outp.push_back(ss.str());
      }
    }
    for(int el_idx = 0; el_idx<outp.size(); el_idx++){
      cout<<outp[el_idx];
      if(el_idx!=outp.size()-1) cout<<" × ";
    }
    cout<<endl;
  }
  cout<<"F = "<<F<<endl;
  cout<<"-------------"<<endl;
  return;
}

void print_operators(vector<pair<int,int>> op_list){
  for(int i = op_list.size()-1; i>=0;i--){
    if(op_list[i].first==-1){
      cout<<"Id ";
    }
    if(op_list[i].first==1){
      cout<<"XMod_"<<op_list[i].second<<" ";
    }
    if(op_list[i].first==0){
      cout<<"CNOT_["<<op_list[i].second<<"] ";
    }
  }
  cout<<endl;
  return;
}

void apply(int &F, vector<bool> &gamma, vector<vector<bool>> &p, vector<vector<bool>> &q, int op, int idx){
  // op = 0: CX
  // op = 1: X_mod
  if(op == -1) return;
  if(op == 0){// CX operator
    int l = idx;
    int r = (idx+1)%N;

    for(int i = 0; i < N; i++){
      if(p[i][l] == 1) p[i][r] = !p[i][r];
      if(q[i][r] == 1) q[i][l] = !q[i][l];
    }
  }
  if(op == 1){//Xmod = (1+X)/2
    vector<int> anti_comm;
    for(int i = 0; i < N; i++){
      if(q[i][idx]==1) anti_comm.push_back(i);
    }
    if(anti_comm.size()==1){
      gamma[anti_comm[0]] = 0;
      for(int j = 0; j <N; j++){
        p[anti_comm[0]][j] = 0;
        q[anti_comm[0]][j] = 0;
      }
      p[anti_comm[0]][idx] = 1;
      F++;
    }
    if(anti_comm.size()>1){
      for(int i = 1; i < anti_comm.size(); i++){
        gamma[anti_comm[i]] = gamma[anti_comm[0]]^gamma[anti_comm[i]];
        for(int j = 0; j < N; j++){
          p[anti_comm[i]][j] = p[anti_comm[i]][j]^p[anti_comm[0]][j];
          q[anti_comm[i]][j] = q[anti_comm[i]][j]^q[anti_comm[0]][j];
        }
      }
      gamma[anti_comm[0]] = 0;
      for(int j = 0; j <N; j++){
        p[anti_comm[0]][j] = 0;
        q[anti_comm[0]][j] = 0;
      }
      p[anti_comm[0]][idx] = 1;
      F++;
    }
  }
  return;
}

void init_state(int &F, vector<bool> &gamma, vector<vector<bool>> &p, vector<vector<bool>> &q){
  //sets all qubits to the 0 state
  F = 0;
  for(int i = 0; i < N;i++){
    for(int j = 0; j < N; j++){
      p[i][j] = 0;
      q[i][j] = 0;
    }
    q[i][i] = 1;
    gamma[i] = 0;
  }
  return;
}

int get_not_Id(vector<bool>gamma,vector<vector<bool>>p,vector<vector<bool>>q){ 
  //finds the number of Identity operators among the generators
  int not_id = 0;
  for(int gen = 0; gen < N; gen++){
    bool is_id = true;
    for(int i =0 ; i < N; i++){
      if(p[gen][i]!=0 || q[gen][i]!=0){
        is_id = false;
        break;
      }
    }
    if(is_id){
      if(gamma[gen] == 1) return inf;
    }
    else not_id++;
  }
  return not_id;
}

int findIndp(int F,vector<bool>gamma,vector<vector<bool>>p,vector<vector<bool>>q, vector<bool>gamma_0){
  for(int i = 0; i < N;i++){
    for(int gen = 0; gen < N; gen++){
      if(q[gen][i]==1){
        gamma[gen] = gamma[gen] ^ gamma_0[i];
        q[gen][i] = 0;
      }
    }
  }
  int power = get_not_Id(gamma,p,q);
  if(power==inf) return inf;
  else return F+power;
}


double get_weight(int &F,vector<bool>&gamma,vector<vector<bool>>&p,vector<vector<bool>>&q, vector<bool>&gamma_0){
  int power = findIndp(F,gamma,p,q,gamma_0);
  if(power==inf) return 0.0;
  else return pow(0.5,((double) power/2));
}

void apply_all(int &F, vector<bool>&gamma, vector<vector<bool>>&p, vector<vector<bool>>&q, vector<pair<int,int>> op_list){
  for(int i = op_list.size()-1; i >=0 ; i--){
    apply(F,gamma,p,q,op_list[i].first,op_list[i].second);
  }
  return;
}

vector<bool> random_gamma(){
  vector<bool> gamma(N);
  for(int i = 0; i < N; i++){
    if(dist(eng)<0.5) gamma[i] = 0;
    else gamma[i] = 1;
  }
  return gamma;
}

double calc_weight(vector<bool> gamma_0,vector<pair<int,int>> op_list){
  int F;
  vector<vector<bool>> p(N, vector<bool> (N,1)),q(N, vector<bool> (N,1));
  vector<bool> gamma(N); 
  init_state(F, gamma,p,q);
  gamma = gamma_0;
  apply_all(F,gamma,p,q,op_list);
  double ans = get_weight(F,gamma,p,q,gamma_0);
  return ans;
}



int main(){
  vector<bool> gamma_0(N);
  vector<pair<int,int>> op_list(L);

  for(int idx = 0; idx < op_list.size(); idx++){
    op_list[idx] = {-1,-1};
  }
  for(double T = 10.1; T >= 0.1; T-=0.4){
    double avg_trace = 0.0;
    int n_ops = 0;
    for(int i = 0; i < op_list.size(); i++) if(op_list[i].first != -1) n_ops++;

    double n_avg = 0.0;
    int acc = 0;
    int rej = 0;
    const int rep_tot = 10000;
    const int rep_therm = rep_tot*1/2;
    double w_cur = calc_weight(gamma_0,op_list);
    for(int rep = 0; rep < rep_tot; rep++){
      //////////////////
      // State Update //
      //////////////////

      vector<bool> gamma_0_new = random_gamma();
      double w_2 = calc_weight(gamma_0_new,op_list);
      if(dist(eng) < w_2/w_cur){ //unless w_2 = 0, we have w_1 = w_2; so always make the transition in latter case
        gamma_0 = gamma_0_new;
        w_cur = w_2;
        acc++;
      }
      else {
        rej++;
      }

      
      /////////////////////
      // Operator Update //
      /////////////////////
      for(int i = 0; i < L; i++){
        if (op_list[i].first==-1){
          double p_op_type = h/(h+1.0);
          if(dist(eng) < p_op_type){
            int idx = int_dist(eng);
            
            double prob = ((double) (N));

            op_list[i].first = 1; op_list[i].second = idx;
            double w_CX = calc_weight(gamma_0, op_list);
            prob *= w_CX/w_cur;
            prob *= (h+1.0);
            prob /= T;
            prob /= ((double) (L - n_ops));
            n_ops++;
            if(prob < dist(eng)){
              
              op_list[i].first = -1; op_list[i].second = -1;
              n_ops--;
            }
            else{
              w_cur = w_CX;
            }
          }
          else{
            int idx = int_dist(eng);
            
            double prob = ((double) (N));

            op_list[i].first = 0; op_list[i].second = idx;
            double w_CX = calc_weight(gamma_0, op_list);
            prob *= w_CX/w_cur;
            prob *= (h+1.0);
            prob /= T;
            prob /= ((double) (L - n_ops));
            n_ops++;
            if(prob < dist(eng)){
              
              op_list[i].first = -1; op_list[i].second = -1;
              n_ops--;
            }
            else{
              w_cur = w_CX;
            }
          }
        }
        else if (op_list[i].first!=-1){
          int type = op_list[i].first;
          int idx = op_list[i].second;
          op_list[i].first = -1; op_list[i].second = -1;
          double w_I = calc_weight(gamma_0, op_list);


          double prob = w_I/w_cur;
          prob *= ((double) (L- n_ops+1));
          prob *= T;
          prob /= ((double) (N));
          prob /= (h+1.0);

          n_ops--;

          if(prob < dist(eng)){
            op_list[i].first = type; op_list[i].second = idx;
            n_ops++;
          }
          else{
            w_cur = w_I;
          }
        }
      }
      if(rep_therm < rep) n_avg += ((double) n_ops);
    }
    n_avg /= ((double) (rep_tot-rep_therm));
    avg_trace /= ((double) rep_tot);
    cout<<"{"<<T<<", "<<-T*n_avg<<"}, ";
    cout.flush();
  }
  return 0;
}