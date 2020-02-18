#pragma once
#include<bits/stdc++.h>
using namespace std;
namespace UTIL{
  void printVector(vector<int>& v){
    for(int i = 0; i < (int)v.size(); i++){
      cout << v[i] << " ";
    }
    cout << endl;
  }

  vector<int> mult_k(vector<int>&v, int k){
    int l = (int)v.size();
    vector<int> res (l);
    for(int i = 0; i < l; i++) res[i] = v[i] * k;
    return res;
  }

  vector<int> div_k(vector<int>&v, int k){
    int l = (int)v.size();
    vector<int> res (l);
    for(int i = 0; i < l; i++) res[i] = v[i] / k;
    return res;
  }

  vector<int> sub(vector<int>&v1, vector<int>&v2){
    int l = (int)v1.size();
    vector<int> res (l);
    for(int i = 0; i < l; i++) res[i] = v1[i] - v2[i];
    return res;
  }

  vector<int> add(vector<int>&v1, vector<int>&v2){
    int l = (int)v1.size();
    vector<int> res (l);
    for(int i = 0; i < l; i++) res[i] = v1[i] + v2[i];
    return res;
  }
}
