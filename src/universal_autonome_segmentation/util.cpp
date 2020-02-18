#pragma once
#include<bits/stdc++.h>
using namespace std;
namespace UTIL{
  template <class number>
  void printVector(vector<number>& v){
    for(int i = 0; i < (int)v.size(); i++){
      cout << v[i] << " ";
    }
    cout << endl;
  }

  template <class number1, class number2>
  vector<number1> mult_k(vector<number1>&v, number2 k){
    int l = (int)v.size();
    vector<number1> res (l);
    for(int i = 0; i < l; i++) res[i] = v[i] * k;
    return res;
  }

  template <class number1, class number2>
  vector<number1> div_k(vector<number1>&v, number2 k){
    int l = (int)v.size();
    vector<number1> res (l);
    for(int i = 0; i < l; i++) res[i] = v[i] / k;
    return res;
  }

  template <class number>
  vector<number> sub(vector<number>&v1, vector<number>&v2){
    int l = (int)v1.size();
    vector<number> res (l);
    for(int i = 0; i < l; i++) res[i] = v1[i] - v2[i];
    return res;
  }

  template <class number1, class number2>
  vector<number1> add(vector<number1>&v1, vector<number2>&v2){
    int l = (int)v1.size();
    vector<number1> res (l);
    for(int i = 0; i < l; i++) res[i] = v1[i] + v2[i];
    return res;
  }

  template <class number>
  vector<double> unitVector(vector<number>&a1){
    return {};
  }
}
