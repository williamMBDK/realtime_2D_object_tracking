#pragma once
#include<bits/stdc++.h>
using namespace std;
namespace UTIL{

  template <class number>
  int map(number num, number l1, number r1, number l2, number r2){
    return (num - l1) * (r2 - l2) / (r1 - l1) + l2;
  }

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
  double length (vector<number>&v){
    int l = (int)v.size();
    number s = 0;
    for(int i = 0; i < l; i++) s += v[i] * v[i];
    return sqrt((double)s);
  }

  template <class number>
  vector<double> unitVector(vector<number>&v){
    int l = (int)v.size();
    vector<double> res (l);
    double len = length(v);
    if(len == 0) return {0.0, 0.0, 0.0};
    for(int i = 0; i < l; i++) res[i] = v[i] / len;
    return res;
  }
}
