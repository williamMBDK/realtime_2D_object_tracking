#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"
#include"util.cpp"
namespace SEG2{
  void evaluateRegions1(DATA::pixel_graph& g){
    vector<vector<int>> res (g.N);
    for(int i = 0; i < g.N; i++){
      int len = g.adjacency_list[i].size();
      int idx = -1;
      int mi = INT_MAX;
      for(int j = 0; j < len; j++){
        int nb = g.adjacency_list[i][j];
        int diff = DATA::squaredDifference(g.averagePixel[i], g.averagePixel[nb]);
        if(diff < mi){
          mi = diff;
          idx = j;
        }
      }
      vector<int>& best = g.averagePixel[g.adjacency_list[i][idx]];
      vector<int> t = UTIL::sub(best, g.averagePixel[i]);
      res[i] = UTIL::div_k(t, 2);
    }
    for(int i = 0; i < g.N; i++){
      g.averagePixel[i] = UTIL::add(g.averagePixel[i], res[i]);
    }
  }

  void evaluateRegions2(DATA::pixel_graph& g){
    vector<vector<double>> res (g.N);
    for(int i = 0; i < g.N; i++){
      int len = g.adjacency_list[i].size();
      int idx = -1;
      int mi = INT_MAX;
      for(int j = 0; j < len; j++){
        int nb = g.adjacency_list[i][j];
        int diff = DATA::squaredDifference(g.averagePixel[i], g.averagePixel[nb]);
        if(diff < mi){
          mi = diff;
          idx = j;
        }
      }
      vector<int>& best = g.averagePixel[g.adjacency_list[i][idx]];
      vector<int> t = UTIL::sub(best, g.averagePixel[i]);
      vector<double> r = UTIL::unitVector(t);
      res[i] = UTIL::mult_k(r, g.MAX_RGB / 10);
    }
    for(int i = 0; i < g.N; i++){
      g.averagePixel[i] = UTIL::add(g.averagePixel[i], res[i]);
    }
  }
}
