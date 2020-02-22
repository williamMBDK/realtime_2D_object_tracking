#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"
#include"util.cpp"
namespace SEG2{

  // returns boolean representing whether node1 should merge into node2
    // save is acceptable in dp.
    // factors
      // size of segment
      // is a segment acceptable
      //
  bool shouldPerformMerge(DATA::graph& g, int node1, int node2){
    if(g.pixel_count[node1] * 100 / (g.W * g.H) < 1) return true;
    vector<int> rgb1 = {g.mean_vector[node1][0], g.mean_vector[node1][1], g.mean_vector[node1][2]};
    vector<int> rgb2 = {g.mean_vector[node2][0], g.mean_vector[node2][1], g.mean_vector[node2][2]};
    vector<int> v1 = UTIL::mult_k(
      rgb1,
      g.pixel_count[node1]
    );
    vector<int> v2 = UTIL::mult_k(
      rgb2,
      g.pixel_count[node2]
    );
    vector<int> sum = UTIL::add(v1, v2);
    vector<int> merge = UTIL::div_k(sum, g.pixel_count[node1] + g.pixel_count[node2]);
    int diff = DATA::squaredDifference(
      rgb1,
      sum
    );
    if(diff > g.MAX_RGB * g.MAX_RGB){
      return false;
    }
    if((int)g.adjacency_list[node1].size() == 1 && diff > g.MAX_RGB * g.MAX_RGB / 2){
      return false;
    }
    if(g.pixel_count[node1] * 10 < g.pixel_count[node2]){
      return false;
    }

    // check if node1 is acceptable as a final segment, we call this acceptable.
    // false: if node1 is surrounded by node2 and node1 is acceptable.
    //        Even if node2 tries to merge towards node1 it will not be able to reach node1's rgbxy vector.
    // false: if node1 is acceptable and node1 gets lost in node2 - meaning the following:
    //          how big is the difference between the size of node1 and node2.
    //          if node1 is not similar enough to node2 and the size difference is big then false.
    // false: if node1 is acceptable and node2 is acceptable and their difference is big
    // true : otherwhise
    return true;
  }

  // returns a value representing the difference between node1 and node2
  int getDiff(DATA::graph& g, int node1, int node2){
    int res = DATA::squaredDifference(g.mean_vector[node1], g.mean_vector[node2]);
    return res;
  }

  // optimize: all difference are calulated twice.
  void evaluateRegions1(DATA::graph& g, bool withValidation){
    vector<vector<int>> res (g.N);
    for(int i = 0; i < g.N; i++){
      int len = g.adjacency_list[i].size();
      int idx = -1;
      int mi = INT_MAX;
      for(int j = 0; j < len; j++){
        int nb = g.adjacency_list[i][j];
        int diff = getDiff(g, i, nb);
        if(diff < mi){
          mi = diff;
          idx = j;
        }
      }
      if(withValidation && shouldPerformMerge(g, i, g.adjacency_list[i][idx])){
        vector<int>& best = g.mean_vector[g.adjacency_list[i][idx]];
        vector<int> t = UTIL::sub(best, g.mean_vector[i]);
        res[i] = UTIL::div_k(t, 2);
      }else if(!withValidation){
        vector<int>& best = g.mean_vector[g.adjacency_list[i][idx]];
        vector<int> t = UTIL::sub(best, g.mean_vector[i]);
        res[i] = UTIL::div_k(t, 2);
      }
    }
    for(int i = 0; i < g.N; i++){
      if(res[i].size() != 0) g.mean_vector[i] = UTIL::add(g.mean_vector[i], res[i]);
    }
  }

  /*double getDeltaPixel(double MAX_RGB, double N){
    return ((N/MAX_RGB*4)*(N/MAX_RGB*4)) + 1;
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
      res[i] = UTIL::mult_k(r, min(getDeltaPixel(g.MAX_RGB, g.N), UTIL::length(t) / 2));
    }
    for(int i = 0; i < g.N; i++){
      g.averagePixel[i] = UTIL::add(g.averagePixel[i], res[i]);
    }
  }*/
}
