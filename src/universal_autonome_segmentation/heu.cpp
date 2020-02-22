#pragma once
#include<bits/stdc++.h>
#include"util.cpp"
#include"data.cpp"
namespace HEU{
  int getSquaredDifference(graph& g, int node1, int node2){
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
    return diff;
  }

  double getPercentileSize(graph& g, int node){
    return (double)g.pixel_count[node] * 100.0 / (double)(g.W * g.H);
  }

  // assumes equal spread of pixels througout the conex hull or alpha shape which the d pixels cover.
  double getPseudoVariance(graph& g, ind node1, int node2){
    return 0.0;
  }
}
