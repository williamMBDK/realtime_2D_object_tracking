#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"

namespace SEG1{
  // modifies g with pixel merging using runge cutta
  void evaluateRegions(DATA::pixel_graph& g, int initialAmountOfSegments){
    int N_MAX = 20;
    //double dT = (double) g.N / 10.0 / (double)initialAmountOfSegments;
    double dT = 0.1;
    vector<vector<double>> dp (N_MAX, vector<double> (g.N, 0.0));
    for(int i = 0; i < g.N; i++){
      dp[0][i] = (double)pow(DATA::brightness(g.averagePixel[i]), 2);
    }
    for(int n = 1; n < N_MAX; n++){
      vector<vector<double>> ks (4, vector<double> (g.N));
      for(int i = 0; i < g.N; i++){
        double m = DBL_MAX;
        double val = dp[n-1][i];
        for(int j = 0; j < (int)g.adjacency_list[i].size(); j++){
          double nbVal = dp[n-1][g.adjacency_list[i][j]];
          if(abs(nbVal - val) < abs(m)){
            m = nbVal - val;
          }
        }
        ks[0][i] = m;
      }
      for(int i = 0; i < g.N; i++){
        double m = DBL_MAX;
        double val = dp[n-1][i] + ks[0][i] * dT / 2.0;
        for(int j = 0; j < (int)g.adjacency_list[i].size(); j++){
          double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[0][g.adjacency_list[i][j]] * dT / 2.0;
          if(abs(nbVal - val) < abs(m)){
            m = nbVal - val;
          }
        }
        ks[1][i] = m;
      }
      for(int i = 0; i < g.N; i++){
        double m = DBL_MAX;
        double val = dp[n-1][i] + ks[1][i] * dT / 2.0;
        for(int j = 0; j < (int)g.adjacency_list[i].size(); j++){
          double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[1][g.adjacency_list[i][j]] * dT / 2.0;
          if(abs(nbVal - val) < abs(m)){
            m = nbVal - val;
          }
        }
        ks[2][i] = m;
      }
      for(int i = 0; i < g.N; i++){
        double m = DBL_MAX;
        double val = dp[n-1][i] + ks[2][i] * dT;
        for(int j = 0; j < (int)g.adjacency_list[i].size(); j++){
          double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[2][g.adjacency_list[i][j]] * dT;
          if(abs(nbVal - val) < abs(m)){
            m = nbVal - val;
          }
        }
        ks[3][i] = m;
      }
      for(int i = 0; i < g.N; i++){
        dp[n][i] = dp[n-1][i] + dT/6.0*(ks[0][i] + 2*ks[1][i] + 2*ks[2][i] + ks[3][i]);
      }
    }
    for(int i = 0; i < g.N; i++){
      dp[N_MAX - 1][i] = sqrt(dp[N_MAX - 1][i]);
      g.averagePixel[i] = {(int)dp[N_MAX - 1][i], (int)dp[N_MAX - 1][i], (int)dp[N_MAX - 1][i]};
    }
  }

  // returns a pixelgraph where regions of same color in g has been merged
  DATA::pixel_graph mergeRegions(DATA::pixel_graph& g){
    DATA::pixel_graph res;
    int component = -1;
    vector<int> components (g.N, -1);
    vector<vector<int>> pixelMap (g.W, vector<int> (g.H));
    for(int i = 0; i < g.N; i++){
      if(components[i] == -1){
        queue<int> q; q.push(i);
        component++;
        while(!q.empty()){
          int node = q.front(); q.pop();
          if(components[node] != -1) continue;
          components[node] = component;
          for(int j = 0; j < (int)g.pixels[node].size(); j++){
            pixelMap[g.pixels[node][j].first][g.pixels[node][j].second] = component;
          }
          for(int j = 0; j < (int)g.adjacency_list[node].size(); j++){
            if(DATA::brightness(g.averagePixel[g.adjacency_list[node][j]]) == DATA::brightness(g.averagePixel[node])){
              q.push(g.adjacency_list[node][j]);
            }
          }
        }
      }
    }
    res.N = component + 1;
    res.MAX_RGB = g.MAX_RGB; res.W = g.W; res.H = g.H;
    res.adjacency_list = vector<vector<int>> (res.N);
    res.pixels = vector<vector<pair<int, int>>> (res.N);
    res.averagePixel = vector<vector<int>> (res.N, {0, 0, 0});
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, -1},
      {0, 1}
    };
    vector<unordered_set<int>> nbs (res.N);
    for(int i = 0; i < g.N; i++){
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        pair<int, int> node = g.pixels[i][j];
        res.pixels[components[i]].push_back(node);
        res.averagePixel[components[i]][0] += g.averagePixel[i][0];
        res.averagePixel[components[i]][1] += g.averagePixel[i][1];
        res.averagePixel[components[i]][2] += g.averagePixel[i][2];
        for(int l = 0; l < 4; l++){
          int newX = node.first + dirs[l].first;
          int newY = node.second + dirs[l].second;
          if(DATA::isInSideGrid(newX, newY, g.W, g.H) && pixelMap[newX][newY] != components[i]){
            nbs[components[i]].insert(pixelMap[newX][newY]);
          }
        }
      }
    }
    for(int i = 0; i < res.N; i++){
      res.adjacency_list[i] = vector<int> (nbs[i].size());
      int idx = 0;
      for(int nb : nbs[i]){
        res.adjacency_list[i][idx++] = nb;
      }
      res.averagePixel[i][0] /= res.pixels[i].size();
      res.averagePixel[i][1] /= res.pixels[i].size();
      res.averagePixel[i][2] /= res.pixels[i].size();
    }
    return res;
  }
}
