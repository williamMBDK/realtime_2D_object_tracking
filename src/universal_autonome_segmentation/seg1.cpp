#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"

namespace SEG1{
  /*// modifies g with pixel merging using runge cutta
  void evaluateRegions(DATA::pixel_graph& g, int initialAmountOfSegments){
    int N_MAX = 20;
    //double dT = (double) g.N / 10.0 / (double)initialAmountOfSegments;
    double dT = 0.01;
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
  }*/
}
