#pragma once
#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

namespace PRESEG{
  /*
  METHOD 1: system of recursive partial differential equations
    2 variatiants of the system
      variant 1
        I(x, y, t) denotes the intensity (0-255) of the pixel at position (x, y) in an image of size (W, H) at time t.
        I'(x, y, t) = max(
          |I(x + 1, y, t - 1) - I(x, y, t - 1)|,
          |I(x - 1, y, t - 1) - I(x, y, t - 1)|,
          |I(x, y + 1, t - 1) - I(x, y, t - 1)|,
          |I(x, y - 1, t - 1) - I(x, y, t - 1)|
        ) / k
        constraints:
          1 <= x <= W
          1 <= Y <= H
          0 <= t <= infinite
          k e +R
    METHOD 1.1:
      Use variant 1
      Solve using dymanic programming and a discrete domain of the functions
    METHOD 1.2
      Use variant 1
      Solve using numerical methods (euler or runge-kutta*).
    METHOD 1.3:
      Use variant 2
      Solve using dymanic programming and a discrete domain of the functions
    METHOD 1.4
      Use variant 2
      Solve using numerical methods (euler or runge-kutta*).
  */

  int brightness(vector<int> pixel){
    return (int)(0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2]);
  }

  bool isInGrid(int x, int y, int W, int H){
    return x > 0 && x < W-1 && y > 0 && y < H-1; // does not hit boundary
  }

  // can be optimized by removing previous dp
  void method1_1(IO::image& img){
    int tMax = 200, K = 5;
    vector<vector<vector<int>>> dp (img.W, vector<vector<int>> (img.H, vector<int> (tMax)));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[i][j][0] = brightness(img.getPixel(i, j));
      }
    }
    vector<pair<int, int>> dirs = {
      {-1, 0}, {0, -1}, {1, 0}, {0, 1}
    };
    for(int t = 1; t < tMax; t++){
      for(int i = 0; i < img.W; i++){
        for(int j = 0; j < img.H; j++){
          dp[i][j][t] = dp[i][j][t-1];
          int m = INT_MAX;
          for(int k = 0; k < 4; k++){
            int newI = i + dirs[k].first;
            int newJ = j + dirs[k].second;
            if(isInGrid(newI, newJ, img.W, img.H) &&
               abs(dp[newI][newJ][t-1] - dp[i][j][t-1]) < abs(m)
            ){
              m = dp[newI][newJ][t-1] - dp[i][j][t-1];
            }
          }
          dp[i][j][t] += m;
        }
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[i][j][tMax - 1] = min(dp[i][j][tMax - 1], img.MAX_RGB);
        dp[i][j][tMax - 1] = max(dp[i][j][tMax - 1], 0);
        img.setPixel(i, j, {dp[i][j][tMax - 1], dp[i][j][tMax - 1], dp[i][j][tMax - 1]});
      }
    }
  }

  void experimential(IO::image& img){
    srand(time(NULL));
    int tMax = 200, K = 5;
    vector<vector<vector<int>>> dp (img.W, vector<vector<int>> (img.H, vector<int> (tMax)));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[i][j][0] = brightness(img.getPixel(i, j));
      }
    }
    vector<pair<int, int>> dirs = {
      {-1, 0}, {0, -1}, {1, 0}, {0, 1}
    };
    for(int t = 1; t < tMax; t++){
      for(int i = 0; i < img.W; i++){
        for(int j = 0; j < img.H; j++){
          dp[i][j][t] = dp[i][j][t-1];
          int m = INT_MAX;
          for(int k = 0; k < 4; k++){
            int newI = i + dirs[k].first;
            int newJ = j + dirs[k].second;
            if(isInGrid(newI, newJ, img.W, img.H) &&
               abs(dp[newI][newJ][t-1] - dp[i][j][t-1]) + ((rand() % 20) - 10) < abs(m)
            ){
              m = dp[newI][newJ][t-1] - dp[i][j][t-1];
            }
          }
          dp[i][j][t] += m;
        }
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[i][j][tMax - 1] = min(dp[i][j][tMax - 1], img.MAX_RGB);
        dp[i][j][tMax - 1] = max(dp[i][j][tMax - 1], 0);
        img.setPixel(i, j, {dp[i][j][tMax - 1], dp[i][j][tMax - 1], dp[i][j][tMax - 1]});
      }
    }
  }



}
