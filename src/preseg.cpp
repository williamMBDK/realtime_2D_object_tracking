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
          I(x + 1, y, t - 1),
          I(x - 1, y, t - 1),
          I(x, y + 1, t - 1),
          I(x, y - 1, t - 1)
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

  // can be optimized by removing previous dp
  void method1_1(IO::image& img){
    int tMax = 5, k = 5;
    vector<vector<vector<int>>> dp (img.W + 2, vector<vector<int>> (img.H + 2, vector<int> (tMax)));
    for(int i = 1; i <= img.W; i++){
      for(int j = 1; j <= img.H; j++){
        dp[i][j][0] = brightness(img.getPixel(i-1, j-1));
      }
    }
    for(int t = 1; t < tMax; t++){
      for(int i = 1; i <= img.W; i++){
        for(int j = 1; j <= img.H; j++){
          dp[i][j][t] = dp[i][j][t-1] + max(
            max(
              dp[i-1][j][t-1],
              dp[i+1][j][t-1]
            ),
            max(
              dp[i][j-1][t-1],
              dp[i][j+1][t-1]
            )
          ) / k;
        }
      }
    }
    for(int i = 1; i <= img.W; i++){
      for(int j = 1; j <= img.H; j++){
        img.setPixel(i - 1, j - 1, {dp[i][j][tMax - 1], dp[i][j][tMax - 1], dp[i][j][tMax - 1]});
      }
    }
  }
}
