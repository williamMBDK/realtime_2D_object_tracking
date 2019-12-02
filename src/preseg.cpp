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

        I'(x, y, t) = max(
          |I(x + 1, y) - I(x, y)|,
          |I(x - 1, y) - I(x, y)|,
          |I(x, y + 1) - I(x, y)|,
          |I(x, y - 1) - I(x, y)|
        ) / k

        I.x.y_(n+1) = I_n + h/6*(k1 + 2*k2 + 2*k3 + k4)
        k1_x_y = f(t_n, I.x+1.y_n, I.x-1.y_n, I.x.y+1_n, I.x.y-1_n)
        k2_x_y = f(t_n + h/2, I.x+1.y_n + , I.x-1.y_n, I.x.y+1_n, I.x.y-1_n)

        run n iterations
          calculate k1 for all pixels
          calucate k2 for all pixels
          calculate k3 for all pixels
          calculate k4 for all pixels
          calucate the n'th value for all pixels
        the results of the last iteration is the resulting image

        constraints:
          1 <= x <= W
          1 <= Y <= H
          0 <= t <= infinite
          k e +R
    METHOD 1.1:
      Use variant 1
      Solve using dymanic programming and a discrete domain of the functions
      problemet med DP er at to pixels kan komme til at ligge og skifte mellem den ene og den anden. Den ene bliver til den anden farve og den anden blier den førstes hele tiden.
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

  // memory optimize dp
  void method1_3(IO::image& img){
    int N_MAX = 1;
    int dT = 0.1;
    vector<pair<int, int>> dirs = {
      {-1, 0}, {0, -1}, {1, 0}, {0, 1}
    };
    vector<vector<vector<int>>> dp (N_MAX, vector<vector<int>> (img.W, vector<int> (img.H)));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[0][i][j] = brightness(img.getPixel(i, j));
      }
    }
    for(int n = 1; n < N_MAX; n++){ // slight memory optimization on ks.
      vector<vector<vector<int>>> ks (4, vector<vector<int>> (img.W, vector<int> (img.H)));
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        int m = INT_MAX;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H) &&
             abs(dp[n-1][newI][newJ] - dp[n-1][i][j]) < abs(m)
          ){
            m = dp[n-1][newI][newJ] - dp[n-1][i][j];
          }
        }
        ks[0][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        int m = INT_MAX;
        int val = dp[n-1][i][j] + ks[0][i][j] * dT / 2.0;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            int nbVal = dp[n-1][newI][newJ] + ks[0][newI][newJ] * dT / 2.0;
            if(abs(nbVal - val) < abs(m)){
              m = nbVal - val;
            }
          }
        }
        ks[1][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        int m = INT_MAX;
        int val = dp[n-1][i][j] + ks[1][i][j] * dT / 2.0;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            int nbVal = dp[n-1][newI][newJ] + ks[1][newI][newJ] * dT / 2.0;
            if(abs(nbVal - val) < abs(m)){
              m = nbVal - val;
            }
          }
        }
        ks[2][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        int m = INT_MAX;
        int val = dp[n-1][i][j] + ks[2][i][j] * dT;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            int nbVal = dp[n-1][newI][newJ] + ks[2][newI][newJ] * dT;
            if(abs(nbVal - val) < abs(m)){
              m = nbVal - val;
            }
          }
        }
        ks[3][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        dp[n][i][j] = dp[n-1][i][j] + dT/6.0*(ks[0][i][j] + 2*ks[1][i][j] + 2*ks[2][i][j] + ks[3][i][j]);
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[N_MAX - 1][i][j] = min(dp[N_MAX - 1][i][j], img.MAX_RGB);
        dp[N_MAX - 1][i][j] = max(dp[N_MAX - 1][i][j], 0);
        img.setPixel(i, j, {dp[N_MAX - 1][i][j], dp[N_MAX - 1][i][j], dp[N_MAX - 1][i][j]});
      }
    }
  }

  // can be optimized by removing previous dp
  void method1_1(IO::image& img){
    int tMax = 200, K = 5;
    vector<vector<vector<int>>> dp (img.W, vector<vector<int>> (img.H, vector<int> (tMax)));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[i][j][0] = pow(brightness(img.getPixel(i, j)), 2);
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
        dp[i][j][tMax - 1] = sqrt(dp[i][j][tMax - 1]);
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
