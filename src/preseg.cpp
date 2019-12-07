#pragma once
#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

namespace PRESEG{
  /*
  METHOD 1: system of recursive partial differential equations
    God hvis man ikke ved hvor mange segmenter der skal / kan være.
    det er på en graf af pixels.
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
    METHOD 1.5
      udvidelse. kør en metode flere gange og byg en ny graf hver gang, der kigger på alle naboer til en knude på mappet. O(logn * n * STOR_CONSTANT) antallet af knuder vil altid mindst halveres.
      Når den nye graf bygges svarer det til et graphcut.
    METHOD 2:
      k-means clusering with k-means++
      dimensions:
        x,y,r,g,b
  */

  int brightness(vector<int> pixel){
    return (int)(0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2]);
  }

  bool isInGrid(int x, int y, int W, int H){
    return x > 0 && x < W-1 && y > 0 && y < H-1; // does not hit boundary
  }

  // memory optimize dp
  void method1_3(IO::image& img){
    int N_MAX = 2000;
    double dT = 0.1;
    vector<pair<int, int>> dirs = {
      {-1, 0}, {0, -1}, {1, 0}, {0, 1}
    };
    vector<vector<vector<double>>> dp (N_MAX, vector<vector<double>> (img.W, vector<double> (img.H)));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        dp[0][i][j] = (double)pow(brightness(img.getPixel(i, j)), 2);
      }
    }
    for(int n = 1; n < N_MAX; n++){ // slight memory optimization on ks.
      vector<vector<vector<double>>> ks (4, vector<vector<double>> (img.W, vector<double> (img.H)));
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        double m = DBL_MAX;
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
        double m = DBL_MAX;
        double val = dp[n-1][i][j] + ks[0][i][j] * dT / 2.0;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            double nbVal = dp[n-1][newI][newJ] + ks[0][newI][newJ] * dT / 2.0;
            if(abs(nbVal - val) < abs(m)){
              m = nbVal - val;
            }
          }
        }
        ks[1][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        double m = DBL_MAX;
        double val = dp[n-1][i][j] + ks[1][i][j] * dT / 2.0;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            double nbVal = dp[n-1][newI][newJ] + ks[1][newI][newJ] * dT / 2.0;
            if(abs(nbVal - val) < abs(m)){
              m = nbVal - val;
            }
          }
        }
        ks[2][i][j] = m;
      }
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        double m = DBL_MAX;
        double val = dp[n-1][i][j] + ks[2][i][j] * dT;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInGrid(newI, newJ, img.W, img.H)){
            double nbVal = dp[n-1][newI][newJ] + ks[2][newI][newJ] * dT;
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
        dp[N_MAX - 1][i][j] = sqrt(dp[N_MAX - 1][i][j]);
        dp[N_MAX - 1][i][j] = min(dp[N_MAX - 1][i][j], (double)img.MAX_RGB);
        dp[N_MAX - 1][i][j] = max(dp[N_MAX - 1][i][j], 0.0);
        img.setPixel(i, j, {(int)dp[N_MAX - 1][i][j], (int)dp[N_MAX - 1][i][j], (int)dp[N_MAX - 1][i][j]});
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

  int dist(vector<int> &v1, vector<int> &v2){
    int s = 0;
    for(int i = 0; i < v1.size(); i++){
      s += (v1[i] - v2[i])*(v1[i] - v2[i]);
    }
    return (int)sqrt(s);
  }

  int mapRange(int v, int r1, int r2){
    return (int)((double)v/(double)r1*(double)r2);
  }

  void method2(IO::image& img){
    srand(time(NULL)); // seeding rng
    int thres = 10;
    int D = 5; // dimensions
    int MAX_VALUE = 300;
    int POSWEIGTH = 5;
    vector<vector<int>> clusters = vector<vector<int>> ();
    vector<vector<int>> clusterColors = vector<vector<int>> ();
    int size = 30;
    for(int i = size/2; i < img.W; i+=size){
      for(int j = size/2; j < img.H; j+=size){
        vector<int> rgb = vector<int> (3);
        int c = 0;
        for(int k = -1; k < 2; k++){
          for(int l = -1; l < 2; l++){
            int newI = i + k;
            int newJ = j + l;
            if(isInGrid(newI, newJ, img.W, img.H)){
              vector<int> pixel = img.getPixel(newI, newJ);
              rgb[0] += pixel[0];
              rgb[1] += pixel[1];
              rgb[2] += pixel[2];
              c++;
            }
          }
        }
        clusters.push_back({
          i * POSWEIGTH,
          j * POSWEIGTH,
          rgb[0] / c, rgb[1] / c, rgb[2] / c
        });
        clusterColors.push_back({rgb[0] / c, rgb[1] / c, rgb[2] / c});
      }
    }
    int K = clusters.size();
    vector<vector<int>> clusterSums = vector<vector<int>> (K, vector<int> (D));
    vector<int> clusterCount = vector<int> (K, 0);
    vector<vector<vector<int>>> imgColors = vector<vector<vector<int>>> (img.W, vector<vector<int>> (img.H));
    while(true){
      cout << "iteration" << endl;
      for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
        pair<int, int> m = {INT_MAX, -1};
        vector<int> pixel = img.getPixel(i, j);
        vector<int> p = {
          i * POSWEIGTH, j * POSWEIGTH, pixel[0], pixel[1], pixel[2]
        };
        for(int k = 0; k < K; k++){
          int d = dist(p, clusters[k]);
          if(d < m.first){
            m = {d, k};
          }
        }
        //cout << m.first << " " << clusters[m.second][0] << " "<< clusters[m.second][1] << " " << p[0] << " " << p[1] << endl;
        clusterCount[m.second]++;
        for(int d = 0; d < D; d++){
          clusterSums[m.second][d] += p[d];
        }
        imgColors[i][j] = clusterColors[m.second];
      }
      vector<vector<int>> preClusters = clusters;
      for(int i = 0; i < K; i++){
        if(clusterCount[i] == 0) continue;
        for(int j = 0; j < D; j++){
          clusters[i][j] = clusterSums[i][j] / clusterCount[i];
        }
      }
      clusterCount = vector<int> (K);
      clusterSums = vector<vector<int>> (K, vector<int> (D));
      bool b = false;
      for(int i = 0; i < K; i++){
        if(dist(preClusters[i], clusters[i]) > thres){
          cout << dist(preClusters[i], clusters[i]) << endl;
          b = true;
          break;
        }
      }
      if(!b) break;
    }
    for(int i = 0; i < img.W; i++) for(int j = 0; j < img.H; j++){
      img.setPixel(i, j, imgColors[i][j]);
    }
  }

  struct ComplexNumber{
    double a;
    double b;
    ComplexNumber(){};
    ComplexNumber(double a, double b){
      this->a = a;
      this->b = b;
    }
  };

  ComplexNumber FT2D(vector<vector<int>>& v, int k, int l, int W, int H){
    ComplexNumber res (0.0, 0.0);
    for(int i = 0; i < W; i++){
      for(int j = 0; j < H; j++){
        double a = v[i][j]*cos(2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        double b = v[i][j]*sin(2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        res.a += a;
        res.b += b;
      }
    }
    return res;
  }

  void FFT_spectrum(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<int>> v (img.W, vector<int> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = brightness(img.getPixel(i, j));
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre (img.W, vector<ComplexNumber> (img.H));
    for(int k = 0; k < img.W; k++){
      for(int l = 0; l < img.H; l++){
        fre[k][l] = FT2D(v, k, l, img.W, img.H);
        int length = sqrt(fre[k][l].a * fre[k][l].a + fre[k][l].b * fre[k][l].b);
        cout << fre[k][l].a << " " << fre[k][l].b << endl;
        img.setPixel(k, l, {length, length, length});
      }
    }
  }
}
