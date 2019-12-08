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
    METHOD 1.6
      experimential
    METHOD 2:
      k-means clusering with k-means++
      dimensions:
        x,y,r,g,b
    METHOD 3:
      3.1:
        band-pass-filter
      3.2
        gaussian low pass filter
      3.3
        fft band pass filter and contrast enhancement
    */

  // utility
    int brightness(vector<int> pixel){
      return (int)(0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2]);
    }

    bool isInGrid(int x, int y, int W, int H){
      return x > 0 && x < W-1 && y > 0 && y < H-1; // does not hit boundary
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

    struct ComplexNumber{
      double a;
      double b;
      ComplexNumber(){};
      ComplexNumber(double a, double b){
        this->a = a;
        this->b = b;
      }
    };

    double doubleDist(vector<double> v1, vector<double> v2){
      double s = 0;
      for(int i = 0; i < v1.size(); i++){
        s += (v1[i] - v2[i])*(v1[i] - v2[i]);
      }
      return sqrt(s);
    }

    ComplexNumber complexMult(ComplexNumber a, ComplexNumber b){
      ComplexNumber res (
        (a.a * b.a - a.b * b.b),
        (a.a * b.b + a.b * b.a)
      );
      return res;
    }

    ComplexNumber complexAdd(ComplexNumber a, ComplexNumber b){
      return ComplexNumber(
        a.a + b.a,
        a.b + b.b
      );
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

  // memory optimize dp
  void method1_3(IO::image& img){
    int N_MAX = 20;
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

  void method1_6(IO::image& img){
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
          //cout << dist(preClusters[i], clusters[i]) << endl;
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

  ComplexNumber ft1d_naive(vector<vector<double>>& v, int k, int l, int W, int H){
    ComplexNumber res (0.0, 0.0);
    for(int i = 0; i < W; i++){
      for(int j = 0; j < H; j++){
        double a = v[i][j]*cos(-2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        double b = v[i][j]*sin(-2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        res.a += a;
        res.b += b;
      }
    }
    return res;
  }

  void FT_spectrum(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j))/(double)img.MAX_RGB;
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre (img.W, vector<ComplexNumber> (img.H));
    for(int k = 0; k < img.W; k++){
      for(int l = 0; l < img.H; l++){
        fre[k][l] = ft1d_naive(v, k, l, img.W, img.H);
      }
    }
    for(int i = 0; i < img.W/2; i++){
      for(int j = 0; j < img.H/2; j++){
        int length = sqrt(fre[i][j].a * fre[i][j].a + fre[i][j].b * fre[i][j].b);
        img.setPixel(img.W/2 + i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 + i, img.H/2 - j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 - j, {length, length, length});
      }
    }
  }

  double ft1d_naive_inverse(vector<vector<ComplexNumber>> &fre, int i, int j, int W, int H){
    ComplexNumber res (0.0, 0.0);
    for(int k = 0; k < W; k++){
      for(int l = 0; l < H; l++){
        double a = cos(2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        double b = sin(2.0*M_PI*((double)k*(double)i/(double)W + (double)l*(double)j/(double)H));
        double c = fre[k][l].a;
        double d = fre[k][l].b;
        res.a += (a * c - b * d);
        res.b += (a * d + b * c);
      }
    }
    return res.a / (W * H);
  }

  void method3_1(IO::image& img){
    int FMAX = 2;
    int DMAX = min(img.W, img.H) / FMAX;
    int FMIN = 6;
    int DMIN = min(img.W, img.H) / FMIN;
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j));
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre (img.W, vector<ComplexNumber> (img.H));
    for(int k = 0; k < img.W; k++){
      for(int l = 0; l < img.H; l++){
        fre[k][l] = ft1d_naive(v, k, l, img.W, img.H);
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int d = sqrt(pow(i - img.W / 2, 2) + pow(j - img.H / 2, 2));
        if(d > DMAX || d < DMIN){
          fre[i][j] = ComplexNumber(0.0, 0.0);
        }
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int v = min(img.MAX_RGB, max((int)ft1d_naive_inverse(fre, i, j, img.W, img.H), 0));
        img.setPixel(i, j, {v, v, v});
      }
    }
  }

  void method3_2(IO::image& img){
    double theta = 100;
    double a = 0, b = 1;
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j));
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre (img.W, vector<ComplexNumber> (img.H));
    for(int k = 0; k < img.W; k++){
      for(int l = 0; l < img.H; l++){
        fre[k][l] = ft1d_naive(v, k, l, img.W, img.H);
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        double d = sqrt(i * i + j * j);
        double f = a + b * (1.0 - exp(-1*(d*d)/(2.0*theta*theta)));
        fre[i][j].a *= f;
        fre[i][j].b *= f;
      }
    }
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int v = min(img.MAX_RGB, max((int)ft1d_naive_inverse(fre, i, j, img.W, img.H), 0));
        img.setPixel(i, j, {v, v, v});
      }
    }
  }

  vector<ComplexNumber> ft1D_DP(vector<ComplexNumber>& v, int SIZE){
    vector<ComplexNumber> res (SIZE, ComplexNumber (0.0, 0.0));
    for(int k = 0; k < SIZE; k++){
      for(int i = 0; i < SIZE; i++){
        ComplexNumber t = complexMult(
          ComplexNumber(
            cos(-2.0*M_PI*(double)k*(double)i/(double)SIZE),
            sin(-2.0*M_PI*(double)k*(double)i/(double)SIZE)
          ),
          v.at(i)
        );
        res[k].a += t.a;
        res[k].b += t.b;
      }
    }
    return res;
  }

  vector<ComplexNumber> ft1D_DP_inverse(vector<ComplexNumber>& v, int SIZE){
    vector<ComplexNumber> res (SIZE, ComplexNumber (0.0, 0.0));
    for(int k = 0; k < SIZE; k++){
      for(int i = 0; i < SIZE; i++){
        ComplexNumber t = complexMult(
          ComplexNumber(
            cos(2.0*M_PI*(double)k*(double)i/(double)SIZE),
            sin(2.0*M_PI*(double)k*(double)i/(double)SIZE)
          ),
          v.at(i)
        );
        res[k].a += t.a;
        res[k].b += t.b;
      }
    }
    //for(int i = 0; i < SIZE; i++) res[i] = complexMult(res[i], ComplexNumber(1/(double)SIZE, 0.0));
    return res;
  }

  vector<vector<ComplexNumber>> ft2D_DP(vector<vector<double>>& v, int W, int H){
    vector<vector<ComplexNumber>> coloums (H, vector<ComplexNumber> (W)); // index[v][x]
    for(int i = 0; i < W; i++){
      vector<ComplexNumber> vZ = vector<ComplexNumber> (H);
      for(int j = 0; j < H; j++){
        vZ[j] = ComplexNumber (v[i][j], 0.0);
      }
      vector<ComplexNumber> t = ft1D_DP(vZ, H);
      for(int j = 0; j < H; j++){
        coloums[j][i] = t[j];
      }
    }
    vector<vector<ComplexNumber>> res (W, vector<ComplexNumber> (H));
    for(int i = 0; i < H; i++){
      vector<ComplexNumber> t = ft1D_DP(coloums[i], W);
      for(int j = 0; j < W; j++){
        res[j][i] = t[j];
      }
    }
    return res;
  }

  void FT_spectrum_dp(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j))/(double)img.MAX_RGB;
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre = ft2D_DP(v, img.W, img.H);
    for(int i = 0; i < img.W/2; i++){
      for(int j = 0; j < img.H/2; j++){
        int length = sqrt(fre[i][j].a * fre[i][j].a + fre[i][j].b * fre[i][j].b);
        img.setPixel(img.W/2 + i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 + i, img.H/2 - j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 - j, {length, length, length});
      }
    }
  }

  // compute res. vary y and v or vary x and u.
  vector<ComplexNumber> fft1D(vector<ComplexNumber> &v, int SIZE){
    if(SIZE == 1){
      return {v[0]};
    }
    vector<ComplexNumber> res (SIZE);
    vector<ComplexNumber> even (SIZE / 2);
    vector<ComplexNumber> uneven (SIZE / 2);
    for(int i = 0; i < SIZE; i+=2){
      even[i/2] = v[i];
      uneven[i/2] = v[i+1];
    }
    /*ComplexNumber W (
      1.0, 0.0
    );
    ComplexNumber WN (
      cos(-2.0*M_PI/(double)SIZE),
      sin(-2.0*M_PI/(double)SIZE)
    );*/
    vector<ComplexNumber> res1 = fft1D(even, SIZE / 2);
    vector<ComplexNumber> res2 = fft1D(uneven, SIZE / 2);
    for(int i = 0; i < SIZE/2; i++){
      ComplexNumber c (
        cos(-2.0*M_PI*i/(double)SIZE),
        sin(-2.0*M_PI*i/(double)SIZE)
      );
      ComplexNumber t = complexMult(c, res2[i]);
      //ComplexNumber t = complexMult(W, res2[i]);
      res[i] = complexAdd(res1[i], t);
      res[i + SIZE/2] = complexAdd(res1[i], complexMult(ComplexNumber(-1.0, 0.0), t));
      //W = complexMult(W, WN);
    }
    return res;
  }

  vector<vector<ComplexNumber>> fft2D(vector<vector<double>>& v, int W, int H){
    vector<vector<ComplexNumber>> coloums (H, vector<ComplexNumber> (W)); // index[v][x]
    for(int i = 0; i < W; i++){
      vector<ComplexNumber> vZ = vector<ComplexNumber> (H);
      for(int j = 0; j < H; j++){
        vZ[j] = ComplexNumber (v[i][j], 0.0);
      }
      vector<ComplexNumber> t = fft1D(vZ, H);
      for(int j = 0; j < H; j++){
        coloums[j][i] = t[j];
      }
    }
    vector<vector<ComplexNumber>> res (W, vector<ComplexNumber> (H));
    for(int i = 0; i < H; i++){
      vector<ComplexNumber> t = fft1D(coloums[i], W);
      for(int j = 0; j < W; j++){
        res[j][i] = t[j];
      }
    }
    return res;
  }

  vector<ComplexNumber> fft1D_inverse(vector<ComplexNumber> &v, int SIZE){
    if(SIZE == 1){
      return {v[0]};
    }
    vector<ComplexNumber> res (SIZE);
    vector<ComplexNumber> even (SIZE / 2);
    vector<ComplexNumber> uneven (SIZE / 2);
    for(int i = 0; i < SIZE; i+=2){
      even[i/2] = v[i];
      uneven[i/2] = v[i+1];
    }
    /*ComplexNumber W (
      1.0, 0.0
    );
    ComplexNumber WN (
      cos(2.0*M_PI/(double)SIZE),
      sin(2.0*M_PI/(double)SIZE)
    );*/
    vector<ComplexNumber> res1 = fft1D_inverse(even, SIZE / 2);
    vector<ComplexNumber> res2 = fft1D_inverse(uneven, SIZE / 2);
    for(int i = 0; i < SIZE/2; i++){
      ComplexNumber c (
        cos(2.0*M_PI*i/(double)SIZE),
        sin(2.0*M_PI*i/(double)SIZE)
      );
      ComplexNumber t = complexMult(c, res2[i]);
      //ComplexNumber t = complexMult(W, res2[i]);
      res[i] = complexAdd(res1[i], t);
      res[i + SIZE/2] = complexAdd(res1[i], complexMult(ComplexNumber(-1.0, 0.0), t));
      /*res[i] = complexMult(res[i], ComplexNumber(
        1.0/(double)SIZE, 0.0
      ));
      res[i + SIZE/2] = complexMult(res[i + SIZE/2], ComplexNumber(
        1.0/(double)SIZE, 0.0
      ));*/
      //W = complexMult(W, WN);
    }
    return res;
  }

  vector<vector<double>> fft2D_inverse(vector<vector<ComplexNumber>>& v, int W, int H){
    vector<vector<ComplexNumber>> coloums (H, vector<ComplexNumber> (W)); // index[y][u]
    for(int u = 0; u < W; u++){
      vector<ComplexNumber> vZ = vector<ComplexNumber> (H);
      for(int y = 0; y < H; y++){
        vZ[y] = v[u][y];
      }
      vector<ComplexNumber> t = fft1D_inverse(vZ, H);
      /*vector<ComplexNumber> t2 = ft1D_DP_inverse(vZ, H);
      for(int i = 0; i < H; i++){
        cout << t1[i].a << " " << t1[i].b << "    " << t2[i].a << " " << t2[i].b << endl;
      }*/
      for(int y = 0; y < H; y++){
        //coloums[y][u] = t[y];
        coloums[y][u] = complexMult(t[y], ComplexNumber(1.0/((double) H), 0.0));
      }
    }
    vector<vector<double>> res (W, vector<double> (H));
    for(int y = 0; y < H; y++){
      vector<ComplexNumber> t = fft1D_inverse(coloums[y], W);
      for(int x = 0; x < W; x++){
        res[x][y] = t[x].a / ((double) (W));
        //res[x][y] = t[x].a;
      }
    }
    return res;
  }

  void FFT_spectrum(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j))/(double)img.MAX_RGB;
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre = fft2D(v, img.W, img.H);
    for(int i = 0; i < img.W/2; i++){
      for(int j = 0; j < img.H/2; j++){
        int length = sqrt(fre[i][j].a * fre[i][j].a + fre[i][j].b * fre[i][j].b);
        img.setPixel(img.W/2 + i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 + j, {length, length, length});
        img.setPixel(img.W/2 + i, img.H/2 - j, {length, length, length});
        img.setPixel(img.W/2 - i, img.H/2 - j, {length, length, length});
      }
    }
  }

  void FFT_test(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j))/(double)img.MAX_RGB;
        //cout << brightness(img.getPixel(i, j)) << " " << img.getPixel(i, j)[0] << " " << img.getPixel(i, j)[1] << " " << img.getPixel(i, j)[2] << " " << endl;
      }
    }
    vector<vector<ComplexNumber>> fre = fft2D(v, img.W, img.H);
    vector<vector<double>> res = fft2D_inverse(fre, img.W, img.H);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int t = min(img.MAX_RGB, max((int)(v[i][j] * (double)img.MAX_RGB), 0));
        //cout << t << " " << res[i][j] << endl;
        img.setPixel(i, j, {t, t, t});
      }
    }
  }

  void method3_3(IO::image& img){
    double FMAX = 0.7;
    double DMAX = min(img.W, img.H) * FMAX;
    double FMIN = 0.3;
    double DMIN = min(img.W, img.H) * FMIN;
    DMIN = 40;
    DMAX = 60;
    double contrastFactor = 2.0;
    //cout << DMIN << " " << DMAX << endl;
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j));
      }
    }
    // now calculate the fourier transform for all frequencies
    vector<vector<ComplexNumber>> fre = fft2D(v, img.W, img.H);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        //double d = sqrt(pow(i, 2) + pow(j, 2));
        //double d = sqrt(pow(fre[i][j].a, 2.0) + pow(fre[i][j].b, 2.0));
        /*double d = min(
          min(i, j), min(img.W - i, img.H - j)
        );*/
        double d = doubleDist({0.0, 0.0}, {img.W/2.0, img.W/2.0}) - doubleDist({(double)i, (double)j}, {img.W/2.0, img.W/2.0});
        if(d > DMAX || d < DMIN){
          //cout << d << endl;
          fre[i][j] = ComplexNumber(0.0, 0.0);
        }else{
          fre[i][j] = complexMult(fre[i][j], ComplexNumber(contrastFactor, 0));
        }
      }
    }
    vector<vector<double>> res = fft2D_inverse(fre, img.W, img.H);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int t = min(img.MAX_RGB, max((int)res[i][j], 0));
        img.setPixel(i, j, {t, t, t});
      }
    }
  }
}
