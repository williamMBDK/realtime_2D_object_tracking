#pragma once
#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

namespace PRESEG{
  /*
  OVERVIEW AND NOTES:
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
    METHOD 1.6
      experimential
    METHOD 2:
      k-means clusering with k-means++
      dimensions:
        x,y,r,g,b
    METHOD 3:
      METHOD 3.1:
        band-pass-filter
      METHOD 3.2
        gaussian low pass filter
      METHOD 3.3
        fft band pass filter and contrast enhancement
      METHOD 3.4
        fft gaussian high pass filter and contract enhancement
    */

  // utility function
    // returns brightness of a pixel
    int brightness(vector<int> pixel){
      return (int)(0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2]);
    }

    // return boolean representating whether a point (x, y) is inside the grid (0, 0) to (W, H).
    bool isInGrid(int x, int y, int W, int H){
      return x > 0 && x < W-1 && y > 0 && y < H-1;
    }

    // return the L2 norm of to vectors
    int dist(vector<int> &v1, vector<int> &v2){
      int s = 0;
      for(int i = 0; i < v1.size(); i++){
        s += (v1[i] - v2[i])*(v1[i] - v2[i]);
      }
      return (int)sqrt(s);
    }

    // return mapping of v from range [0, r1] to range [0, r2]
    int mapRange(int v, int r1, int r2){
      return (int)((double)v/(double)r1*(double)r2);
    }

    // structure for storing a complex number, z = a + ib;
    struct ComplexNumber{
      double a;
      double b;
      ComplexNumber(){};
      ComplexNumber(double a, double b){
        this->a = a;
        this->b = b;
      }
    };

    // return the L2 norm of to vectors (double)
    double doubleDist(vector<double> v1, vector<double> v2){
      double s = 0;
      for(int i = 0; i < v1.size(); i++){
        s += (v1[i] - v2[i])*(v1[i] - v2[i]);
      }
      return sqrt(s);
    }

    // returns a ComplexNumber as the result of multiplying two ComplexNumber's a and b.
    ComplexNumber complexMult(ComplexNumber a, ComplexNumber b){
      ComplexNumber res (
        (a.a * b.a - a.b * b.b),
        (a.a * b.b + a.b * b.a)
      );
      return res;
    }

    // returns the addition of ComplexNumber a and ComplexNumber b as a new ComplexNumber
    ComplexNumber complexAdd(ComplexNumber a, ComplexNumber b){
      return ComplexNumber(
        a.a + b.a,
        a.b + b.b
      );
    }

  // modifies img with a pixel merging using dymanic programming and a discrete space
  // OPTIMIZATION : can be optimized by removing previous dp
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

  // modifies img with pixel merging using runge cutta
  // OPTIMIZATION: memory optimize dp
  void method1_3(IO::image& img){
    int N_MAX = 20;
    double dT = 0.1;
    vector<pair<int, int>> dirs = {
      {-1, 0}, {0, -1}, {1, 0}, {0, 1}
    };
    vector<vector<vector<double>>> dp (N_MAX, vector<vector<double>> (img.W, vector<double> (img.H, 0.0)));
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

  // modifies img with an experimential version of pixel merging
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

  // modifies img with k-means algorithm for finding superpixels
  void method2(IO::image& img){
    srand(time(NULL)); // seeding rng
    int thres = 10;
    int D = 5; // dimensions
    int MAX_VALUE = 300;
    int POSWEIGTH = 5;
    vector<vector<int>> clusters = vector<vector<int>> ();
    vector<vector<int>> clusterColors = vector<vector<int>> ();
    int size = 10;
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

  // returns the 1d discrete fourier transform of a specific frequency
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

  // modifies img to make it illustrate the fourier spectrum of the original img
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

  // returns the inverse fourier transform for a single value
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

  // modifies img with the naive (slow) ft and applys a band-pass-filter on the fourier transform matrix
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

  // modifies img with gaussian high pass filter applied to the fourier transform matrix (using naive method)
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

  // return a vector representing the 1D discrete fourier transform by using dynamic programming
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

  // returns a vector representing the inverse 1D discrete fourier transform by using dynamic programming
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
    return res;
  }

  // return a 2D vector representing the 2D discrete fourier transform by using dynamic programming (using  above method)
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

  // modifies img to become the fourier spectrum of the original img using dynamic programming
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

  // return a vector representing the 1D discrete fourier transform by using the fast fourier transform algorithm
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
    vector<ComplexNumber> res1 = fft1D(even, SIZE / 2);
    vector<ComplexNumber> res2 = fft1D(uneven, SIZE / 2);
    for(int i = 0; i < SIZE/2; i++){
      ComplexNumber c (
        cos(-2.0*M_PI*i/(double)SIZE),
        sin(-2.0*M_PI*i/(double)SIZE)
      );
      ComplexNumber t = complexMult(c, res2[i]);
      res[i] = complexAdd(res1[i], t);
      res[i + SIZE/2] = complexAdd(res1[i], complexMult(ComplexNumber(-1.0, 0.0), t));
    }
    return res;
  }

  // returns wheter x is a power of two
  // OBS: zero gives true
  bool isPowerOfTwo(int x){
    return (x ^ 0) == x - 1;
  }

  // returns the most significant bit in x
  int MSB(int x){
    for(int i = 31; i >= 0; i--){
      if(((x >> i) & 1) == 1) return i;
    }
    return -1;
  }

  // modifies W and H to represent a square of size MxM where M is a power of two
  void findProperDimensions(int& W, int& H){
    if(isPowerOfTwo(W) && isPowerOfTwo(H)) return;
    cout << MSB(W) << endl;
    W = (1 << (MSB(W) + 1));
    H = (1 << (MSB(W) + 1));
    if(W > H) H = W;
    else W = H;
  }

  // modifies v to be extended to fit dimensions W and H
  template <typename T>
  void applyProperDimensions(int W, int H, vector<vector<T>>& v){
    vector<vector<T>> res = vector<vector<T>>(W, vector<T> (H));
    for(int i = 0; i < v.size(); i++){
      for(int j = 0; j < v[i].size(); j++){
        res[i][j] = v[i][j];
      }
    }
    v = res;
  }

  // return a vector representing the 2D discrete fourier transform by using the fast fourier transform algorithm
  vector<vector<ComplexNumber>> fft2D(vector<vector<double>> v, int W, int H){
    int oriW = W;
    int oriH = H;
    findProperDimensions(W, H);
    applyProperDimensions<double>(W, H, v);
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
    vector<vector<ComplexNumber>> toReturn = vector<vector<ComplexNumber>>(oriW, vector<ComplexNumber> (oriH));
    for(int i = 0; i < oriW; i++){
      for(int j = 0; j < oriH; j++){
        toReturn[i][j] = res[i][j];
      }
    }
    return toReturn;
  }

  // return a vector representing the inverse 1D discrete fourier transform by using the fast fourier transform algorithm
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
    vector<ComplexNumber> res1 = fft1D_inverse(even, SIZE / 2);
    vector<ComplexNumber> res2 = fft1D_inverse(uneven, SIZE / 2);
    for(int i = 0; i < SIZE/2; i++){
      ComplexNumber c (
        cos(2.0*M_PI*i/(double)SIZE),
        sin(2.0*M_PI*i/(double)SIZE)
      );
      ComplexNumber t = complexMult(c, res2[i]);
      res[i] = complexAdd(res1[i], t);
      res[i + SIZE/2] = complexAdd(res1[i], complexMult(ComplexNumber(-1.0, 0.0), t));
    }
    return res;
  }

  // return a vector representing the inverse 2D discrete fourier transform by using the fast fourier transform algorithm
  vector<vector<double>> fft2D_inverse(vector<vector<ComplexNumber>>& v, int W, int H){
    int oriW = W;
    int oriH = H;
    findProperDimensions(W, H);
    applyProperDimensions<ComplexNumber>(W, H, v);
    vector<vector<ComplexNumber>> coloums (H, vector<ComplexNumber> (W)); // index[y][u]
    for(int u = 0; u < W; u++){
      vector<ComplexNumber> vZ = vector<ComplexNumber> (H);
      for(int y = 0; y < H; y++){
        vZ[y] = v[u][y];
      }
      vector<ComplexNumber> t = fft1D_inverse(vZ, H);
      for(int y = 0; y < H; y++){
        coloums[y][u] = complexMult(t[y], ComplexNumber(1.0/((double) H), 0.0));
      }
    }
    vector<vector<double>> res (W, vector<double> (H));
    for(int y = 0; y < H; y++){
      vector<ComplexNumber> t = fft1D_inverse(coloums[y], W);
      for(int x = 0; x < W; x++){
        res[x][y] = t[x].a / ((double) (W));
      }
    }
    vector<vector<double>> toReturn = vector<vector<double>>(oriW, vector<double> (oriH));
    for(int i = 0; i < oriW; i++){
      for(int j = 0; j < oriH; j++){
        toReturn[i][j] = res[i][j];
      }
    }
    return toReturn;
  }

  // modifies img to represents the original img's fourier spectrum using FFT
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

  // TEST FUNCTION
  void FFT_test(IO::image& img){
    // build 2d vector of brightness values
    vector<vector<double>> v (img.W, vector<double> (img.H));
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        v[i][j] = (double)brightness(img.getPixel(i, j))/(double)img.MAX_RGB;
      }
    }
    vector<vector<ComplexNumber>> fre = fft2D(v, img.W, img.H);
    vector<vector<double>> res = fft2D_inverse(fre, img.W, img.H);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int t = min(img.MAX_RGB, max((int)(v[i][j] * (double)img.MAX_RGB), 0));
        img.setPixel(i, j, {t, t, t});
      }
    }
  }

  // modifies img using band pass filter on the fourier transform matrix calculated using FFT
  void method3_3(IO::image& img){
    double FMAX = 0.7;
    double DMAX = min(img.W, img.H) * FMAX;
    double FMIN = 0.3;
    double DMIN = min(img.W, img.H) * FMIN;
    DMIN = 20;
    DMAX = 200;
    double contrastFactor = 5.0;
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
        double d = doubleDist({0.0, 0.0}, {img.W/2.0, img.W/2.0}) - doubleDist({(double)i, (double)j}, {img.W/2.0, img.W/2.0});
        if(d > DMAX || d < DMIN){
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

  // modifies img using gaussian high pass filter on the fourier transform matrix calculated using FFT
  void method3_4(IO::image& img){
    double contrastFactor = 50.0;
    double theta = 50.0;
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
        double d = doubleDist({0.0, 0.0}, {img.W/2.0, img.W/2.0}) - doubleDist({(double)i, (double)j}, {img.W/2.0, img.W/2.0});
        double f = 1.0 - exp(-1*(d*d)/(2.0*theta*theta));
        fre[i][j] = complexMult(fre[i][j], ComplexNumber(f, 0));
        fre[i][j] = complexMult(fre[i][j], ComplexNumber(contrastFactor, 0));
      }
    }
    vector<vector<double>> res = fft2D_inverse(fre, img.W, img.H);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int t = min(img.MAX_RGB, max((int)res[i][j], 0));
        if(t < img.MAX_RGB/2.0){
          img.setPixel(i, j, {0, 0, 0});
        }else{
          img.setPixel(i, j, {img.MAX_RGB, img.MAX_RGB, img.MAX_RGB});
        }
      }
    }
  }
}
