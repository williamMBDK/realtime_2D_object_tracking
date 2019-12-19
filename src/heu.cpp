#pragma once
#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

namespace HEU{
  // return boolean representating whether a point (x, y) is inside the grid (0, 0) to (W, H).
  bool isInSideGrid(int x, int y, int W, int H){
    return x > -1 && x < W && y > -1 && y < H;
  }

  // returns the average brightness difference of pairs of pixels that neigbours in img
  int avgDiff(IO::image img){
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, -1},
      {0, 1}
    };
    double sum = 0;
    int c = 0;
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInSideGrid(newI, newJ, img.W, img.H)){
            c++;
            vector<int> p1 = img.getPixel(i, j);
            vector<int> p2 = img.getPixel(newI, newJ);
            double brightness1 = 0.2126*(double)p1[0] + 0.7152*(double)p1[1] + 0.0722*(double)p1[2];
            double brightness2 = 0.2126*(double)p2[0] + 0.7152*(double)p2[1] + 0.0722*(double)p2[2];
            sum += abs(brightness1 - brightness2);
          }
        }
      }
    }
    return sum / (double)c;
  }
}
