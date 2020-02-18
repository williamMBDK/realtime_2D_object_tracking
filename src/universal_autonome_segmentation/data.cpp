#pragma once
#include"IO.cpp"

namespace DATA{
  // return boolean representating whether a point (x, y) is inside the grid (0, 0) to (W, H).
  bool isInSideGrid(int x, int y, int W, int H){
    return x > -1 && x < W && y > -1 && y < H;
  }

  // returns brightness of a pixel
  double brightness(vector<int> pixel){
    return 0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2];
  }

  // returns a heuristic for the difference between two pixels
  int deltaPixel_v1(vector<int> p1, vector<int> p2){
    int sum = 0;
    for(int i = 0; i < 3; i++){
      sum += (p1[i] - p2[i])*(p1[i] - p2[i]);
    }
    int ma = 255 * 255 * 3 + 1;
    return ma - sum;
  }

  // returns a heuristic for the difference between two pixels
  int deltaPixel_v2(vector<int> p1, vector<int> p2){
    double k = 40;
    double brightness1 = 0.2126*(double)p1[0] + 0.7152*(double)p1[1] + 0.0722*(double)p1[2];
    double brightness2 = 0.2126*(double)p2[0] + 0.7152*(double)p2[1] + 0.0722*(double)p2[2];
    double res = 100 * exp(-1 * pow(brightness1-brightness2, 2)/(2*k*k));
    return (int) res;
  }

  // structure to represent pixel graphs (superpixels)
  struct pixel_graph{
    int N = 0; // number of nodes
    int W, H, MAX_RGB;
    vector<vector<int>> adjacency_list;
    vector<vector<pair<int, int>>> pixels;
    vector<vector<int>> averagePixel;
  };

  //modifies g to contains neccesary heuristics about img
  void applyHeuristicsToPixelGraph(pixel_graph& g, IO::image& img){
    g.averagePixel = vector<vector<int>> (g.N, vector<int> (3));
    for(int i = 0; i < g.N; i++){
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        vector<int> pixel = img.getPixel(g.pixels[i][j].first, g.pixels[i][j].second);
        g.averagePixel[i][0] += pixel[0];
        g.averagePixel[i][1] += pixel[1];
        g.averagePixel[i][2] += pixel[2];
      }
      g.averagePixel[i][0] /= g.pixels[i].size();
      g.averagePixel[i][1] /= g.pixels[i].size();
      g.averagePixel[i][2] /= g.pixels[i].size();
    }
  }

  // checks if pixel1 equals pixel2
  bool compare(vector<int> pixel1, vector<int> pixel2){
    return pixel1[0] == pixel2[0] && pixel1[1] == pixel2[1] && pixel1[2] == pixel2[2];
  }

  // modifies img to contains the image represented by g
  void pixelGraphToIMG(pixel_graph& g, IO::image& img){
    for(int i = 0; i < g.N; i++){
      vector<int> color = g.averagePixel[i];
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
      }
    }
  }

  // modifies img to represent the pixelGraph but with regions having the average color of the original image
  void pixelGraphToIMG_averageColor(DATA::pixel_graph& g, IO::image& img){
    for(int i = 0; i < g.N; i++){
      vector<int> color = {0, 0, 0};
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        vector<int> pixel = img.getPixel(g.pixels[i][j].first, g.pixels[i][j].second);
        color[0] += pixel[0];
        color[1] += pixel[1];
        color[2] += pixel[2];
      }
      color[0] /= g.pixels[i].size();
      color[1] /= g.pixels[i].size();
      color[2] /= g.pixels[i].size();
      color[0] = min(max(color[0], 0), img.MAX_RGB);
      color[1] = min(max(color[1], 0), img.MAX_RGB);
      color[2] = min(max(color[2], 0), img.MAX_RGB);
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
      }
    }
  }

  // modifies img contain the regions from g but with random colors
  void pixelGraphToIMG_random(pixel_graph& g, IO::image& img){
    srand(time(NULL));
    for(int i = 0; i < g.N; i++){
      vector<int> color = {
        rand() % img.MAX_RGB,
        rand() % img.MAX_RGB,
        rand() % img.MAX_RGB
      };
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
      }
    }
  }

  // modifies img to represent the pixelGraph with grayscale
  void pixelGraphToIMG_greyscale(DATA::pixel_graph& g, IO::image& img){
    for(int i = 0; i < g.N; i++){
      vector<int> color = g.averagePixel[i];
      color[0] = min(max(color[0], 0), img.MAX_RGB);
      color[1] = min(max(color[1], 0), img.MAX_RGB);
      color[2] = min(max(color[2], 0), img.MAX_RGB);
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
      }
    }
  }

  // returns a pixelGraph that represents the image, img
  DATA::pixel_graph imageToPixelGraph(IO::image& img){
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, -1},
      {0, 1}
    };
    DATA::pixel_graph g;
    g.N = img.W * img.H;
    g.MAX_RGB = img.MAX_RGB; g.W = img.W; g.H = img.H;
    g.adjacency_list = vector<vector<int>> (g.N);
    g.pixels = vector<vector<pair<int, int>>> (g.N);
    g.averagePixel = vector<vector<int>> (g.N);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int node = i * img.H + j;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          int nb = newI * img.H + newJ;
          if(isInSideGrid(newI, newJ, img.W, img.H)){
            g.adjacency_list[node].push_back(nb);
          }
        }
        g.averagePixel[node] = img.getPixel(i, j);
        g.pixels[node].push_back({i, j});
      }
    }
    return g;
  }
}
