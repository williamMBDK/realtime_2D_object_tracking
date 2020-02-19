#pragma once
#include"IO.cpp"
#include"util.cpp"

namespace DATA{

  struct image_to_graph_translator_object{
    int W, H;
    vector<vector<int>> pixel_map;
    image_to_graph_translator_object(int W, int H){
      this->W = W;
      this->H = H;
      this->pixel_map = vector<vector<int>> (W, vector<int> (H));
    }
  };

  struct graph{
      int N;
      int MAX_RGB;
      vector<vector<int>> derived_nodes;
      vector<vector<int>> adjacency_list;
      vector<vector<int>> mean_vector;
      vector<int> pixel_count;
      graph(int N, int MAX_RGB){
        this->N = N;
        this->MAX_RGB = MAX_RGB;
        derived_nodes = vector<vector<int>>(N);
        adjacency_list = vector<vector<int>>(N);
        mean_vector = vector<vector<int>>(N, vector<int> (3));
        pixel_count = vector<int> (N);
      };
      graph(){};
  };

  // difference between two vectors
  int squaredDifference(vector<int>& a, vector<int> &b){
    int s = 0;
    int l = (int)a.size();
    for(int i = 0; i < l; i++) s += (a[i] - b[i])*(a[i] - b[i]);
    return s;
  }

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

  // checks if pixel1 equals pixel2
  bool compare(vector<int> pixel1, vector<int> pixel2){
    return pixel1[0] == pixel2[0] && pixel1[1] == pixel2[1] && pixel1[2] == pixel2[2];
  }

  // modifies img to contains the image represented by g
  void pixelGraphToIMG_averageColor(
    graph& g,
    image_to_graph_translator_object& itg,
    IO::image& img
  ){
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        img.setPixel(i, j, g.mean_vector[itg.pixel_map[i][j]]);
      }
    }
  }

  // modifies img to represent the pixelGraph but with regions having the average color of the original image
  void pixelGraphToIMG_averageColor_image(
    graph& g,
    image_to_graph_translator_object& itg,
    IO::image& img
  ){
    vector<vector<int>> average_colors (g.N, vector<int> (3));
    for(int i = 0; i < itg.W; i++){
      for(int j = 0; j < itg.H; j++){
        vector<int> pixel = img.getPixel(i, j);
        average_colors[itg.pixel_map[i][j]][0] += pixel[0];
        average_colors[itg.pixel_map[i][j]][1] += pixel[1];
        average_colors[itg.pixel_map[i][j]][2] += pixel[2];
      }
    }
    for(int i = 0; i < g.N; i++){
      average_colors[i] = UTIL::div_k(average_colors[i], g.pixel_count[i]);
    }
    for(int i = 0; i < itg.W; i++){
      for(int j = 0; j < itg.H; j++){
        img.setPixel(i, j, average_colors[itg.pixel_map[i][j]]);
      }
    }
  }

  /*// modifies img contain the regions from g but with random colors
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
  }*/

  // returns a graph that represents the image, img
  graph imageToGraph(IO::image& img){
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, -1},
      {0, 1}
    };
    graph g (img.W * img.H, img.MAX_RGB);
    for(int i = 0; i < g.N; i++) g.pixel_count[i] = 1;
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int node = i * img.H + j;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInSideGrid(newI, newJ, img.W, img.H)){
            int nb = newI * img.H + newJ;
            g.adjacency_list[node].push_back(nb); // SPEEDUP: not ideal
          }
        }
        g.mean_vector[node] = img.getPixel(i, j);
      }
    }
    return g;
  }
}
