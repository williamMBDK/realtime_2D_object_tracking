#pragma once
// this file contains function to convert an image to a flow graph
#include"IO.cpp"

namespace DATA{
  bool isInSideGrid(int x, int y, int W, int H){
    return x > -1 && x < W && y > -1 && y < H;
  }

  double brightness(vector<int> pixel){
    return 0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2];
  }

  // BEGIN NETWORK FLOW

  struct flow_graph{
    int s;
    int t;
    vector<vector<int>> matrix;
    flow_graph(int s, int t, vector<vector<int>>& matrix){
      this->matrix = matrix;
      this->s = s;
      this->t = t;
    }
  };

  int deltaPixel_v1(vector<int> p1, vector<int> p2){
    int sum = 0;
    for(int i = 0; i < 3; i++){
      sum += (p1[i] - p2[i])*(p1[i] - p2[i]);
    }
    int ma = 255 * 255 * 3 + 1;
    return ma - sum;
  }

  int deltaPixel_v2(vector<int> p1, vector<int> p2){
    double k = 6;
    double brightness1 = 0.2126*(double)p1[0] + 0.7152*(double)p1[1] + 0.0722*(double)p1[2];
    double brightness2 = 0.2126*(double)p2[0] + 0.7152*(double)p2[1] + 0.0722*(double)p2[2];
    double res = 100 * exp(-1 * pow(brightness1-brightness2, 2)/(2*k*k));
    cout << brightness1 << " " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
    //cout << (int) res + 1 << endl;
    return (int) res;
  }

  flow_graph imageToFlowGraph(IO::image &img,
  pair<int, int> component1Node, pair<int, int> component2Node){
    vector<vector<int>> matrix (img.W * img.H + 2, vector<int> (img.W * img.H + 2, -1));
    vector<pair<int, int>> dirs = {
      {1, 0},
      {-1, 0},
      {0, 1},
      {0, -1}
    };
    int ma = 0;
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        int idx = img.H * i + j;
        for(int k = 0; k < 4; k++){
          int newI = i + dirs[k].first;
          int newJ = j + dirs[k].second;
          if(isInSideGrid(newI, newJ, img.W, img.H)){
            int newIdx = img.H * newI + newJ; // mistake found here
            matrix[idx][newIdx] = deltaPixel_v2(
              img.getPixel(i, j),
              img.getPixel(newI, newJ)
            );
            ma = max(matrix[idx][newIdx], ma);
            //cout << i << " "<<j << " --> " << newI <<" "<< newJ << " of " << matrix[idx][newIdx] << endl;
            //cout << newI <<" "<< newJ<<" " << newIdx << "         ";
            //cout << matrix[idx][newIdx] << " ";
          }else{
            //cout << 0 << " ";
          }
        }
        //cout << endl;
      }
    }
    int s = img.W*img.H;
    int t = img.W*img.H + 1;
    int component1NodeIdx = component1Node.first * img.H + component1Node.second;
    int component2NodeIdx = component2Node.first * img.H + component2Node.second;
    matrix[s][component1NodeIdx] = INT_MAX;
    matrix[component2NodeIdx][t] = INT_MAX;
    return flow_graph(s, t, matrix);
  }

  void printFlowGraph(flow_graph g){
    for(int i = 0; i < g.matrix.size(); i++){
      for(int j = 0; j < g.matrix[i].size(); j++){
        cout << g.matrix[i][j] << " ";
      }
      cout << endl;
    }
  }

  void applyMinCutOnImage(IO::image &img, vector<pair<int, int>> &min_cut){
    for(int i = 0; i < min_cut.size(); i++){
      int x1 = min_cut[i].first / img.H;
      int y1 = min_cut[i].first % img.H;
      int x2 = min_cut[i].second / img.H;
      int y2 = min_cut[i].second % img.H;
      //cout << x1 << " " << y1 << endl;
      //cout << x2 << " " << y2 << endl;
      if(isInSideGrid(x1, y1, img.W, img.H)) img.setPixel(x1, y1, {0, img.MAX_RGB, 0}); // green
      if(isInSideGrid(x2, y2, img.W, img.H)) img.setPixel(x2, y2, {0, 0, img.MAX_RGB}); // blue
      // green is outer boundary and blue is inner
    }
  }

  void printMinCutImage(vector<pair<int, int>>& min_cut, IO::image &img){
    for(int i = 0; i < min_cut.size(); i++){
      cout << min_cut[i].first/img.H << " " << min_cut[i].first%img.H << "     " << min_cut[i].second/img.H << " " << min_cut[i].second%img.H << endl;
    }
  }

  //END NETWORK FLOW

  // BEGIN PIXEL GRAPHS
  struct pixel_graph{
    int N = 0; // number of nodes
    int W, H, MAX_RGB;
    vector<vector<int>> adjacency_list;
    vector<vector<pair<int, int>>> pixels;
    /*pixel_graph(int N){
      this->N = N;
      this->adjacency_list = vector<vector<int>> (N);
      this->pixels = vector<vector<int>> (N);
    }*/
  };

  pixel_graph getPixelGraph_binary(IO::image& img){
    int MAX_SIZE = 300;
    pixel_graph g;
    g.W = img.W;
    g.H = img.H;
    g.MAX_RGB = img.MAX_RGB;
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, 1},
      {0, -1}
    };
    vector<vector<int>> components = vector<vector<int>> (img.W, vector<int> (img.H));
    int component = 0;
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        if(components[i][j] == 0){ // unvisited
          component++;
          queue<pair<int, int>> q;
          q.push({i, j});
          int count = 0;
          while(count <= MAX_SIZE && !q.empty()){
            pair<int, int> node = q.front(); q.pop();
            if(components[node.first][node.second] != 0) continue;
            count++;
            components[node.first][node.second] = component;
            for(int l = 0; l < 4; l++){
              int newX = node.first + dirs[l].first;
              int newY = node.second + dirs[l].second;
              if(isInSideGrid(newX, newY, img.W, img.H) && brightness(img.getPixel(newX, newY)) == brightness(img.getPixel(node.first, node.second))){
                q.push({newX, newY});
              }
            }
          }
        }
      }
    }
    g.N = component;
    g.pixels = vector<vector<pair<int, int>>> (g.N);
    g.adjacency_list = vector<vector<int>> (g.N);
    for(int i = 0; i < img.W; i++){
      for(int j = 0; j < img.H; j++){
        g.pixels[components[i][j] - 1].push_back({i, j});
      }
    }
    for(int i = 0; i < g.N; i++){
      unordered_set<int> nbs; // neighbours
      for(int j = 0; j < g.pixels[i].size(); j++){
        pair<int, int> node = g.pixels[i][j];
        for(int l = 0; l < 4; l++){
          int newX = node.first + dirs[l].first;
          int newY = node.second + dirs[l].second;
          if(isInSideGrid(newX, newY, img.W, img.H) && components[newX][newY] - 1 != i){
            nbs.insert(components[newX][newY] - 1);
          }
        }
      }
      g.adjacency_list[i] = vector<int> (nbs.size());
      int idx = 0;
      for(int nb : nbs){
        g.adjacency_list[i][idx++] = nb;
      }
    }
    return g;
  }

  IO::image pixelGraphToIMG(pixel_graph& g){
    srand(time(NULL));
    IO::image img (g.W, g.H, g.MAX_RGB);
    for(int i = 0; i < g.N; i++){
      vector<int> color = {
        rand() % g.MAX_RGB,
        rand() % g.MAX_RGB,
        rand() % g.MAX_RGB
      };
      for(int j = 0; j < g.pixels[i].size(); j++){
        img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
      }
    }
    return img;
  }
}
