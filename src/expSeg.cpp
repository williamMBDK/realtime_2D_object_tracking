#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"

int brightness(vector<int> pixel){
  return (int)(0.2126*(double)pixel[0] + 0.7152*(double)pixel[1] + 0.0722*(double)pixel[2]);
}

void evaluateRegions(DATA::pixel_graph& g){
  int N_MAX = 10;
  double dT = 0.2;
  vector<vector<double>> dp (N_MAX, vector<double> (g.N, 0.0));
  for(int i = 0; i < g.N; i++){
    dp[0][i] = (double)pow(brightness(g.averagePixel[i]), 2);
  }
  for(int n = 1; n < N_MAX; n++){
    //cout << "iteration" << endl;
    vector<vector<double>> ks (4, vector<double> (g.N));
    for(int i = 0; i < g.N; i++){
      double m = DBL_MAX;
      double val = dp[n-1][i];
      for(int j = 0; j < g.adjacency_list[i].size(); j++){
        double nbVal = dp[n-1][g.adjacency_list[i][j]];
        if(abs(nbVal - val) < abs(m)){
          m = nbVal - val;
        }
      }
      ks[0][i] = m;
    }
    for(int i = 0; i < g.N; i++){
      double m = DBL_MAX;
      double val = dp[n-1][i] + ks[0][i] * dT / 2.0;
      for(int j = 0; j < g.adjacency_list[i].size(); j++){
        double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[0][g.adjacency_list[i][j]] * dT / 2.0;
        if(abs(nbVal - val) < abs(m)){
          m = nbVal - val;
        }
      }
      ks[1][i] = m;
    }
    for(int i = 0; i < g.N; i++){
      double m = DBL_MAX;
      double val = dp[n-1][i] + ks[1][i] * dT / 2.0;
      for(int j = 0; j < g.adjacency_list[i].size(); j++){
        double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[1][g.adjacency_list[i][j]] * dT / 2.0;
        if(abs(nbVal - val) < abs(m)){
          m = nbVal - val;
        }
      }
      ks[2][i] = m;
    }
    for(int i = 0; i < g.N; i++){
      double m = DBL_MAX;
      double val = dp[n-1][i] + ks[2][i] * dT;
      for(int j = 0; j < g.adjacency_list[i].size(); j++){
        double nbVal = dp[n-1][g.adjacency_list[i][j]] + ks[2][g.adjacency_list[i][j]] * dT;
        if(abs(nbVal - val) < abs(m)){
          m = nbVal - val;
        }
      }
      ks[3][i] = m;
    }
    for(int i = 0; i < g.N; i++){
      dp[n][i] = dp[n-1][i] + dT/6.0*(ks[0][i] + 2*ks[1][i] + 2*ks[2][i] + ks[3][i]);
    }
  }
  for(int i = 0; i < g.N; i++){
    dp[N_MAX - 1][i] = sqrt(dp[N_MAX - 1][i]);
    //dp[N_MAX - 1][i] = min(dp[N_MAX - 1][i], (double)g.MAX_RGB);
    //dp[N_MAX - 1][i] = max(dp[N_MAX - 1][i], 0.0);
    g.averagePixel[i] = {(int)dp[N_MAX - 1][i], (int)dp[N_MAX - 1][i], (int)dp[N_MAX - 1][i]};
    //if(g.N % 128 == 0) cout << "color" << endl;
  }
}

bool isInSideGrid(int x, int y, int W, int H){
  return x > -1 && x < W && y > -1 && y < H;
}

DATA::pixel_graph mergeRegions(DATA::pixel_graph& g){
  DATA::pixel_graph res;
  int component = -1;
  vector<int> components (g.N, -1);
  vector<vector<int>> pixelMap (g.W, vector<int> (g.W));
  for(int i = 0; i < g.N; i++){
    if(components[i] == -1){
      queue<int> q; q.push(i);
      component++;
      while(!q.empty()){
        int node = q.front(); q.pop();
        if(components[node] != -1) continue;
        components[node] = component;
        for(int j = 0; j < g.pixels[node].size(); j++){
          pixelMap[g.pixels[node][j].first][g.pixels[node][j].second] = component;
        }
        for(int j = 0; j < g.adjacency_list[node].size(); j++){
          if(brightness(g.averagePixel[g.adjacency_list[node][j]]) == brightness(g.averagePixel[node])){
            q.push(g.adjacency_list[node][j]);
          }
        }
      }
    }
  }
  res.N = component + 1;
  res.MAX_RGB = g.MAX_RGB; res.W = g.W; res.H = g.H;
  res.adjacency_list = vector<vector<int>> (res.N);
  res.pixels = vector<vector<pair<int, int>>> (res.N);
  res.averagePixel = vector<vector<int>> (res.N, {0, 0, 0});
  vector<pair<int, int>> dirs = {
    {-1, 0},
    {1, 0},
    {0, -1},
    {0, 1}
  };
  vector<unordered_set<int>> nbs (res.N);
  for(int i = 0; i < g.N; i++){
    for(int j = 0; j < g.pixels[i].size(); j++){
      pair<int, int> node = g.pixels[i][j];
      res.pixels[components[i]].push_back(node);
      res.averagePixel[components[i]][0] += g.averagePixel[i][0];
      res.averagePixel[components[i]][1] += g.averagePixel[i][1];
      res.averagePixel[components[i]][2] += g.averagePixel[i][2];
      for(int l = 0; l < 4; l++){
        int newX = node.first + dirs[l].first;
        int newY = node.second + dirs[l].second;
        if(isInSideGrid(newX, newY, g.W, g.H) && pixelMap[newX][newY] != components[i]){
          nbs[components[i]].insert(pixelMap[newX][newY]);
        }
      }
    }
  }
  for(int i = 0; i < res.N; i++){
    res.adjacency_list[i] = vector<int> (nbs[i].size());
    int idx = 0;
    for(int nb : nbs[i]){
      res.adjacency_list[i][idx++] = nb;
    }
    res.averagePixel[i][0] /= res.pixels[i].size();
    res.averagePixel[i][1] /= res.pixels[i].size();
    res.averagePixel[i][2] /= res.pixels[i].size();
  }
  return res;
}

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

void pixelGraphToImage_greyscale(DATA::pixel_graph& g, IO::image& img){
  for(int i = 0; i < g.N; i++){
    vector<int> color = g.averagePixel[i];
    color[0] = min(max(color[0], 0), img.MAX_RGB);
    color[1] = min(max(color[1], 0), img.MAX_RGB);
    color[2] = min(max(color[2], 0), img.MAX_RGB);
    for(int j = 0; j < g.pixels[i].size(); j++){
      img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
    }
  }
}

void pixelGraphToImage_color(DATA::pixel_graph& g, IO::image& img){
  for(int i = 0; i < g.N; i++){
    vector<int> color = {0, 0, 0};
    for(int j = 0; j < g.pixels[i].size(); j++){
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
    for(int j = 0; j < g.pixels[i].size(); j++){
      img.setPixel(g.pixels[i][j].first, g.pixels[i][j].second, color);
    }
  }
}

int main(int argc, char const *argv[]){
  if(argc < 3){
    cerr << "missing argument file" << endl;
    return 1;
  }
  IO::image img;
  IO::readPPM(argv[1], img);
  auto start = chrono::high_resolution_clock::now();

  DATA::pixel_graph g = imageToPixelGraph(img);
  int ITERATIONS = 23;
  for(int i = 0; i < ITERATIONS; i++){
    evaluateRegions(g);
    g = mergeRegions(g);
  }
  pixelGraphToImage_greyscale(g, img);

  auto stop = chrono::high_resolution_clock::now();
  IO::writePPM(argv[2], img);
  double duration = ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
  cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
