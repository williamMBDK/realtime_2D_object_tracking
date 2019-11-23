// this file contains function to convert an image to a flow graph

#include"IO.cpp"

namespace DATA{
  bool isInSideGrid(int x, int y, int W, int H){
    return x > -1 && x < W && y > -1 && y < H;
  }

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
    double k = 10.0;
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
}
