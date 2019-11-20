#include<bits/stdc++.h>
#include"IO.cpp"
#include"data.cpp"
#include"flow.cpp"
using namespace std;

int main(int argc, char const *argv[]) {
  //ios_base::sync_with_stdio(false);
  //cin.tie(NULL);
  if(argc < 3){
    cerr << "missing argument file" << endl;
    return 1;
  }
  // v1 tests
  /*IO::image img;
  IO::readPPM(argv[1], img);
  img.setPixel(2, 0, {100, 100, 100});
  IO::writePPM(argv[2], img);*/

  // v2 tests
  IO::image img;
  IO::readPPM(argv[1], img);
  pair<int, int> backgroundPixel, foregroundPixel;
  cout << "background pixel (0-indexed, x, y): ";
  cin >> backgroundPixel.first >> backgroundPixel.second;
  cout << "foreground pixel (0-indexed, x, y): ";
  cin >> foregroundPixel.first >> foregroundPixel.second;
  DATA::flow_graph flowGraph = DATA::imageToFlowGraph(img, backgroundPixel, foregroundPixel);
  //printFlowGraph(flowGraph);
  vector<pair<int, int>> min_cut = minCut(flowGraph.matrix, flowGraph.s, flowGraph.t);
  //printMinCut(min_cut);
  DATA::applyMinCutOnImage(img, min_cut);
  IO::writePPM(argv[2], img);
}
