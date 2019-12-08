#include<bits/stdc++.h>
#include"IO.cpp"
#include"data.cpp"
#include"flow.cpp"
#include"heu.cpp"
#include"preseg.cpp"
using namespace std;

int main(int argc, char const *argv[]) {
  //ios_base::sync_with_stdio(false);
  //cin.tie(NULL);
  if(argc < 3){
    cerr << "missing argument file" << endl;
    return 1;
  }
  IO::image img;
  IO::readPPM(argv[1], img);
  auto start = chrono::high_resolution_clock::now();

  /*// v1 tests
  IO::image img;
  IO::readPPM(argv[1], img);
  img.setPixel(2, 0, {100, 100, 100});
  IO::writePPM(argv[2], img);*/

  // v2 tests
  /*IO::image img;
  IO::readPPM(argv[1], img);
  pair<int, int> backgroundPixel, foregroundPixel;
  cout << "background pixel (0-indexed, x, y): ";
  cin >> backgroundPixel.first >> backgroundPixel.second;
  cout << "foreground pixel (0-indexed, x, y): ";
  cin >> foregroundPixel.first >> foregroundPixel.second;
  DATA::flow_graph flowGraph = DATA::imageToFlowGraph(img, backgroundPixel, foregroundPixel);
  //printFlowGraph(flowGraph);
  vector<pair<int, int>> min_cut = minCut(flowGraph.matrix, flowGraph.s, flowGraph.t);
  DATA::printMinCutImage(min_cut, img);
  DATA::applyMinCutOnImage(img, min_cut);
  IO::writePPM(argv[2], img);*/

  /*// heuristic: avgDiff test
  IO::image img;
  IO::readPPM(argv[1], img);
  cout << HEU::avgDiff(img) << endl;*/

  /*// presegmentation method 1.1 test
  IO::image img;
  IO::readPPM(argv[1], img);
  PRESEG::method1_1(img);
  IO::writePPM(argv[2], img);*/

  /*// presegmentation experimential test
  IO::image img;
  IO::readPPM(argv[1], img);
  PRESEG::experimential(img);
  IO::writePPM(argv[2], img);*/

  /*// presegmentation method 1.3 test
  IO::image img;
  IO::readPPM(argv[1], img);
  PRESEG::method1_3(img);
  IO::writePPM(argv[2], img);*/

  // k means (method 2)
  PRESEG::method2(img);

  // FT spectrum test
  //PRESEG::FT_spectrum(img);

  /*// methods 3.1 : FT
  IO::image img;
  IO::readPPM(argv[1], img);
  PRESEG::method3_1(img);
  IO::writePPM(argv[2], img);*/

  /*// methods 3.2 : FT
  IO::image img;
  IO::readPPM(argv[1], img);
  PRESEG::method3_2(img);
  IO::writePPM(argv[2], img);*/

  // FT spectrum dp test
  //PRESEG::FT_spectrum_dp(img);

  // FFT spectrum dp test
  //PRESEG::FFT_spectrum(img);

  auto stop = chrono::high_resolution_clock::now();
  IO::writePPM(argv[2], img);
  double duration = ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
  cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
