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

  // v1 tests
  //img.setPixel(2, 0, {100, 100, 100});

  // v2 tests
  /*pair<int, int> backgroundPixel, foregroundPixel;
  cout << "background pixel (0-indexed, x, y): ";
  cin >> backgroundPixel.first >> backgroundPixel.second;
  cout << "foreground pixel (0-indexed, x, y): ";
  cin >> foregroundPixel.first >> foregroundPixel.second;
  start = chrono::high_resolution_clock::now();
  DATA::flow_graph flowGraph = DATA::imageToFlowGraph(img, backgroundPixel, foregroundPixel);
  //printFlowGraph(flowGraph);
  vector<pair<int, int>> min_cut = minCut(flowGraph.matrix, flowGraph.s, flowGraph.t);
  //DATA::printMinCutImage(min_cut, img);
  DATA::applyMinCutOnImage(img, min_cut);*/

  /*// heuristic: avgDiff test
  cout << HEU::avgDiff(img) << endl;
  return 0;*/

  /*// presegmentation method 1.1 test
  PRESEG::method1_1(img);*/

  /*// presegmentation experimential test
  PRESEG::experimential(img);*/

  // presegmentation method 1.3 test
  //PRESEG::method1_3(img);

  // k means (method 2)
  //PRESEG::method2(img);

  // FT spectrum test
  //PRESEG::FT_spectrum(img);

  /*// methods 3.1 : FT
  PRESEG::method3_1(img);*/

  /*// methods 3.2 : FT
  PRESEG::method3_2(img);*/

  // FT spectrum dp test
  //PRESEG::FT_spectrum_dp(img);

  // FFT spectrum dp test
  //PRESEG::FFT_spectrum(img);

  // FFT general test
  //PRESEG::FFT_test(img);

  // method 3.3
  //PRESEG::method3_3(img);

  // method 1.3 + 3.3
  /*PRESEG::method1_3(img);
  PRESEG::method3_3(img);*/

  /*// method 3.4
  PRESEG::method3_4(img);*/

  /*// method 1.3 + 3.4
  PRESEG::method1_3(img);
  PRESEG::method3_4(img);*/

  // method 1.3 + 3.4 + visualisation
  /*IO::image temp = img;
  PRESEG::method1_3(temp);
  PRESEG::method3_4(temp);
  DATA::pixel_graph g = DATA::getPixelGraph(temp, false);
  DATA::applyHeuristicsToPixelGraph(g, img);
  DATA::pixelGraphToIMG(g, img);*/

  /*// method 2 + min cut
  pair<int, int> backgroundPixel, foregroundPixel;
  cout << "background pixel (0-indexed, x, y): ";
  cin >> backgroundPixel.first >> backgroundPixel.second;
  cout << "foreground pixel (0-indexed, x, y): ";
  cin >> foregroundPixel.first >> foregroundPixel.second;
  start = chrono::high_resolution_clock::now();
  IO::image temp = img;
  PRESEG::method2(temp);
  DATA::pixel_graph g = DATA::getPixelGraph(temp, true);
  DATA::applyHeuristicsToPixelGraph(g, img);
  //DATA::pixelGraphToIMG(g, img);
  DATA::flow_graph fg = DATA::getFlowGraphFromPixelGraph(g, backgroundPixel, foregroundPixel);
  vector<pair<int, int>> min_cut = minCut(fg.matrix, fg.s, fg.t);
  DATA::applyMinCutOnImageFromPixelGraph(img, min_cut, g, fg);
  //img = temp;*/

  // method 3.4 + 1-3 + min cut
  /*pair<int, int> backgroundPixel, foregroundPixel;
  cout << "background pixel (0-indexed, x, y): ";
  cin >> backgroundPixel.first >> backgroundPixel.second;
  cout << "foreground pixel (0-indexed, x, y): ";
  cin >> foregroundPixel.first >> foregroundPixel.second;
  start = chrono::high_resolution_clock::now();
  IO::image temp = img;
  PRESEG::method1_3(temp);
  PRESEG::method3_4(temp);
  DATA::pixel_graph g = DATA::getPixelGraph(temp, false);
  DATA::applyHeuristicsToPixelGraph(g, img);
  //DATA::pixelGraphToIMG(g, img);
  //DATA::pixelGraphToIMG_random(g, img);
  DATA::flow_graph fg = DATA::getFlowGraphFromPixelGraph(g, backgroundPixel, foregroundPixel);
  vector<pair<int, int>> min_cut = minCut(fg.matrix, fg.s, fg.t);
  DATA::applyMinCutOnImageFromPixelGraph(img, min_cut, g, fg);
  //img = temp;*/

  auto stop = chrono::high_resolution_clock::now();
  IO::writePPM(argv[2], img);
  double duration = ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
  cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
