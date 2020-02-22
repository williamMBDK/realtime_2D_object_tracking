#include<bits/stdc++.h>
#include"IO.cpp"
#include"data.cpp"
#include"seg1.cpp"
#include"seg2.cpp"
#include"merge.cpp"
#include"segmentation.cpp"
#include"util.cpp"
#define current_time chrono::high_resolution_clock::now()
using namespace std;

double getDuration(auto start, auto stop){
  return ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
}

void seg3(IO::image& img, auto& start){

  bool withPosition = true;
  bool randomized = false;

  start = current_time;

  DATA::graph currentGraph = DATA::imageToGraph(img, withPosition);
  vector<DATA::graph> graphs;

  while(currentGraph.N != 1){
    SEG2::evaluateRegions1(currentGraph);
    graphs.push_back(currentGraph);
    currentGraph = MERGE::mergeRegions(currentGraph);
  }
  cout << (int)graphs.size() + 1 << " " << img.W << " " << img.H << endl;
  DATA::graph& resGraph = currentGraph;
  DATA::image_to_graph_translator_object itg = SEGMENT::getDerivedPixelsFromGraph(resGraph, graphs, img);
  if(randomized) DATA::pixelGraphToIMG_random_image(resGraph, itg, img);
  else DATA::pixelGraphToIMG_averageColor_image(resGraph, itg, img);
}

void seg2(IO::image& img, auto& start){

  cout << "Approximate amount of segments: ";
  int SEGMENTS; cin >> SEGMENTS;
  cout << "With position as parameters (1 or 0): ";
  bool withPosition; cin >> withPosition;
  cout << "Randomized color in output (1 or 0): ";
  bool randomized; cin >> randomized;
  cout << endl;

  start = current_time;

  DATA::graph currentGraph = DATA::imageToGraph(img, withPosition);
  vector<DATA::graph> graphs;

  while(currentGraph.N > SEGMENTS){
    cout << "Found segmentation with " << to_string(currentGraph.N) << " segments" << endl;
    auto st = current_time;
    SEG2::evaluateRegions1(currentGraph);
    auto en = current_time;
    graphs.push_back(currentGraph);
    cout << "segmentation generation time: " << getDuration(st, en) << endl;
    st = current_time;
    currentGraph = MERGE::mergeRegions(currentGraph);
    en = current_time;
    cout << "segmentation merge time: " << getDuration(st, en) << endl;
  }
  DATA::graph& resGraph = currentGraph;
  cout << "Found segmentation with " << to_string(resGraph.N) << " segments" << endl << endl;

  cout << "Amount of segments in final segmentation: "  << to_string(resGraph.N) << endl;
  cout << "Amount of segmentations: " << (int)graphs.size() + 1 << endl;

  auto st = current_time;
  DATA::image_to_graph_translator_object itg = SEGMENT::getDerivedPixelsFromGraph(resGraph, graphs, img);
  auto en = current_time;
  cout << "get derived pixels time: " << getDuration(st, en) << endl;
  if(randomized) DATA::pixelGraphToIMG_random_image(resGraph, itg, img);
  else DATA::pixelGraphToIMG_averageColor_image(resGraph, itg, img);
}

int main(int argc, char const *argv[]){
  if(argc < 4){
    cerr << "missing arguments" << endl;
    return 1;
  }

  if(argc == 4) cout << "Initializing segmentation engine..." << endl;

  IO::image img;
  IO::readPPM(argv[2], img, argc == 4);

  int type = stoi(argv[1]);

  auto start = current_time;

  switch (type) {
    case 2: seg2(img, start); break;
    case 3: seg3(img, start); return 0;
    default:
      cout << "invalid segmentation type" << endl;
  }

  auto stop = chrono::high_resolution_clock::now();

  IO::writePPM(argv[3], img, argc == 4);

  double duration = getDuration(start, stop);
  if(argc == 4) cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
