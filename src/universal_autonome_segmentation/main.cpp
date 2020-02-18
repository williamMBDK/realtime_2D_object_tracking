#include<bits/stdc++.h>
#include"IO.cpp"
#include"data.cpp"
#include"seg1.cpp"
#define current_time chrono::high_resolution_clock::now()
using namespace std;
void seg1(IO::image& img, auto& start){
  cout << "Approximate amount of segments: ";
  int SEGMENTS; cin >> SEGMENTS;
  start = current_time;
  DATA::pixel_graph g = DATA::imageToPixelGraph(img);
  int initialAmountOfSegments = g.N;
  while(g.N > SEGMENTS){
    cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
    SEG1::evaluateRegions(g, initialAmountOfSegments);
    g = SEG1::mergeRegions(g);
  }
  DATA::pixelGraphToIMG_averageColor(g, img);
  cout << "Amount of segments: "  << to_string(g.N) << endl;
}
int main(int argc, char const *argv[]){
  if(argc < 4){
    cerr << "missing arguments" << endl;
    return 1;
  }

  cout << "Initializing segmentation engine..." << endl;

  IO::image img;
  IO::readPPM(argv[2], img);

  int type = stoi(argv[1]);

  auto start = current_time;

  switch (type) {
    case 1: seg1(img, start); break;
    default:
      cout << "invalid segmentation type" << endl;
  }

  auto stop = chrono::high_resolution_clock::now();

  IO::writePPM(argv[3], img);

  double duration = ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
  cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
