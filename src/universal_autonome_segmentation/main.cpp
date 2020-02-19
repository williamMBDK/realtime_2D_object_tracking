#include<bits/stdc++.h>
#include"IO.cpp"
#include"data.cpp"
#include"seg1.cpp"
#include"seg2.cpp"
#include"merge.cpp"
#define current_time chrono::high_resolution_clock::now()
using namespace std;
double getDuration(auto start, auto stop){
  return ((double)(chrono::duration_cast<chrono::microseconds>(stop - start)).count())/1000.0;
}
void seg3(IO::image& img, auto& start){
  cout << "Approximate amount of segments: ";
  int SEGMENTS; cin >> SEGMENTS;
  start = current_time;
  DATA::pixel_graph g = DATA::imageToPixelGraph(img);
  cout << endl;
  while(g.N > SEGMENTS){
    cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
    SEG2::evaluateRegions2(g);
    g = MERGE::mergeRegions(g, true);
  }
  cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
  cout << endl;
  DATA::pixelGraphToIMG_averageColor(g, img);
  cout << "Amount of segments: "  << to_string(g.N) << endl;
}
void seg2(IO::image& img, auto& start){
  cout << "Approximate amount of segments: ";
  int SEGMENTS; cin >> SEGMENTS;
  start = current_time;
  DATA::pixel_graph g = DATA::imageToPixelGraph(img);
  cout << endl;
  int cnt = 1;
  while(g.N > SEGMENTS){
    cnt++;
    cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
    auto st = current_time;
    SEG2::evaluateRegions1(g);
    auto en = current_time;
    cout << "segmentation generation time: " << getDuration(st, en) << endl;
    st = current_time;
    g = MERGE::mergeRegions(g, true);
    en = current_time;
    cout << "segmentation merge time: " << getDuration(st, en) << endl;
  }
  cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
  cout << endl;
  DATA::pixelGraphToIMG_averageColor(g, img);
  cout << "Amount of segments in final segmentation: "  << to_string(g.N) << endl;
  cout << "Amount of segmentations: " << cnt << endl;
}
void seg1(IO::image& img, auto& start){
  cout << "Approximate amount of segments: ";
  int SEGMENTS; cin >> SEGMENTS;
  start = current_time;
  DATA::pixel_graph g = DATA::imageToPixelGraph(img);
  int initialAmountOfSegments = g.N;
  while(g.N > SEGMENTS){
    cout << "Found segmentation with " << to_string(g.N) << " segments" << endl;
    SEG1::evaluateRegions(g, initialAmountOfSegments);
    g = MERGE::mergeRegions(g, false);
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
    case 2: seg2(img, start); break;
    case 3: seg3(img, start); break;
    default:
      cout << "invalid segmentation type" << endl;
  }

  auto stop = chrono::high_resolution_clock::now();

  IO::writePPM(argv[3], img);

  double duration = getDuration(start, stop);
  cout << "Execution time: " << duration << " ms, excluding IO" << endl;
}
