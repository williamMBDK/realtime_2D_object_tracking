#include<bits/stdc++.h>
using namespace std;
int main(){
  ifstream in;
  in.open("out/universal_autonome_segmentation/results.txt");
  //in.open("out/universal_autonome_segmentation/test.txt");
  int N; in >> N;
  double sumDiff = 0.0;
  for(int i = 0; i < N; i++){
    double s, w, h; in >> s >> w >> h; s--;
    double real = log2(w*h);
    sumDiff += (real - s);
  }
  cout << sumDiff / N << endl;
}
