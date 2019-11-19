#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

int main(int argc, char const *argv[]) {
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);
  if(argc < 3){
    cerr << "missing argument file" << endl;
    return 1;
  }
  IO::image img;
  IO::readPPM(argv[1], img);
  img.setPixel(2, 0, {100, 100, 100});
  IO::writePPM(argv[2], img);
}
