#include<bits/stdc++.h>
#include"IO.cpp"
using namespace std;

// returns an image with random brightness values
IO::image random(int W, int H, int MAX_RGB){
  srand(time(NULL)); // seconds from 1970
  IO::image img (W, H, MAX_RGB);
  for(int i = 0; i < W; i++){
    for(int j = 0; j < H; j++){
      int v = rand() % (MAX_RGB+1);
      img.setPixel(i, j,
        {v, v, v}
      );
    }
  }
  return img;
}

// return boolean representating whether a point (x, y) is inside the grid (0, 0) to (W, H).
bool check(int x, int y, int W, int H){
  return x > 0 && x < W-1 && y > 0 && y < H-1; // does not hit boundary
}

// returns and image with an object that has grown from starting point randomly
IO::image uniformly_random_object(int W, int H, int MAX_RGB){
  cout << "size (1-100): ";
  double size; cin >> size;
  srand(time(NULL)); // seconds from 1970
  priority_queue<pair<int, int>> q;
  q.push({0, (rand() % W) * H + rand() % H});
  IO::image img (W, H, MAX_RGB);
  for(int i = 0; i < W; i++) for(int j = 0; j < H; j++){
    img.setPixel(i, j, {255, 255, 255});
  }
  vector<pair<int, int>> dirs = {
    {-1, 0},
    {1, 0},
    {0, -1},
    {0, 1}
  };
  vector<bool> visited (W*H, false);
  for(int i = 0; i < (size/100.0)*W*H; i++){
    pair<int, int> curr = q.top(); q.pop();
    if(visited[curr.second]){
      i--;
      continue;
    }
    visited[curr.second] = true;
    int x = curr.second / H;
    int y = curr.second % H;
    img.setPixel(x, y, {0, 0, 0});
    for(int j = 0; j < 4; j++){
      int newX = x + dirs[j].first, newY = y + dirs[j].second;
      if(check(newX, newY, W, H)) q.push({rand(), newX * H + newY});
    }
  }
  return img;
}

// return distance between two points
int dist(pair<int, int> p1, pair<int, int> p2){
  return abs(p1.first - p2.first) + abs(p1.second - p2.second);
}

// returns an image with an object with the specified properties
IO::image object(int W, int H, int MAX_RGB){
  cout << "noise (0-100): ";
  double noise; cin >> noise;
  cout << "size (1-100): ";
  double size; cin >> size;
  srand(time(NULL)); // seconds from 1970
  priority_queue<pair<int, int>> q;
  pair<int, int> start = {W/2, H/2}; // {rand() % W, rand() % H}; // middle
  q.push({0, start.first * H + start.second});
  IO::image img (W, H, MAX_RGB);
  for(int i = 0; i < W; i++) for(int j = 0; j < H; j++){
    img.setPixel(i, j, {255, 255, 255});
  }
  vector<pair<int, int>> dirs = {
    {-1, 0},
    {1, 0},
    {0, -1},
    {0, 1}
  };
  vector<bool> visited (W*H, false);
  for(int i = 0; i < (size/100.0)*W*H; i++){
    pair<int, int> curr = q.top(); q.pop();
    if(visited[curr.second]){
      i--;
      continue;
    }
    visited[curr.second] = true;
    int x = curr.second / H;
    int y = curr.second % H;
    img.setPixel(x, y, {0, 0, 0});
    for(int j = 0; j < 4; j++){
      int newX = x + dirs[j].first, newY = y + dirs[j].second;
      if(check(newX, newY, W, H)){
        int spread = noise == 0.0 ? 0 : rand() % (int)(((double)W+(double)H)*noise/100.0); // when noise is 100 spread and distance are equally weigthed
        int distance = W+H - dist({newX, newY}, start); // we want to draw close to the start point
        int distanceToEdge = min(min(newX, newY), min(W-newX, H-newY)); // we dont want to draw close to the edges of the image
        q.push({spread + distance + distanceToEdge, newX * H + newY});
      }
    }
  }
  return img;
}

// main function being run at execution
int main(){
  cout << "Type (random/randomobject/object): ";
  string type; cin >> type;
  cout << "width : ";
  int W; cin >> W;
  cout << "height : ";
  int H; cin >> H;
  cout << "max rgb value : ";
  int MAX_RGB; cin >> MAX_RGB;
  cout << "output file: ";
  string out; cin >> out;
  IO::image img;
  if(type == "random") img = random(W, H, MAX_RGB);
  else if(type == "randomobject") img = uniformly_random_object(W, H, MAX_RGB);
  else if(type == "object") img = object(W, H, MAX_RGB);
  else{
    cerr << "not a valid type" << endl;
    return 1;
  }
  cout << "Succesfully generated data. Press to continue.";
  string t; getline(cin, t); getline(cin, t);
  IO::writePPM(out, img);
}
