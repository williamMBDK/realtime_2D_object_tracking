#pragma once
#include<bits/stdc++.h>
using namespace std;
namespace IO{
  struct image{
    int H, W, MAX_RGB;
    vector<vector<vector<int>>> data;
    image(){};
    image(int W, int H, int MAX_RGB){
      this->data = vector<vector<vector<int>>>(W, vector<vector<int>> (H, vector<int> (3)));
      this->H = H;
      this->W = W;
      this->MAX_RGB = MAX_RGB;
    }
    void setPixel(int i, int j, vector<int> rgb){
      this->data.at(i).at(j) = rgb;
    }

    vector<int> getPixel(int i, int j){
      return this->data.at(i).at(j);
    }
  };

  void __test1(vector<vector<vector<int>>> &imageData){
    stringstream ss;
    cout << "test1" << endl;
    for(int i = 0; i < imageData.size(); i++){
      for(int j = 0; j < imageData[i].size(); j++){
        for(int k = 0; k < 3; k++){
          ss << imageData[i][j][k] << " ";
        }
        ss << "    ";
      }
      ss << endl;
    }
    cout << ss.str() << endl;
  }

  void readPPM(string file, image& img) {
    ifstream fileIn (file);
    if(fileIn.is_open()){
      string type; fileIn >> type;
      int W, H; fileIn >> W >> H;
      int MAX_RGB; fileIn >> MAX_RGB;
      cout << "Opened file of type " << type << " and dimensions " << W << " * " << H << " and max rgb value of " << MAX_RGB << endl;
      img = image(W, H, MAX_RGB);
      for(int i = 0; i < W; i++){
        for(int j = 0; j < H; j++){
          int r,g,b; fileIn >> r >> g >> b;
          img.setPixel(i, j, {r, g, b});
        }
      }
      cout << "Succesfully read image data from " << file << endl;
      fileIn.close();
      //__test1(img.data);
    }else{
      cerr << "could not open file" << endl;
      return;
    }
  }

  void writePPM(string file, image& img){
    ofstream fileOut (file);
    if(fileOut.is_open()){
      fileOut << "P3" << endl;
      fileOut << img.W << " " << img.H << endl;
      fileOut << img.MAX_RGB << endl;
      for(int i = 0; i < img.W; i++){
        for(int j = 0; j < img.H; j++){
          vector<int> pixel = img.getPixel(i, j);
          fileOut << pixel.at(0) << " " << pixel.at(1) << " " << pixel.at(2) << endl;
        }
      }
      //__test1(img.data);
      fileOut.close();
      cout << "Succesfully wrote image data to file " << file << endl;
    }else{
      cerr << "could not write to file" << endl;
    }
  }
}
