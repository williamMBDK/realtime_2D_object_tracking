#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"IO.cpp"

namespace SEGMENT{
  // returns amount of pixels a given node contains and modifies pixelMap
  int traverseSegmentationTree(
    vector<int>& derived_nodes,
    vector<DATA::graph>& previousGraphs,
    vector<vector<int>>& pixelMap,
    IO::image& img,
    int node
  ){
    stack<pair<int, int>> s;
    int amountOfDerivedNodes = (int)derived_nodes.size(); // removes redundant consecutive .size() calls
    int amountOfGraphs = (int)previousGraphs.size();
    for(int i = 0; i < amountOfDerivedNodes; i++){
      s.push({amountOfGraphs - 1, derived_nodes[i]});
    }
    while(!s.empty()){
      pair<int, int> curr = s.top();
      s.pop();
      if(curr.first == 0){
        // this is a pixel
        int x = curr.second / img.H;
        int y = curr.second % img.H;
        pixelMap[x][y] = node;
      }else{
        int l = (int)previousGraphs[curr.first].derived_nodes[curr.second].size();
        for(int i = 0; i < l; i++){
          s.push({
            curr.first - 1,
            previousGraphs[curr.first].derived_nodes[curr.second][i]
          });
        }
      }
    }
  }

  // returns a graph_to_image_translator_object that describes the pixels that
  DATA::image_to_graph_translator_object getDerivedPixelsFromGraph(
    DATA::graph& currentGraph,
    vector<DATA::graph>& previousGraphs,
    IO::image& img
  ){
    DATA::image_to_graph_translator_object res (img.W, img.H);
    for(int i = 0; i < currentGraph.N; i++){
      traverseSegmentationTree(currentGraph.derived_nodes[i], previousGraphs, res.pixel_map, img, i);
    }
    return res;
  }
}
