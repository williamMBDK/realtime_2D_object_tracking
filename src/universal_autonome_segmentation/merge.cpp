#pragma once
#include<bits/stdc++.h>
#include"data.cpp"
#include"util.cpp"

namespace MERGE{
  DATA::graph mergeRegions(DATA::graph &g){
    int componentCounter = -1;
    vector<int> components (g.N, -1);
    vector<int> nodeCounts (g.N);
    for(int i = 0; i < g.N; i++){
      if(components[i] == -1){
        queue<int> q; q.push(i);
        componentCounter++;
        while(!q.empty()){
          int curr = q.front(); q.pop();
          if(components[curr] != -1) continue;
          nodeCounts[componentCounter]++;
          components[curr] = componentCounter;
          int l = (int)g.adjacency_list[curr].size();
          for(int j = 0; j < l; j++){
            if(
              //components[g.adjacency_list[curr][j]] == -1 &&
              DATA::squaredDifference(g.mean_vector[g.adjacency_list[curr][j]], g.mean_vector[curr]) < 2*2
            ){
               q.push(g.adjacency_list[curr][j]);
            }
          }
        }
      }
    }
    componentCounter++;
    DATA::graph res (componentCounter, g.MAX_RGB);
    vector<unordered_set<int>> nbs (componentCounter); // avoid hashsets - could be slow..
    for(int i = 0; i < g.N; i++){
      int l = g.adjacency_list[i].size();
      for(int j = 0; j < l; j++){
        int nb = g.adjacency_list[i][j];
        if(components[i] != components[nb]){
          nbs[components[i]].insert(components[nb]);
        }
      }
    }
    for(int i = 0; i < componentCounter; i++){
      res.adjacency_list[i] = vector<int> (nbs[i].size());
      res.derived_nodes[i] = vector<int> (nodeCounts[i]);
      int idx = 0;
      for(int nb : nbs[i]){
        res.adjacency_list[i][idx++] = nb;
      }
    }
    for(int i = 0; i < g.N; i++){
      res.derived_nodes[components[i]][--nodeCounts[components[i]]] = i;
      res.pixel_count[components[i]] += g.pixel_count[i];
      res.mean_vector[components[i]][0] += g.mean_vector[i][0] * g.pixel_count[i];
      res.mean_vector[components[i]][1] += g.mean_vector[i][1] * g.pixel_count[i];
      res.mean_vector[components[i]][2] += g.mean_vector[i][2] * g.pixel_count[i];
      /*vector<int> t = UTIL::mult_k(g.mean_vector[i], g.pixel_count[i]); // this could be optimized in terms of time but also it will overflow if to many pixels
      res.mean_vector[components[i]] = UTIL::add(
        res.mean_vector[components[i]],
        t
      );*/
    }
    for(int i = 0; i < componentCounter; i++){
      res.mean_vector[i][0] /= res.pixel_count[i];
      res.mean_vector[i][1] /= res.pixel_count[i];
      res.mean_vector[i][2] /= res.pixel_count[i];
      //res.mean_vector[i] = UTIL::div_k(res.mean_vector[i], res.pixel_count[i]);
    }
    return res;
  }

  /*// returns a pixelgraph where regions of same color in g has been merged
  DATA::pixel_graph mergeRegions(DATA::pixel_graph& g, bool withColor){
    DATA::pixel_graph res;
    int component = -1;
    vector<int> components (g.N, -1);
    vector<vector<int>> pixelMap (g.W, vector<int> (g.H));
    for(int i = 0; i < g.N; i++){
      if(components[i] == -1){
        queue<int> q; q.push(i);
        component++;
        while(!q.empty()){
          int node = q.front(); q.pop();
          if(components[node] != -1) continue;
          components[node] = component;
          for(int j = 0; j < (int)g.pixels[node].size(); j++){
            pixelMap[g.pixels[node][j].first][g.pixels[node][j].second] = component;
          }
          for(int j = 0; j < (int)g.adjacency_list[node].size(); j++){
            if(!withColor && DATA::brightness(g.averagePixel[g.adjacency_list[node][j]]) == DATA::brightness(g.averagePixel[node])){
              q.push(g.adjacency_list[node][j]);
            }else if(withColor && DATA::squaredDifference(g.averagePixel[g.adjacency_list[node][j]], g.averagePixel[node]) < 2*2){
              q.push(g.adjacency_list[node][j]);
            }
          }
        }
      }
    }
    res.N = component + 1;
    res.MAX_RGB = g.MAX_RGB; res.W = g.W; res.H = g.H;
    res.adjacency_list = vector<vector<int>> (res.N);
    res.pixels = vector<vector<pair<int, int>>> (res.N);
    res.averagePixel = vector<vector<int>> (res.N, {0, 0, 0});
    vector<pair<int, int>> dirs = {
      {-1, 0},
      {1, 0},
      {0, -1},
      {0, 1}
    };
    vector<unordered_set<int>> nbs (res.N);
    for(int i = 0; i < g.N; i++){
      for(int j = 0; j < (int)g.pixels[i].size(); j++){
        pair<int, int> node = g.pixels[i][j];
        res.pixels[components[i]].push_back(node);
        res.averagePixel[components[i]][0] += g.averagePixel[i][0];
        res.averagePixel[components[i]][1] += g.averagePixel[i][1];
        res.averagePixel[components[i]][2] += g.averagePixel[i][2];
        for(int l = 0; l < 4; l++){
          int newX = node.first + dirs[l].first;
          int newY = node.second + dirs[l].second;
          if(DATA::isInSideGrid(newX, newY, g.W, g.H) && pixelMap[newX][newY] != components[i]){
            nbs[components[i]].insert(pixelMap[newX][newY]);
          }
        }
      }
    }
    for(int i = 0; i < res.N; i++){
      res.adjacency_list[i] = vector<int> (nbs[i].size());
      int idx = 0;
      for(int nb : nbs[i]){
        res.adjacency_list[i][idx++] = nb;
      }
      res.averagePixel[i][0] /= res.pixels[i].size();
      res.averagePixel[i][1] /= res.pixels[i].size();
      res.averagePixel[i][2] /= res.pixels[i].size();
    }
    return res;
  }*/
}
