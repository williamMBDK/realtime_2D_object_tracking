#include<bits/stdc++.h>
using namespace std;

#define debug(x) cout << #x << " : " << x << endl

void adjustResidualGraph(vector<int> &path, vector<vector<int>> &matrix, int flow){
  for(int i = path.size() - 1; i > 0; i--){
    matrix[path[i]][path[i-1]] -= flow;
    matrix[path[i-1]][path[i]] += flow;
  }
}

int getMinFlowOnPath(vector<int> &path, vector<vector<int>> &matrix){
  int maxEdgeCapacity = INT_MAX;
  for(int i = path.size() - 1; i > 0; i--){
    int edgeCapacity = matrix[path[i]][path[i-1]];
    maxEdgeCapacity = min(maxEdgeCapacity, edgeCapacity);
  }
  return maxEdgeCapacity;
}

void dfs_parent(vector<vector<int>> &matrix, vector<int>& parents, int node){
  for(int i = 0; i < matrix.size(); i++){
    //cout << node << " " << i << endl;
    //cout << node / 15 << " "<< node % 15 << " --> " << i / 15 << " "<< i % 15 << endl;
    if(matrix[node][i] != 0 && parents[i] == -1){
      //cout << node / 15 << " "<< node % 15 << " --> " << i / 15 << " "<< i % 15 << endl;
      parents[i] = node;
      dfs_parent(matrix, parents, i);
    }
  }
}

/*
  optimering: bruge bfs
  optimering: stoppe ligeså snart t er noget i bfs metoden
*/
vector<int> getPath(vector<vector<int>> &matrix, int s, int t){
  vector<int> parents (matrix.size(), -1);
  dfs_parent(matrix, parents, s);
  if(parents[t] == -1) return vector<int> (0);
  vector<int> path = {t};
  int node = t;
  while(node != s){
    node = parents[node];
    path.push_back(node);
  }
  return path;
}

void dfs_markVisited(vector<vector<int>> &matrix, int node, vector<bool> &visited){
  if(visited[node]) return;
  visited[node] = true;
  for(int i = 0; i < matrix.size(); i++){
    if(matrix[node][i] != 0){
      dfs_markVisited(matrix, i, visited);
    }
  }
}

int maxFlow(vector<vector<int>> &matrix, int s, int t){
  int res = 0;
  int flow = -1;
  while(flow != 0){
    vector<int> path = getPath(matrix, s, t);
    //debug(path.size());
    if(path.size() != 0){
      flow = getMinFlowOnPath(path, matrix);
      adjustResidualGraph(path, matrix, flow);
    }else{
      flow = 0;
    }
    /*debug(path.size());
    for(int i = 0; i < path.size(); i++) cout << (char)(path[i] + (int)'A') << " ";
    cout << endl;
    debug(flow);*/
    res += flow;
  }
  return res;
}

vector<pair<int, int>> minCut(vector<vector<int>> &matrix, int s, int t){
  vector<vector<int>> original = matrix; // copy
  maxFlow(matrix, s, t);
  vector<bool> visited (matrix.size());
  dfs_markVisited(matrix, s, visited);
  vector<pair<int, int>> cut;
  for(int i = 0; i < original.size(); i++){
    if(visited[i]){
      for(int j = 0; j < original[i].size(); j++){
        if(!visited[j] && original[i][j] != 0){
          cut.push_back({i, j});
        }
      }
    }
  }
  return cut;
}

void validateMaxFlowTestData(string infileString, string answerfileString, int test){
  cout << "testing " << infileString << endl;
  ifstream inFile (infileString);
  int M; inFile >> M; // connections
  char sC, tC; inFile >> sC >> tC;
  int s = (int)sC - (int)'A';
  int t = (int)tC - (int)'A';
  vector<vector<int>> matrix (26, vector<int> (26, 0));
  for(int i = 0; i < M; i++){
    char aC, bC; inFile >> aC >> bC;
    int flow; inFile >> flow;
    int a = (int)aC - (int)'A';
    int b = (int)bC - (int)'A';
    matrix[a][b] = flow;
    matrix[b][a] = flow;
  }
  int ans = maxFlow(matrix, s, t);
  ifstream answerFile (answerfileString);
  int corr; answerFile >> corr;
  if(ans == corr) cout << "\033[1;32mpassed test \033[0m" << test << " with answer = " << ans << " and correct = " << corr << endl;
  else cout << "\033[1;31mfailed test \033[0m" << test << " with answer = " << ans  << " and correct = " << corr << endl;
}

void validateAllMaxFlowTestData(){
  cout << "MAX FLOW VALIDATION" << endl;
  for(int i = 1; i <= 10; i++){
    validateMaxFlowTestData("./tests/maxflowtest/in" + to_string(i) + ".txt", "./tests/maxflowtest/corr" + to_string(i) + ".txt", i);
  }
}

int countComponents(vector<vector<int>>& matrix){
  int count = 0;
  vector<bool> v (matrix.size());
  for(int i = 0; i < matrix.size(); i++){
    if(!v[i]){
      count++;
      dfs_markVisited(matrix, i, v);
    }
  }
  return count;
}

bool isBinaryPartition(vector<vector<int>>& matrix){
  vector<vector<int>> copy = matrix;
  for(int i = 0; i < copy.size(); i++){
    for(int j = 0; j < copy.size(); j++){
      if(copy[i][j] != 0 && copy[j][i] == 0){
        copy[j][i] = copy[i][j];
      }
    }
  }
  //cout << countComponents(matrix) << " " << countComponents(copy) << endl;
  return countComponents(copy) == 2;
}

void printMinCut(vector<pair<int, int>>& min_cut){
  for(int i = 0; i < min_cut.size(); i++){
    cout << (char)(min_cut[i].first + (int)'A') << " " << (char)(min_cut[i].second + (int)'A') << endl;
  }
}

// assumes max flow algorithm is correct
void validateMinCut(string infileString, string answerfileString, int test){
  cout << "testing " << infileString << endl;
  ifstream inFile (infileString);
  int M; inFile >> M; // connections
  char sC, tC; inFile >> sC >> tC;
  int s = (int)sC - (int)'A';
  int t = (int)tC - (int)'A';
  vector<vector<int>> matrix (26, vector<int> (26, 0));
  set<int> nodes;
  for(int i = 0; i < M; i++){
    char aC, bC; inFile >> aC >> bC;
    int flow; inFile >> flow;
    int a = (int)aC - (int)'A';
    int b = (int)bC - (int)'A';
    nodes.insert(a);
    nodes.insert(b);
    matrix[a][b] = flow;
    matrix[b][a] = flow;
  }
  /*if(countComponents(matrix) != 1){
    cerr << "invalid graph" << endl;
    return;
  }*/
  vector<vector<int>> original = matrix;
  //cout << countComponents(original) - (26 - nodes.size()) << endl;
  vector<pair<int, int>> min_cut = minCut(matrix, s, t);
  int min_cut_size = 0;
  for(int i = 0; i < min_cut.size(); i++){
    int flow = original[min_cut[i].first][min_cut[i].second];
    min_cut_size += flow;
  }
  for(int i = 0; i < min_cut.size(); i++){
    original[min_cut[i].first][min_cut[i].second] = 0;
    original[min_cut[i].second][min_cut[i].first] = 0;
  }
  //printMinCut(min_cut);
  /*bool isPartitioned = isBinaryPartition(original);
  cout << isPartitioned << endl;*/
  int amountOfComponents = countComponents(original) - (26 - nodes.size());
  ifstream answerFile (answerfileString);
  int max_flow; answerFile >> max_flow;
  if(max_flow == min_cut_size && amountOfComponents == 2) cout << "\033[1;32mpassed test \033[0m" << test << " with answer = " << min_cut_size << " and correct = " << max_flow << endl;
  else cout << "\033[1;31mfailed test \033[0m" << test << " with answer = " << min_cut_size << " and correct = " << max_flow << endl;
}

void validateAllMinCut(){
  cout << "MIN CUT VALIDATION" << endl;
  for(int i = 1; i <= 10; i++){
    validateMinCut("./tests/maxflowtest/in" + to_string(i) + ".txt", "./tests/maxflowtest/corr" + to_string(i) + ".txt", i);
  }
}

/*int main(){
  validateAllMaxFlowTestData();
  cout << endl;
  validateAllMinCut();
}*/
