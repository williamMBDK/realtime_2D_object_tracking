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

void dfs(vector<vector<int>> &matrix, vector<int>& parents, int node){
  for(int i = 0; i < matrix.size(); i++){
    if(matrix[node][i] != 0 && parents[i] == -1){
      parents[i] = node;
      dfs(matrix, parents, i);
    }
  }
}

/*
  optimering: bruge bfs
  optimering: stoppe ligeså snart t er noget i bfs metoden
*/
vector<int> getPath(vector<vector<int>> &matrix, int s, int t){
  vector<int> parents (matrix.size(), -1);
  dfs(matrix, parents, s);
  if(parents[t] == -1) return vector<int> (0);
  vector<int> path = {t};
  int node = t;
  while(node != s){
    node = parents[node];
    path.push_back(node);
  }
  return path;
}

int maxFlow(vector<vector<int>> &matrix, int s, int t){
  int res = 0;
  int flow = -1;
  while(flow != 0){
    vector<int> path = getPath(matrix, s, t);
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

void validateTestData(string infileString, string answerfileString, int test){
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
  if(ans == corr) cout << "passed test " << test << " with answer = " << ans << " and correct = " << corr << endl;
  else cout << "failed test " << test << " with answer = " << ans  << " and correct = " << corr << endl;
}

void validateAllTestData(){
  for(int i = 1; i <= 10; i++){
    validateTestData("./tests/maxflowtest/in" + to_string(i) + ".txt", "./tests/maxflowtest/corr" + to_string(i) + ".txt", i);
  }
}

int main(){
  validateAllTestData();
}
