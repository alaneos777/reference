#include <bits/stdc++.h>
using namespace std;

template<typename T>
struct flowEdge{
	int dest;
	T flow, capacity, cost;
	flowEdge *res;

	flowEdge(): dest(0), flow(0), capacity(0), cost(0), res(NULL){}
	flowEdge(int dest, T flow, T capacity, T cost = 0): dest(dest), flow(flow), capacity(capacity), cost(cost), res(NULL){}

	void addFlow(T flow){
		this->flow += flow;
		this->res->flow -= flow;
	}
};

template<typename T>
struct flowGraph{
	T inf = numeric_limits<T>::max();
	vector<vector<flowEdge<T>*>> adjList;
	vector<int> dist, pos;
	int V;
	flowGraph(int V): V(V), adjList(V), dist(V), pos(V){}
	~flowGraph(){
		for(int i = 0; i < V; ++i)
			for(int j = 0; j < adjList[i].size(); ++j)
				delete adjList[i][j];
	}
	void addEdge(int u, int v, T capacity, T cost = 0){
		flowEdge<T> *uv = new flowEdge<T>(v, 0, capacity, cost);
		flowEdge<T> *vu = new flowEdge<T>(u, capacity, capacity, -cost);
		uv->res = vu;
		vu->res = uv;
		adjList[u].push_back(uv);
		adjList[v].push_back(vu);
	}

	//Maximun Flow using Dinic Algorithm O(EV^2)
	T blockingFlow(int u, int t, T flow){
		if(u == t) return flow;
		for(int &i = pos[u]; i < adjList[u].size(); ++i){
			flowEdge<T> *v = adjList[u][i];
			if(v->capacity > v->flow && dist[u] + 1 == dist[v->dest]){
				T fv = blockingFlow(v->dest, t, min(flow, v->capacity - v->flow));
				if(fv > 0){
					v->addFlow(fv);
					return fv;
				}
			}
		}
		return 0;
	}
	T dinic(int s, int t){
		T maxFlow = 0;
		dist[t] = 0;
		while(dist[t] != -1){
			fill(dist.begin(), dist.end(), -1);
			queue<int> Q;
			Q.push(s);
			dist[s] = 0;
			while(!Q.empty()){
				int u = Q.front(); Q.pop();
				for(flowEdge<T> *v : adjList[u]){
					if(dist[v->dest] == -1 && v->flow != v->capacity){
						dist[v->dest] = dist[u] + 1;
						Q.push(v->dest);
					}
				}
			}
			if(dist[t] != -1){
				T f;
				fill(pos.begin(), pos.end(), 0);
				while(f = blockingFlow(s, t, inf))
					maxFlow += f;
			}
		}
		return maxFlow;
	}

	//Maximun Flow using Edmonds-Karp Algorithm O(VE^2)
	T edmondsKarp(int s, int t){
		T maxFlow = 0;
		vector<flowEdge<T>*> parent(V);
		while(true){
			fill(parent.begin(), parent.end(), nullptr);
			queue<int> Q;
			Q.push(s);
			while(!Q.empty() && !parent[t]){
				int u = Q.front(); Q.pop();
				for(flowEdge<T> *v : adjList[u]){
					if(!parent[v->dest] && v->capacity > v->flow){
						parent[v->dest] = v;
						Q.push(v->dest);
					}
				}
			}
			if(!parent[t]) break;
			T f = inf;
			for(int u = t; u != s; u = parent[u]->res->dest)
				f = min(f, parent[u]->capacity - parent[u]->flow);
			for(int u = t; u != s; u = parent[u]->res->dest)
				parent[u]->addFlow(f);
			maxFlow += f;
		}
		return maxFlow;
	}

	//Max Flow Min Cost
	pair<T, T> maxFlowMinCost(int s, int t){
		vector<bool> inQueue(V);
		vector<T> distance(V), cap(V);
		vector<flowEdge<T>*> parent(V);
		T maxFlow = 0, minCost = 0;
		while(true){
			fill(distance.begin(), distance.end(), inf);
			fill(parent.begin(), parent.end(), nullptr);
			fill(cap.begin(), cap.end(), 0);
			distance[s] = 0;
			cap[s] = inf;
			queue<int> Q;
			Q.push(s);
			while(!Q.empty()){
				int u = Q.front(); Q.pop(); inQueue[u] = 0;
				for(flowEdge<T> *v : adjList[u]){
					if(v->capacity > v->flow && distance[v->dest] > distance[u] + v->cost){
						distance[v->dest] = distance[u] + v->cost;
						parent[v->dest] = v;
						cap[v->dest] = min(cap[u], v->capacity - v->flow);
						if(!inQueue[v->dest]){
							Q.push(v->dest);
							inQueue[v->dest] = true;
						}
					}
				}
			}
			if(!parent[t]) break;
			maxFlow += cap[t];
			minCost += cap[t] * distance[t];
			for(int u = t; u != s; u = parent[u]->res->dest)
				parent[u]->addFlow(cap[t]);
		}
		return {maxFlow, minCost};
	}
};

//Given a n*m cost matrix (n<=m), it finds a minimum cost assignment.
//The actual assignment is in the vector returned.
//To find the maximum, negate the values and the answer.
template<typename T>
pair<T, vector<int>> hungarian(const vector<vector<T>> & a){
	int n = a.size(), m = a[0].size();
	assert(n <= m);
	vector<int> ans(n), pa(n+1, -1), pb(m+1, -1), way(m, -1);
	vector<T> minv(m), u(n+1), v(m+1);
	vector<bool> used(m+1);
	T inf = numeric_limits<T>::max();
	for(int i = 0; i < n; ++i){
		fill(minv.begin(), minv.end(), inf);
		fill(used.begin(), used.end(), false);
		pb[m] = i;
		pa[i] = m;
		int j0 = m;
		do{
			used[j0] = true;
			int i0 = pb[j0];
			T delta = inf;
			int j1 = -1;
			for(int j = 0; j < m; ++j){
				if(used[j]) continue;
				T cur = a[i0][j] - u[i0] - v[j];
				if(cur < minv[j]){
					minv[j] = cur;
					way[j] = j0;
				}
				if(minv[j] < delta){
					delta = minv[j];
					j1 = j;
				}
			}
			for(int j = 0; j <= m; ++j){
				if(used[j]){
					u[pb[j]] += delta;
					v[j] -= delta;
				}else{
					minv[j] -= delta;
				}
			}
			j0 = j1;
		}while(pb[j0] != -1);
		do{
			int j1 = way[j0];
			pb[j0] = pb[j1];
			pa[pb[j0]] = j0;
			j0 = j1;
		}while(j0 != m);
	}
	for(int i = 0; i < n; ++i)
		ans[pb[i]] = i;
	return {-v[m], ans};
}

int main(){
	vector<vector<int>> a = {{108, 125, 150}, {150, 135, 175}, {122, 148, 250}};
	for(auto & row : a){
		for(auto & x : row){
			x = -x;
		}
	}
	auto ans = hungarian<int>(a);
	cout << ans.first << "\n";
	for(int i = 0; i < 3; ++i){
		cout << (i+1) << ", " << (ans.second[i]+1) << "\n";
	}
}