#include <bits/stdc++.h>

using namespace std;

int inf = numeric_limits<int>::max();

struct disjointSet{
    int N;
    vector<short int> rank;
    vector<int> parent;

    disjointSet(int N){
        this->N = N;
        parent.resize(N);
        rank.resize(N);
    }

    void makeSet(int v){
        parent[v] = v;
    }

    int findSet(int v){
        if(v == parent[v]) return v;
        return parent[v] = findSet(parent[v]);
    }

    void unionSet(int a, int b){
        a = findSet(a);
        b = findSet(b);
        if(a == b) return;
        if(rank[a] < rank[b]){
            parent[a] = b;
        }else{
            parent[b] = a;
            if(rank[a] == rank[b]){
                ++rank[a];
            }
        }
    }
};

struct edge{
    int source, dest, cost;
    edge(){
        this->source = this->dest = this->cost = 0;
    }
    edge(int dest, int cost){
        this->dest = dest;
        this->cost = cost;
    }
    edge(int source, int dest, int cost){
        this->source = source;
        this->dest = dest;
        this->cost = cost;
    }
    bool operator==(const edge & b) const{
        return source == b.source && dest == b.dest && cost == b.cost;
    }
    bool operator<(const edge & b) const{
        if(cost == b.cost){
            if(dest == b.dest){
                return source < b.source;
            }else{
                return dest < b.dest;
            }
        }else{
            return cost < b.cost;
        }
    }
    bool operator>(const edge & b) const{
        if(cost == b.cost){
            if(dest == b.dest){
                return source > b.source;
            }else{
                return dest > b.dest;
            }
        }else{
            return cost > b.cost;
        }
    }
};

struct path{
    int cost = inf;
    list<int> vertices;
    int size = 1;
    int previous = -1;
};

struct grafo{
    vector<vector<edge>> adjList;
    vector<vector<bool>> adjMatrix;
    vector<vector<int>> costMatrix;
    vector<edge> edges;
    int V = 0;
    bool dir = false;

    grafo(int n, bool dirigido){
        V = n;
        dir = dirigido;
        adjList.resize(V, vector<edge>());
        edges.resize(V);
        adjMatrix.resize(V, vector<bool>(V, false));
        costMatrix.resize(V, vector<int>(V, inf));
        for(int i = 0; i < V; i++)
            costMatrix[i][i] = 0;
    }

    void anadir_vertice(int source, int dest, int cost){
        adjList[source].push_back(edge(source, dest, cost));
        edges.push_back(edge(source, dest, cost));
        adjMatrix[source][dest] = true;
        costMatrix[source][dest] = cost;
        if(!dir){
            adjList[dest].push_back(edge(dest, source, cost));
            adjMatrix[dest][source] = true;
            costMatrix[dest][source] = cost;
        }
    }

    struct comparador{
        bool operator() (const edge & a, const edge & b) const{
            return a > b;
        }
    };

    void buildPaths(vector<path> & paths){
        for(int i = 0; i < V; i++){
            int actual = i;
            for(int j = 0; j < paths[i].size; j++){
                paths[i].vertices.push_front(actual);
                actual = paths[actual].previous;
            }
        }
    }

    vector<path> dijkstra(int start){
        priority_queue<edge, vector<edge>, comparador> cola;
        vector<path> paths(V, path());
        vector<bool> visited(V, false);
        cola.push(edge(start, 0));
        paths[start].cost = 0;
        while(!cola.empty()){
            int u = cola.top().dest;
            cola.pop();
            if(visited[u]) continue;
            visited[u] = true;
            for(edge & current : adjList[u]){
                int v = current.dest;
                if(visited[v]) continue;
                int nuevo = paths[u].cost + current.cost;
                if(nuevo == paths[v].cost && paths[u].size + 1 < paths[v].size){
                    paths[v].previous = u;
                    paths[v].size = paths[u].size + 1;
                }else if(nuevo < paths[v].cost){
                    paths[v].previous = u;
                    paths[v].size = paths[u].size + 1;
                    cola.push(edge(v, nuevo));
                    paths[v].cost = nuevo;
                }
            }
        }
        buildPaths(paths);
        return paths;
    }

    vector<path> bellmanFord(int start){
        vector<path> paths(V, path());
        paths[start].cost = 0;
        bool cambio = true;
        int j = 1;
        while(cambio){
            cambio = false;
            for(int u = 0; u < V; ++u){
                if(paths[u].cost == inf) continue;
                for(edge & current : adjList[u]){
                    int v = current.dest;
                    int nuevo = paths[u].cost + current.cost;
                    if(nuevo == paths[v].cost && paths[u].size + 1 < paths[v].size){
                        paths[v].previous = u;
                        paths[v].size = paths[u].size + 1;
                    }else if(nuevo < paths[v].cost){
                        if(j == V){
                            cout << "Ciclo negativo\n";
                            return {};
                        }
                        paths[v].previous = u;
                        paths[v].size = paths[u].size + 1;
                        paths[v].cost = nuevo;
                        cambio = true;
                    }
                }
            }
            j++;
        }
        buildPaths(paths);
        return paths;
    }

    vector<vector<int>> floyd(){
        vector<vector<int>> tmp = costMatrix;
        for(int k = 0; k < V; ++k)
            for(int i = 0; i < V; ++i)
                for(int j = 0; j < V; ++j)
                    if(tmp[i][k] != inf && tmp[k][j] != inf)
                        tmp[i][j] = min(tmp[i][j], tmp[i][k] + tmp[k][j]);
        return tmp;
    }

    vector<vector<bool>> transitiveClosure(){
        vector<vector<bool>> tmp = adjMatrix;
        for(int k = 0; k < V; ++k)
            for(int i = 0; i < V; ++i)
                for(int j = 0; j < V; ++j)
                    tmp[i][j] = tmp[i][j] || (tmp[i][k] && tmp[k][j]);
        return tmp;
    }

    void DFSClosure(int start, int source, vector<vector<bool>> & tmp){
        for(edge & current : adjList[source]){
            int v = current.dest;
            if(!tmp[start][v]){
                tmp[start][v] = true;
                DFSClosure(start, v, tmp);
            }
        }
    }

    vector<vector<bool>> transitiveClosureDFS(){
        vector<vector<bool>> tmp(V, vector<bool>(V, false));
        for(int u = 0; u < V; u++)
            DFSClosure(u, u, tmp);
        return tmp;
    }

    bool isBipartite(){
        vector<int> side(V, -1);
        queue<int> q;
        for (int st = 0; st < V; ++st) {
            if(side[st] != -1) continue;
            q.push(st);
            side[st] = 0;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (edge & current : adjList[u]) {
                    int v = current.dest;
                    if (side[v] == -1) {
                        side[v] = side[u] ^ 1;
                        q.push(v);
                    } else {
                        if(side[v] == side[u]) return false;
                    }
                }
            }
        }
        return true;
    }

    vector<int> topologicalSort(){
        vector<int> order;
        int visited = 0;
        vector<int> indegree(V);
        for(auto & node : adjList){
            for(edge & current : node){
                int v = current.dest;
                ++indegree[v];
            }
        }
        queue<int> Q;
        for(int i = 0; i < V; ++i){
            if(indegree[i] == 0) Q.push(i);
        }
        while(!Q.empty()){
            int source = Q.front();
            Q.pop();
            order.push_back(source);
            ++visited;
            for(edge & current : adjList[source]){
                int v = current.dest;
                --indegree[v];
                if(indegree[v] == 0) Q.push(v);
            }
        }
        if(visited == V) return order;
        else return {};
    }

    void DFSCycle(int u, vector<int> & color, bool & cycle){
        if(color[u] == 0){
            color[u] = 1;
            for(edge & current : adjList[u]){
                int v = current.dest;
                DFSCycle(v, color, cycle);
            }
            color[u] = 2;
        }else if(color[u] == 1){
            cycle = true;
        }
    }

    bool DFSCycle(int u, vector<bool> & visited, int source){
        visited[u] = true;
        for(edge & current : adjList[u]){
            int v = current.dest;
            if(!visited[v]){
                if(DFSCycle(v, visited, u)) return true;
            }else if(v != source){
                return true;
            }
            return false;
        }
    }

    bool hasCycle(){
        if(dir){
            vector<int> color(V);
            bool cycle = false;
            for(int u = 0; u < V; ++u){
                DFSCycle(u, color, cycle);
                if(cycle) return true;
            }
            return false;
        }else{
            vector<bool> visited(V, false);
            for(int u = 0; u < V; ++u){
                if(!visited[u] && DFSCycle(u, visited, -1)) return true;
            }
            return false;
        }
    }

    int articulationBridges(int u, int p, vector<int> & low, vector<int> & label, int & time, vector<bool> & points, vector<edge> & bridges){
        label[u] = low[u] = ++time;
        int hijos = 0, ret = 0;
        for(edge & current : adjList[u]){
            int v = current.dest;
            if(v == p && !ret++) continue;
            if(!label[v]){
                ++hijos;
                articulationBridges(v, u, low, label, time, points, bridges);
                if(label[u] <= low[v])
                    points[u] = true;
                else if(label[u] < low[v])
                    bridges.push_back(current);
                low[u] = min(low[u], low[v]);
            }
            low[u] = min(low[u], label[v]);
        }
        return hijos;
    }

    pair<vector<bool>, vector<edge>> articulationBridges(){
        vector<int> low(V), label(V);
        vector<bool> points(V);
        vector<edge> bridges;
        int time = 0;
        for(int u = 0; u < V; ++u)
            if(!label[u])
                points[u] = articulationBridges(u, -1, low, label, time, points, bridges) > 1;
        return make_pair(points, bridges);
    }

    void scc(int u, vector<int> & low, vector<int> & label, int & time, vector<vector<int>> & ans, stack<int> & S){
        label[u] = low[u] = ++time;
        S.push(u);
        for(edge & current : adjList[u]){
            int v = current.dest;
            if(!label[v]) scc(v, low, label, time, ans, S);
            low[u] = min(low[u], low[v]);
        }
        if(label[u] == low[u]){
            vector<int> comp;
            while(S.top() != u){
                comp.push_back(S.top());
                low[S.top()] = V + 1;
                S.pop();
            }
            comp.push_back(S.top());
            S.pop();
            low[u] = V + 1;
        }
    }

    vector<vector<int>> scc(){
        vector<int> low(V), label(V);
        int time = 0;
        vector<vector<int>> ans;
        stack<int> S;
        for(int u = 0; u < V; ++u)
            if(!label[u]) scc(u, low, label, time, ans, S);
        return ans;
    }

    void DFSComponents(int source, vector<bool> & visited){
        visited[source] = true;
        for(edge & current : adjList[source]){
            int v = current.dest;
            if(!visited[v]) DFSComponents(v, visited);
        }
    }

    int components(){
        int ans = 0;
        vector<bool> visited(V, false);
        for(int u = 0; u < V; ++u){
            if(!visited[u]){
                DFSComponents(u, visited);
                ++ans;
            }
        }
        return ans;
    }

    vector<edge> kruskal(){
        sort(edges.begin(), edges.end());
        vector<edge> MST;
        disjointSet DS(V);
        for(int u = 0; u < V; ++u)
            DS.makeSet(u);
        int i = 0;
        while(i < edges.size() && MST.size() < V - 1){
            edge current = edges[i++];
            int u = current.source, v = current.dest;
            if(DS.findSet(u) != DS.findSet(v)){
                MST.push_back(current);
                DS.unionSet(u, v);
            }
        }
        return MST;
    }

    grafo inducido(vector<int> nuevos){
        int tam = nuevos.size();
        grafo ans(tam, true);
        for(int i = 0; i < tam; i++){
            int v1 = nuevos[i];
            for(int j = 0; j < tam; j++){
                int v2 = nuevos[j];
                if(adjMatrix[v1][v2]) ans.anadir_vertice(i, j, costMatrix[v1][v2]);
            }
        }
        return ans;
    }
};

int main()
{
    /*grafo g(6, false);
    g.anadir_vertice(0, 1, 5);
    g.anadir_vertice(1, 2, 1);
    g.anadir_vertice(2, 5, 1);
    g.anadir_vertice(5, 4, 7);
    g.anadir_vertice(4, 3, 2);
    g.anadir_vertice(3, 0, 1);
    g.anadir_vertice(1, 3, 6);
    g.anadir_vertice(1, 4, 1);
    g.anadir_vertice(1, 5, 2);*/

    /*grafo g(8, false);
    g.anadir_vertice(0, 3, 1);
    g.anadir_vertice(3, 6, 9);
    g.anadir_vertice(6, 7, 1);
    g.anadir_vertice(7, 4, 2);
    g.anadir_vertice(4, 1, 6);
    g.anadir_vertice(1, 0, 7);
    g.anadir_vertice(0, 2, 5);
    g.anadir_vertice(1, 2, 2);
    g.anadir_vertice(2, 4, 1);
    g.anadir_vertice(2, 5, 3);
    g.anadir_vertice(2, 6, 4);
    g.anadir_vertice(5, 6, 2);
    g.anadir_vertice(5, 7, 1);*/

    /*grafo g(9, false);
    g.anadir_vertice(0, 1, 6);
    g.anadir_vertice(1, 2, 5);
    g.anadir_vertice(2, 5, 2);
    g.anadir_vertice(5, 8, 4);
    g.anadir_vertice(8, 7, 1);
    g.anadir_vertice(7, 6, 8);
    g.anadir_vertice(6, 3, 1);
    g.anadir_vertice(3, 0, 1);
    g.anadir_vertice(1, 4, 1);
    g.anadir_vertice(4, 7, 1);
    g.anadir_vertice(3, 4, 7);
    g.anadir_vertice(0, 4, 5);
    g.anadir_vertice(1, 5, 1);
    g.anadir_vertice(4, 6, 1);
    g.anadir_vertice(5, 7, 1);*/

    /*grafo g(9, false);
    g.anadir_vertice(0, 1, 4);
    g.anadir_vertice(1, 2, 8);
    g.anadir_vertice(2, 3, 7);
    g.anadir_vertice(3, 4, 9);
    g.anadir_vertice(4, 5, 10);
    g.anadir_vertice(5, 6, 2);
    g.anadir_vertice(6, 7, 1);
    g.anadir_vertice(7, 0, 8);
    g.anadir_vertice(1, 7, 11);
    g.anadir_vertice(2, 8, 2);
    g.anadir_vertice(2, 5, 4);
    g.anadir_vertice(8, 7, 7);
    g.anadir_vertice(8, 6, 6);
    g.anadir_vertice(3, 5, 14);*/

    grafo g(9, false);
    g.anadir_vertice(0, 1, 3);
    g.anadir_vertice(1, 2, 7);
    g.anadir_vertice(2, 3, 1);
    g.anadir_vertice(3, 8, 5);
    g.anadir_vertice(0, 6, 5);
    g.anadir_vertice(6, 7, 2);
    g.anadir_vertice(7, 8, 4);
    g.anadir_vertice(0, 4, 7);
    g.anadir_vertice(4, 5, 1);
    g.anadir_vertice(1, 4, 1);
    g.anadir_vertice(4, 2, 2);
    g.anadir_vertice(2, 5, 2);
    g.anadir_vertice(5, 3, 3);
    g.anadir_vertice(6, 4, 3);
    g.anadir_vertice(4, 7, 3);
    g.anadir_vertice(7, 5, 3);
    g.anadir_vertice(5, 8, 2);

    /*grafo g(6, true);
    g.anadir_vertice(0, 1, 1);
    g.anadir_vertice(1, 2, 1);
    g.anadir_vertice(2, 3, 1);
    g.anadir_vertice(3, 4, 1);
    g.anadir_vertice(3, 1, 1);*/

    /*grafo g(12, false);
    g.anadir_vertice(0, 1, 1);
    g.anadir_vertice(1, 2, 2);
    g.anadir_vertice(2, 3, 1);
    g.anadir_vertice(4, 5, 8);
    g.anadir_vertice(5, 6, 1);
    g.anadir_vertice(6, 7, 10);
    g.anadir_vertice(8, 9, 1);
    g.anadir_vertice(9, 10, 2);
    g.anadir_vertice(10, 11, 1);
    g.anadir_vertice(0, 4, 9);
    g.anadir_vertice(4, 8, 7);
    g.anadir_vertice(1, 5, 11);
    g.anadir_vertice(5, 9, 4);
    g.anadir_vertice(2, 6, 13);
    g.anadir_vertice(6, 10, 9);
    g.anadir_vertice(3, 7, 2);
    g.anadir_vertice(7, 11, 1);
    g.anadir_vertice(2, 7, 6);
    g.anadir_vertice(5, 8, 1);
    g.anadir_vertice(6, 11, 8);*/

    /*grafo g(6, false);
    g.anadir_vertice(0, 1, 1);
    g.anadir_vertice(1, 2, 1);
    g.anadir_vertice(2, 3, 1);
    g.anadir_vertice(3, 5, 2);
    g.anadir_vertice(0, 4, 3);
    g.anadir_vertice(4, 5, 2);*/

    vector<path> rutas = g.dijkstra(0);
    for(path & p : rutas){
        cout << p.size << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.cost << "\n";
    }
    cout << "\n";

    rutas = g.bellmanFord(0);
    for(path & p : rutas){
        cout << p.size << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.cost << "\n";
    }
    cout << "\n";

    vector< vector<int> > m = g.floyd();
    for(vector<int> & fila : m){
        for(int & valor : fila){
            cout << valor << " ";
        }
        cout << "\n";
    }
    cout << "\n";

    vector<vector<bool>> t = g.transitiveClosure();
    for(vector<bool> & fila : t){
        for(bool valor : fila){
            cout << valor << " ";
        }
        cout << "\n";
    }
    cout << "\n";

    t = g.transitiveClosureDFS();
    for(vector<bool> & fila : t){
        for(bool valor : fila){
            cout << valor << " ";
        }
        cout << "\n";
    }
    cout << "\n";

    cout << "Componentes: " << g.components() << "\n\nMST:\n";

    vector<edge> MST = g.kruskal();
    for(edge & current : MST){
        cout << current.source << " <-> " << current.dest << "; w=" << current.cost << "\n";
    }
    cout << "\n";

    /*grafo G(11, false);
    G.anadir_vertice(0, 4, 1);
    G.anadir_vertice(0, 5, 1);
    G.anadir_vertice(0, 6, 1);
    G.anadir_vertice(1, 3, 1);
    G.anadir_vertice(1, 4, 1);
    G.anadir_vertice(1, 5, 1);
    G.anadir_vertice(2, 5, 1);
    G.anadir_vertice(2, 6, 1);
    G.anadir_vertice(2, 7, 1);
    G.anadir_vertice(3, 4, 1);
    G.anadir_vertice(3, 8, 1);
    G.anadir_vertice(4, 5, 1);
    G.anadir_vertice(4, 8, 1);
    G.anadir_vertice(4, 10, 1);
    G.anadir_vertice(5, 6, 1);
    G.anadir_vertice(5, 10, 1);
    G.anadir_vertice(6, 7, 1);
    G.anadir_vertice(6, 9, 1);
    G.anadir_vertice(6, 10, 1);
    G.anadir_vertice(7, 9, 1);
    G.anadir_vertice(8, 10, 1);
    G.anadir_vertice(9, 10, 1);

    grafo G2(5, true);
    G2.anadir_vertice(0, 1, 4);
    G2.anadir_vertice(0, 2, 2);
    G2.anadir_vertice(1, 2, 3);
    G2.anadir_vertice(2, 1, 1);
    G2.anadir_vertice(1, 3, 2);
    G2.anadir_vertice(2, 4, 5);
    G2.anadir_vertice(1, 4, 3);
    G2.anadir_vertice(2, 3, 4);
    G2.anadir_vertice(4, 3, -5);

    /*grafo G2(5, true);
    G2.anadir_vertice(0, 1, 2);
    G2.anadir_vertice(1, 2, 2);
    G2.anadir_vertice(2, 3, -4);
    G2.anadir_vertice(3, 4, 3);
    G2.anadir_vertice(1, 3, 1);
    G2.anadir_vertice(3, 1, 1);*/

    /*vector<path> rutas = G2.bellmanFord(0);
    for(path & p : rutas){
        cout << p.size << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.cost << endl;
    }
    cout << endl;

    vector< vector<int> > m = G2.floyd();
    for(vector<int> & fila : m){
        for(int & valor : fila){
            cout << valor << " ";
        }
        cout << "\n";
    }*/

    return 0;
}
