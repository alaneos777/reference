#include <bits/stdc++.h>

using namespace std;

int inf = numeric_limits<int>::max();

struct node{
    int value;
    node * parent;
    int rank;

    node(int value){
        this->value = value;
        this->parent = this;
        this->rank = 0;
    }
};

struct disjointSet{
    int N;
    vector<node*> nodes;

    disjointSet(int N){
        this->N = N;
        nodes.resize(N);
    }

    void makeSet(int v){
        nodes[v] = new node(v);
    }

    node * findSet(int v){
        node * current = nodes[v];
        if(current -> value != current -> parent -> value){
            current -> parent = findSet(current -> parent -> value);
        }
        return current -> parent;
    }

    node * unionSet(int a, int b){
        node * p1 = findSet(a);
        node * p2 = findSet(b);
        if(p1 -> rank >= p2 -> rank){
            p2 -> parent = p1;
            if(p1 -> rank == p2 -> rank){
                p1 -> rank++;
            }
            return p1;
        }else{
            p1 -> parent = p2;
            return p2;
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
        adjList[source].push_back(edge(dest, cost));
        edges.push_back(edge(source, dest, cost));
        adjMatrix[source][dest] = true;
        costMatrix[source][dest] = cost;
        if(!dir){
            adjList[dest].push_back(edge(source, cost));
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
            int source = cola.top().dest;
            cola.pop();
            if(visited[source]) continue;
            visited[source] = true;
            for(edge & current : adjList[source]){
                int dest = current.dest;
                if(visited[dest]) continue;
                int nuevo = paths[source].cost + current.cost;
                if(nuevo == paths[dest].cost && paths[source].size + 1 < paths[dest].size){
                    paths[dest].previous = source;
                    paths[dest].size = paths[source].size + 1;
                }else if(nuevo < paths[dest].cost){
                    paths[dest].previous = source;
                    paths[dest].size = paths[source].size + 1;
                    cola.push(edge(dest, nuevo));
                    paths[dest].cost = nuevo;
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
            for(int source = 0; source < V; source++){
                if(paths[source].cost == inf) continue;
                for(edge & current : adjList[source]){
                    int dest = current.dest;
                    int nuevo = paths[source].cost + current.cost;
                    if(nuevo == paths[dest].cost && paths[source].size + 1 < paths[dest].size){
                        paths[dest].previous = source;
                        paths[dest].size = paths[source].size + 1;
                    }else if(nuevo < paths[dest].cost){
                        if(j == V){
                            cout << "Ciclo negativo\n";
                            return {};
                        }
                        paths[dest].previous = source;
                        paths[dest].size = paths[source].size + 1;
                        paths[dest].cost = nuevo;
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
        for(int k = 0; k < V; k++){
            for(int i = 0; i < V; i++){
                for(int j = 0; j < V; j++){
                    if(tmp[i][k] != inf && tmp[k][j] != inf)
                        tmp[i][j] = min(tmp[i][j], tmp[i][k] + tmp[k][j]);
                }
            }
        }
        return tmp;
    }

    vector<vector<bool>> transitiveClosure(){
        vector<vector<bool>> tmp = adjMatrix;
        for(int k = 0; k < V; k++){
            for(int i = 0; i < V; i++){
                for(int j = 0; j < V; j++){
                    tmp[i][j] = tmp[i][j] || (tmp[i][k] && tmp[k][j]);
                }
            }
        }
        return tmp;
    }

    void DFS(int start, int source, vector<vector<bool>> & tmp){
        for(edge & current : adjList[source]){
            if(!tmp[start][current.dest]){
                tmp[start][current.dest] = true;
                DFS(start, current.dest, tmp);
            }
        }
    }

    vector<vector<bool>> transitiveClosureDFS(){
        vector<vector<bool>> tmp(V, vector<bool>(V, false));
        for(int i = 0; i < V; i++){
            DFS(i, i, tmp);
        }
        return tmp;
    }

    void DFS(int source, vector<bool> & visited){
        visited[source] = true;
        for(edge & current : adjList[source]){
            if(!visited[current.dest]){
                DFS(current.dest, visited);
            }
        }
    }

    bool isBipartite(){
        vector<int> side(V, -1);
        bool is_bipartite = true;
        queue<int> q;
        for (int st = 0; st < V; ++st) {
            if (side[st] == -1) {
                q.push(st);
                side[st] = 0;
                while (!q.empty()) {
                    int v = q.front();
                    q.pop();
                    for (edge & u : adjList[v]) {
                        if (side[u.dest] == -1) {
                            side[u.dest] = side[v] ^ 1;
                            q.push(u.dest);
                        } else {
                            is_bipartite &= side[u.dest] != side[v];
                        }
                    }
                }
            }
        }
        return is_bipartite;
    }

    int components(){
        int ans = 0;
        vector<bool> visited(V, false);
        for(int i = 0; i < V; i++){
            if(!visited[i]){
                DFS(i, visited);
                ans++;
            }
        }
        return ans;
    }

    vector<edge> kruskal(){
        sort(edges.begin(), edges.end());
        vector<edge> MST;
        disjointSet DS(V);
        for(int i = 0; i < V; i++)
            DS.makeSet(i);
        int i = 0;
        while(i < edges.size() && MST.size() < V - 1){
            edge current = edges[i++];
            if(DS.findSet(current.source) != DS.findSet(current.dest)){
                MST.push_back(current);
                DS.unionSet(current.source, current.dest);
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
