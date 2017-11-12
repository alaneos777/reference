#include <bits/stdc++.h>

using namespace std;

int inf = numeric_limits<int>::max();

struct arista{
    int v, w;
    arista(int _v, int _w){
        v = _v, w = _w;
    }
};

struct camino{
    int costo = inf;
    list<int> vertices;
    int tamano = 1;
    int anterior;
};

struct grafo{
    vector< vector<arista> > vertices;
    vector< vector<int> > matriz, matriz_pesos;
    int V = 0;
    bool dir = false;

    grafo(int n, bool dirigido){
        V = n;
        dir = dirigido;
        vertices.resize(V, vector<arista>());
        matriz.resize(V, vector<int>(V, 0));
        matriz_pesos.resize(V, vector<int>(V, inf));
        for(int i = 0; i < V; i++)
            matriz_pesos[i][i] = 0;
    }

    void anadir_vertice(int u, int v, int w){
        vertices[u].push_back(arista(v, w));
        matriz[u][v] = 1;
        matriz_pesos[u][v] = w;
        if(!dir){
            vertices[v].push_back(arista(u, w));
            matriz[v][u] = 1;
            matriz_pesos[v][u] = w;
        }
    }

    struct comparador{
        bool operator() (const arista & a, const arista & b) const{
            if(a.w == b.w){
                return a.v > b.v;
            }else{
                return a.w > b.w;
            }
        }
    };

    void reconstruirCaminos(int origen, vector<camino> & caminos){
        for(int i = 0; i < V; i++){
            int actual = i;
            while(true){
                caminos[i].vertices.push_front(actual);
                if(actual == origen) break;
                actual = caminos[actual].anterior;
            }
        }
    }

    vector<camino> dijkstra(int origen){
        priority_queue<arista, vector<arista>, comparador> cola;
        vector<camino> caminos(V, camino());
        vector<bool> visitados(V, false);
        cola.push(arista(origen, 0));
        caminos[origen].costo = 0;
        while(!cola.empty()){
            arista actual = cola.top();
            cola.pop();
            if(visitados[actual.v]) continue;
            visitados[actual.v] = true;
            for(arista & dest : vertices[actual.v]){
                if(visitados[dest.v]) continue;
                int nuevo = caminos[actual.v].costo + dest.w; //actual.w + dest.w;
                if(nuevo == caminos[dest.v].costo && caminos[actual.v].tamano + 1 < caminos[dest.v].tamano){
                    caminos[dest.v].anterior = actual.v;
                    caminos[dest.v].tamano = caminos[actual.v].tamano + 1;
                }else if(nuevo < caminos[dest.v].costo){
                    caminos[dest.v].anterior = actual.v;
                    caminos[dest.v].tamano = caminos[actual.v].tamano + 1;
                    cola.push(arista(dest.v, nuevo));
                    caminos[dest.v].costo = nuevo;
                }
            }
        }
        reconstruirCaminos(origen, caminos);
        return caminos;
    }

    vector<camino> bellmanFord(int origen){
        vector<camino> caminos(V, camino());
        caminos[origen].costo = 0;
        bool cambio = true;
        int j = 1;
        while(cambio){
            cambio = false;
            for(int i = 0; i < V; i++){
                for(arista & dest : vertices[i]){
                    if(caminos[i].costo == inf) continue;
                    int nuevo = caminos[i].costo + dest.w;
                    if(nuevo == caminos[dest.v].costo && caminos[i].tamano + 1 < caminos[dest.v].tamano){
                        caminos[dest.v].anterior = i;
                        caminos[dest.v].tamano = caminos[i].tamano + 1;
                    }else if(nuevo < caminos[dest.v].costo){
                        if(j == V){
                            cout << "Ciclo negativo\n";
                            return {};
                        }
                        caminos[dest.v].anterior = i;
                        caminos[dest.v].tamano = caminos[i].tamano + 1;
                        caminos[dest.v].costo = nuevo;
                        cambio = true;
                    }
                }
            }
            j++;
        }
        reconstruirCaminos(origen, caminos);
        return caminos;
    }

    vector< vector<int> > floyd(){
        vector< vector<int> > tmp = matriz_pesos;
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

    void recorrer(int inicio, vector<bool> & pendientes){
        pendientes[inicio] = false;
        for(arista & dest : vertices[inicio]){
            if(pendientes[dest.v])
                recorrer(dest.v, pendientes);
        }
    }

    int componentes(){
        int ans = 0;
        vector<bool> pendientes(V, true);
        for(int i = 0; i < V; i++){
            if(pendientes[i]){
                recorrer(i, pendientes);
                ans++;
            }
        }
        return ans;
    }

    grafo inducido(vector<int> nuevos){
        int tam = nuevos.size();
        grafo ans(tam, true);
        for(int i = 0; i < tam; i++){
            int v1 = nuevos[i];
            for(int j = 0; j < tam; j++){
                int v2 = nuevos[j];
                if(matriz[v1][v2] == 1) ans.anadir_vertice(i, j, matriz_pesos[v1][v2]);
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

    vector<camino> rutas = g.dijkstra(8);
    for(camino & p : rutas){
        cout << p.tamano << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.costo << endl;
    }
    cout << "\n";

    rutas = g.bellmanFord(8);
    for(camino & p : rutas){
        cout << p.tamano << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.costo << endl;
    }
    cout << endl;

    vector< vector<int> > m = g.floyd();
    for(vector<int> & fila : m){
        for(int & valor : fila){
            cout << valor << " ";
        }
        cout << "\n";
    }

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

    cout << "k(G)=" << G.componentes() << endl;*/

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

    /*vector<camino> rutas = G2.bellmanFord(0);
    for(camino & p : rutas){
        cout << p.tamano << ": ";
        for(int & i : p.vertices){
            cout << i << " ";
        }
        cout << ": " << p.costo << endl;
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
