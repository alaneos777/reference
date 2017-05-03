#include <bits/stdc++.h>

using namespace std;

int inf = 1e9;

struct arista{
    int destino, costo;
    arista(int v, int w){
        destino = v, costo = w;
    }
};

struct camino{
    int costo = inf;
    list<int> vertices;
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
        for(int i = 0; i < V; i++) matriz_pesos[i][i] = 0;
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

    vector<camino> dijkstra(int origen){
        set< pair<int, int> > info;
        info.insert(make_pair(0, origen));
        vector<camino> caminos(V, camino());
        caminos[origen].costo = 0;
        while(!info.empty()){
            pair<int, int> i = *info.begin();
            info.erase(info.begin());
            int actual = i.second;
            for(size_t j = 0; j < vertices[actual].size(); j++){
                int dest = vertices[actual][j].destino;
                int nuevo = caminos[actual].costo + vertices[actual][j].costo; //i.first + vertices[actual][j].costo
                if(nuevo < caminos[dest].costo){
                    info.erase(make_pair(caminos[dest].costo, dest));
                    info.insert(make_pair(nuevo, dest));
                    caminos[dest].costo = nuevo;
                    caminos[dest].anterior = actual;
                }
            }
        }
        for(int i = 0; i < V; i++){
            int actual = i;
            while(true){
                caminos[i].vertices.push_front(actual);
                if(actual == origen) break;
                actual = caminos[actual].anterior;
            }
        }
        return caminos;
    }

    vector< vector<int> > floyd(){
        vector< vector<int> > tmp = matriz_pesos;
        for(int k = 0; k < V; k++){
            for(int i = 0; i < V; i++){
                for(int j = 0; j < V; j++){
                    tmp[i][j] = min(tmp[i][j], tmp[i][k] + tmp[k][j]);
                }
            }
        }
        return tmp;
    }

    void recorrer(int inicio, vector<bool> & pendientes){
        pendientes[inicio] = false;
        for(size_t j = 0; j < vertices[inicio].size(); j++){
            if(pendientes[vertices[inicio][j].destino]) recorrer(vertices[inicio][j].destino, pendientes);
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

    /*grafo g(9, false);
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
    g.anadir_vertice(5, 8, 2);*/

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
    g.anadir_vertice(1, 5, 15);
    g.anadir_vertice(5, 9, 4);
    g.anadir_vertice(2, 6, 13);
    g.anadir_vertice(6, 10, 9);
    g.anadir_vertice(3, 7, 2);
    g.anadir_vertice(7, 11, 1);
    g.anadir_vertice(2, 7, 6);
    g.anadir_vertice(5, 8,1);
    g.anadir_vertice(6, 11, 8);*/

    /*vector<camino> rutas = g.dijkstra(0);
    for(camino p : rutas){
        for(int i:p.vertices){
            cout << (char)(i+97) << " ";
        }
        cout << ": " << p.costo << endl;
    }
    cout << endl;
    vector< vector<int> > m = g.floyd();
    for(vector<int> fila:m){
        for(int valor:fila){
            cout << valor << " ";
        }
        cout << "\n";
    }*/

    grafo G(11, false);
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

    //cout << "k(G)=" << G.componentes() << endl;

    grafo ind = G.inducido({1, 3, 4, 8});

    set< vector<int> > ans;
    for(int a=0;a<11;a++){
        for(int b=0;b<11;b++){
            for(int c=0;c<11;c++){
                for(int d=0;d<11;d++){
                    if(a!=b && a!=c && a!=d && b!=c && b!=d && c!=d){
                        vector<int> v_tmp = {a, b, c, d};
                        grafo temp = G.inducido(v_tmp);
                        bool test = true;
                        for(int i=0;i<4;i++){
                            for(int j=0;j<4;j++){
                                test = test && (ind.matriz[i][j] == temp.matriz[i][j]);
                            }
                        }
                        if(test){
                            sort(v_tmp.begin(), v_tmp.begin()+4);
                            ans.insert(v_tmp);
                        }
                    }
                }
            }
        }
    }

    int contador = 0;
    for(vector<int> y:ans){
        contador++;
        cout << contador << ": ";
        for(int x:y) cout << (char)(x+97) << " ";
        cout << endl;
    }

    return 0;
}
