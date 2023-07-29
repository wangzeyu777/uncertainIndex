#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

// #include <map>

#include "edge.h"
#include "indexEntry.h"

using namespace std;
template<typename T>

double sqr(double t);
void handle_error(const char* msg);
typedef struct{
	int index;
	int value;
}sort_v;
bool compareVertex(sort_v a, sort_v b){
	return a.value>b.value;
}

bool compareEdge(Edge a, Edge b){
	return a.d>b.d;
}

class Graph{
    public:
        int m,n;
        string folder="../dataset/";
        string graph_file;
        vector<int> vertices;
        vector<int> degree, indegree, outdegree;
        vector<map<int, double>> inneighbors, outneighbors;
        vector<Edge> edges;
        vector<Edge> reverseEdges;

        vector<int> sortedVertices;
        vector<Edge> sortedEdges;
        vector<bool> inVertexCover;
        vector<bool> RSSvisit;
        
        vector<vector<IndexEntry>> Lf, Lb;

        double time_R, time_ERD, time_DCR;
        double BSS_R, BSS_ERD, BSS_DCR;
        double RSS_R, RSS_ERD, RSS_DCR;

        void readNM();
        void readGraph();
        void statistics();
        void degreeSort(int k, string str);
        void edgeDegree();
        void edgeSort(int k);
        vector<Edge> randomEdge();
    Graph(string graph_file): graph_file(graph_file)
    {
        readNM();
        for(int i=0;i<n;i++){
            vertices.push_back(i);
        }
        // edges = vector<Edge>(m, Edge());	
        // reverseEdges = vector<Edge>(m, Edge());	
        // sortedEdges = vector<Edge>(m, Edge());	
        sortedVertices = vector<int>(n, 0);
        degree = vector<int>(n, 0);
        indegree = vector<int>(n, 0);
        outdegree = vector<int>(n, 0);	
        inneighbors = vector<map<int, double>>(n, map<int, double>());
        outneighbors = vector<map<int, double>>(n, map<int, double>());
        inVertexCover = vector<bool>(n, false);
        RSSvisit = vector<bool>(m, false);
        Lf = vector<vector<IndexEntry>>(n, vector<IndexEntry>());
        Lb = vector<vector<IndexEntry>>(n, vector<IndexEntry>());
        time_R = 0; time_ERD = 0; time_DCR = 0;
        BSS_R = 0; BSS_ERD = 0; BSS_DCR = 0;
        RSS_R = 0; RSS_ERD = 0; RSS_DCR = 0;
        readGraph();
        cout<<"Finish loading..."<<endl;  
    }
    
};

void Graph::readNM(){
    ifstream file;
    file.open(folder+graph_file+"attribute.txt");
    if(!file)
        handle_error("open");
    int line=0;
    string str;
    while(getline(file,str)){
        if(!line){
            n=atoi(str.substr(2,str.size()).c_str());
            line++;
        }
        else
        m=atoi(str.substr(2,str.size()).c_str());
    }
    file.close();
}

void Graph::readGraph(){
    ifstream file;
    file.open(folder+graph_file+"graph.txt");
    if(!file)
        handle_error("open");
    string str;
// int txtcnt=0;
    while(getline(file,str)){
        int size=str.size();
        int j=0,space=0;
        int start,end;
        double p;
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                if(!space){
                    start=atoi(str.substr(j,i-j).c_str());
                    space++;
                }
                else{
                    end=atoi(str.substr(j,i-j).c_str());
                    p=atof(str.substr(i,size-i).c_str());
                    break;
                }
                j=i+1;
            }
        }
        outdegree[start]++;
        indegree[end]++;
        degree[start]++;
        degree[end]++;
        outneighbors[start][end] = p;
        inneighbors[end][start] = p;
        Edge e1(start,end,p),e2(end,start,p);
        edges.push_back(e1);
        reverseEdges.push_back(e2);
// txtcnt++;
// if(txtcnt%1000==0) 
// cout<<start<<' '<<end<<' '<<p<<endl;
    }
    file.close();
}

void Graph::statistics(){
    for(int i=0;i<n;i++){
        if(degree[i]!=indegree[i]+outdegree[i])
        cout<<"error"<<endl;
    }
    cout<<"true"<<endl;
}

void Graph::degreeSort(int k=0, string str=""){
    vector<sort_v> sort_array(degree.size());
    for(int i=0; i<degree.size(); i++){
		sort_array[i].index=i;
        if(str=="indegree")
            sort_array[i].value=indegree[i];
        else if(str=="outdegree")
            sort_array[i].value=outdegree[i];
        else
            sort_array[i].value=degree[i];
	}
	sort(sort_array.begin(),sort_array.end(),compareVertex);
    for(int i=0; i<k; i++){
        cout<<sort_array[i].index<<' '<<sort_array[i].value<<endl;
    }
    for(auto i: sort_array){
        sortedVertices.push_back(i.index);
    }
}

void Graph::edgeDegree(){
    for(int i=0;i<edges.size();i++){
        Edge e=edges[i];
        edges[i].d=degree[e.start]+degree[e.end];
    }
}

void Graph::edgeSort(int k=0){
    edgeDegree();
    vector<Edge> sort_array=edges;
	sort(sort_array.begin(),sort_array.end(),compareEdge);
    sortedEdges=sort_array;
}

vector<Edge> Graph::randomEdge(){
    srand(time(0));
    vector<Edge> E=edges;
    for (int i = 0; i < E.size(); i++)
        E[i].d=rand()%m;
    sort(E.begin(),E.end(),compareEdge);
    return E;
}

void handle_error(const char* msg) 
{
	perror(msg);
	exit(255);
}