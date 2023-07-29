#include <deque>
#include <numeric>
#include <random>

// #include "graph.h"

#define MAX 1000000

//以vector的形式返回一个map里所有的key
vector<int> KeySet(map<int, double> test)
{
    vector<int> keys;
    for(map<int, double>::iterator it = test.begin(); it != test.end(); ++it)
        keys.push_back(it->first);
    return keys;
}
//找出以node为起点的所有h跳路径
vector<vector<int>> findPaths(Graph &g, int h, int node){
    vector<vector<int>> Path;
    if (h==0){
        vector<int> p;
        p.push_back(node);
        Path.push_back(p);
    }
    else{
        for (auto neighbor:KeySet(g.outneighbors[node])){
            for (auto path:findPaths(g, h-1, neighbor)){
                if (find(path.begin(),path.end(),node)==path.end()){
                    path.insert(path.begin(),node);
                    Path.push_back(path);
                }
            }
        }
        
    }
    return Path;
}

//构建一个vertex cover
vector<int> vertexCover(Graph &g, int h){
    vector<int> vc;
    if(h==0){
        vc=g.vertices;
        g.inVertexCover = vector<bool>(g.n, true);
    }
    else if(h==1){
        // Method 1: Traverse edges
        // g.edgeSort();
        // for(auto e:g.sortedEdges){
        //     if(find(vc.begin(),vc.end(),e.start)==vc.end() && find(vc.begin(),vc.end(),e.end)==vc.end()){
        //         if(g.degree[e.start]>=g.degree[e.end]){
        //             vc.push_back(e.start);
        //             g.inVertexCover[e.start]=true;
        //         }
        //         else{
        //             vc.push_back(e.end);
        //             g.inVertexCover[e.end]=true;
        //         }
        //         // vc.push_back(e.start);g.inVertexCover[e.start]=true;
        //         // vc.push_back(e.end);g.inVertexCover[e.end]=true;
        //     }
        // }
        // Method 2: Random selection && recursions
        vector<Edge> E = g.randomEdge();
        for(auto e:E){
            if(!g.inVertexCover[e.start] && !g.inVertexCover[e.end]){
                double randDouble = rand()%MAX/(double)(MAX+1);
                if(randDouble<=0.5){
                    vc.push_back(e.start);
                    g.inVertexCover[e.start]=true;
                }
                else{
                    vc.push_back(e.end);
                    g.inVertexCover[e.end]=true;
                }
                // if(g.degree[e.start]>=g.degree[e.end]){
                //     vc.push_back(e.start);
                //     g.inVertexCover[e.start]=true;
                // }
                // else{
                //     vc.push_back(e.end);
                //     g.inVertexCover[e.end]=true;
                // }
            }
        }
    }
    else{
        // 这部分首先把所有长为h的路径都存储起来
        vector<vector<int>> Path;
        for (int i = 0; i < g.n; i++){
            int node = i;
            //将生成的新vector插入到Path尾部
            vector<vector<int>> paths=findPaths(g, h, node);
            Path.insert(Path.end(),paths.begin(),paths.end());
        }
        // cout<<"Paths: "<<Path.size()<<endl;
        //随机访问一条边，直到所有边都放问到
        vector<int> v = vector<int>(Path.size(), -1);
        for(int i=0; i<Path.size(); i++) v[i] = i;
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(v.begin(),v.end(), default_random_engine(seed));
        bool flag = false;
        for (auto i:v){
            vector<int> p = Path[i];
            flag = false;
            int maxdeg = 0;
            for (auto node : p){
                if(g.inVertexCover[node]) flag = true;
                // maxdeg = g.degree[node] > maxdeg ? node : maxdeg;
            }
            if (!flag){
                maxdeg = rand()%h;
                g.inVertexCover[p[maxdeg]] = true;
                vc.push_back(p[maxdeg]);
            }
        }
        
    }
    cout<<"Vector cover size: "<<vc.size()<<endl;
    // for(auto e:g.edges){
    //     if(find(vc.begin(),vc.end(),e.start)==vc.end() && find(vc.begin(),vc.end(),e.end)==vc.end())
    //     cout<<"Edge "<<e.start<<", "<<e.end<<"has no vertices in VC, ERROR!!"<<endl;
    // }
    return vc;
}



// 构建索引（展示的是正向索引的构建过程；反向索引只需要对传入参数进行修改即可，flag变量用于标记是正向还是反向索引）
void constructIndex(bool flag, Graph &g, vector<map<int,vector<int>>> &reach, int K, int h, int dc){
    vector<map<int, double>> neighbors = vector<map<int, double>>(g.n, map<int, double>());
    vector<map<int, double>> rneighbors = vector<map<int, double>>(g.n, map<int, double>());
    if(flag) {neighbors = g.outneighbors; rneighbors = g.inneighbors;}
    else {neighbors = g.inneighbors; rneighbors = g.outneighbors;}
    int node;
    double r, erd;
    vector<vector<double>> dcr;
// int cnt = 0;
int dcr_size;
deque<int> q;
vector<bool> visited = vector<bool>(g.n, false);
vector<double> probs = vector<double>(g.n, 0);
vector<double> erds = vector<double>(g.n, 0);
vector<double> dcrs = vector<double>(g.n*h, 0);
vector<int> dis = vector<int>(g.n,0);
    for (int v = 0; v < g.n; v++){
        if(g.inVertexCover[v])dcr_size=dc; else dcr_size=h;
        dcr = vector<vector<double>>(dcr_size, vector<double>(2, 0));
        if (g.inVertexCover[v]){
            for(auto it:reach[v]){
                if(it.first==v)continue;
                for(int i = 0; i<dcr_size; i++){
                    dcr[i][0] = i+1;
                    dcr[i][1] = 0;
                }
                node = it.first;
                r = 1.0*it.second.size()/K;
                // int sum = accumulate(it.second.begin(), it.second.end(), 0);
                int sum = 0;
                // //计算拼接组合dcr
                if (node == v) continue;
                // cout<<"Node: "<<v<<" to "<<node<<", distances: ";
                for (auto d: it.second){
                    sum += d;
                    // cout<<d<<' ';
                    if (d>dcr_size) continue;
                    else dcr[d-1][1] += 1.0/K; 
                }
                // cout<<endl;
                erd = 1.0*sum/it.second.size();
                IndexEntry entry(node, r, erd, dcr);
                if (flag)
                    g.Lf[v].push_back(entry);
                else
                    g.Lb[v].push_back(entry);
            }
        }
        else{
            
            if (h == 1){
                for (auto it:neighbors[v]){
                    // cnt++;
                    node = it.first;
                    r = it.second;
                    erd = 1;
                    // dcr[node].push_back(1);
                    dcr = vector<vector<double>>(1, vector<double>(2, 0));
                    dcr[0][0] = 1;
                    dcr[0][1] = r;
                    IndexEntry entry(node, r, erd, dcr);
                    if (flag)
                        g.Lf[v].push_back(entry);
                    else
                        g.Lb[v].push_back(entry);
                }
            }
            else if(h==2){
                double prob;
                map<int, IndexEntry> entries;
                map<int, IndexEntry>::iterator i;
                for(auto it:neighbors[v]){
                    if(it.first==v) continue;
                    if(g.inVertexCover[it.first]){
                        node = it.first;
                        r = it.second;
                        erd = 1;
                        // dcr[node].push_back(1);
                        dcr = vector<vector<double>>(2, vector<double>(2, 0));
                        dcr[0][0] = 1;
                        dcr[0][1] = r;
                        dcr[1][0] = 2;
                        dcr[1][1] = 0;
                        // for(auto rit:rneighbors[node]){
                        //     if(neighbors[v].count(rit.first)){
                        //         dcr[1][1] = 1-(1-dcr[1][1])*(1-neighbors[v][rit.first]*rit.second);
                        //         erd = (1*it.second + 2*dcr[1][1])/ (it.second+dcr[1][1]);
                        //     }
                        // }
                        IndexEntry entry(node, r, erd, dcr);
                        entries.insert(pair<int, IndexEntry>(node, entry));
                        // entries[node] = entry;
                        // if (flag)
                        //     g.Lf[v].push_back(entry);
                        // else
                        //     g.Lb[v].push_back(entry);
                    }
                    else{
                        for(auto iit:neighbors[it.first]){
                            node = iit.first;
                            if(node==v)continue;
                            i = entries.find(node);
                            if(i==entries.end()){
                                erd = 2;
                                r = it.second*iit.second;
                                dcr = vector<vector<double>>(2, vector<double>(2, 0));
                                dcr[0][0] = 1;
                                dcr[0][1] = 0;
                                dcr[1][0] = 2;
                                dcr[1][1] = r;
                                IndexEntry entry(node, r, erd, dcr);
                                entries.insert(pair<int, IndexEntry>(node, entry));
                            }
                            else{
                                double tmp = (*i).second.reliability;
                                if(it.second*iit.second>tmp) (*i).second.reliability = it.second*iit.second;
                                // (*i).second.reliability = 1 - (1-(*i).second.reliability)*(1-it.second*iit.second);
                                if((*i).second.erd==2) (*i).second.erd = 2;
                                else{
                                    (*i).second.erd = ((*i).second.erd*tmp + 2*((*i).second.reliability-tmp)) / ((*i).second.reliability);
                                }
                                if(it.second*iit.second>dcr[1][1]) dcr[1][1] = it.second*iit.second;
                                // dcr[1][1] = 1 - (1-dcr[1][1])*(1-it.second*iit.second);
                            }
                        }
                    }
                }
                for(auto entry:entries){
                    if(entry.first==v||!g.inVertexCover[entry.first]) continue;
                    if (flag)
                        g.Lf[v].push_back(entry.second);
                    else
                        g.Lb[v].push_back(entry.second);
                }
            }
            else{
                q.clear();
                visited.clear();
                // probs = vector<double>(g.n, 0);
                fill(probs.begin(),probs.end(),0);
                // erds = vector<double>(g.n, -1);
                fill(erds.begin(),erds.end(),0);
                fill(dcrs.begin(),dcrs.end(),0);
                fill(dis.begin(),dis.end(),0);
                q.push_back(v);
                dis[v] = 0;
                probs[v] = 1;
                while(!q.empty()){
                    int node = q.front();
                    q.pop_front();
                    for(auto it:neighbors[node]){
                        if(visited[it.first]){
                            // if(g.inVertexCover[it.first])cout<<"More than once visiting vertex "<<it.first<<" from "<<v<<endl;
                            continue;
                        }
                        visited[it.first] = true;
                        dis[it.first] = dis[node]+1;
                        double p0 = probs[it.first];
                        probs[it.first] = 1-(1-probs[it.first])*(1-probs[node]*it.second);
                        erds[it.first] = (erds[it.first]*p0+dis[it.first]*(probs[it.first]-p0))/probs[it.first];
                        dcrs[(dis[it.first]-1)+h*it.first] = probs[it.first];
                        if(!g.inVertexCover[it.first]&&dis[it.first]<h){
                            q.push_back(it.first);
                        }
                    }
                }
                for(int node=0; node<g.n; node++){
                    if(node!=v&&g.inVertexCover[node]&&probs[node]!=0){
                        for(int i = 0; i<dcr_size; i++){
                            dcr[i][0] = i+1;
                            if(i==0) dcr[i][1] = dcrs[node*h+i];
                            else dcr[i][1] = dcrs[node*h+i]=dcr[i-1][1];
                        }
                        r = probs[node];
                        erd = erds[node];
                        IndexEntry entry(node, r, erd, dcr);
                        if (flag)
                            g.Lf[v].push_back(entry);
                        else
                            g.Lb[v].push_back(entry);
                    }
                }
            }
        }
    }
    // cout<<"Computation times: "<<cnt<<endl;
}



// Monte Calo模拟的过程
vector<map<int,vector<int>>> MonteCalo(bool flag, Graph &g, vector<int> vc, int K, int h, int b = MAX, double r = 0){
    srand(time(NULL));
    vector<map<int,vector<int>>> reach = vector<map<int,vector<int>>> (g.n, map<int,vector<int>>());
    vector<map<int, double>> neighbors = vector<map<int, double>>(g.n, map<int, double>());
    if(flag) neighbors = g.outneighbors;
    else neighbors = g.inneighbors;
    vector<bool> visited;
    deque<int> q;
    // int d = 0;
    // int cnt = 0;
    for(int i = 0; i < K; i++){
        for(auto v:vc){
        // for(int v=0; v<g.n; v++){
            map<int,int> reachable = map<int,int>();
            q.clear();
            q.push_back(v);
            reachable[v] = 0;
            reach[v][v].push_back(0);
            while (!q.empty()){
                int node = q.front();
                q.pop_front();
                for(auto it:neighbors[node]){
                    if (reachable.count(it.first))
                        continue;
                    double randDouble = rand()%MAX/(double)(MAX+1);
                    if (randDouble<=it.second){
                        q.push_back(it.first);
                        reachable[it.first]=reachable[node]+1;
                        if(g.inVertexCover[it.first]){
                            // q.push_back(it.first);
                            reach[v][it.first].push_back(reachable[it.first]);
                        }
                        // else{
                        //     if(reachable[it.first]<h){
                        //         q.push_back(it.first);
                        //         reach[v][it.first].push_back(reachable[it.first]);
                        //     }
                        // }
                    }  
                }
            }
            // cout<<"Finish sampling "<<v<<", sampled "<<++cnt<<" vertices."<<endl;
        }
        // if (i%50==0)
        // {
           cout<<"Sampling "<<i<<" times."<<endl;
        // }
    }
    // constructIndex(flag, g, reach, h);
    return reach;
}






// 一个输出函数，将index输出到文件中，用以查看其存储空间大小、索引条目等信息
void indexOutput(Graph g, int h){
    // ofstream   ofresult( "../output.txt ",ios::app); 
    ofstream ofoutput( g.folder+g.graph_file+to_string(h)+"_Index.txt",ios::ate); 
    for (int v = 0; v < g.n; v++){
        // if (g.inVertexCover[v]){
        for (auto entry: g.Lf[v]){
            ofoutput<<"Lf "<<v<<' '<<entry.node<<' '<<entry.reliability<<' '<<entry.erd<<' ';
            for(auto d: entry.dcr){
                ofoutput<<d[0]<<'~'<<d[1]<<' ';
            }
            ofoutput<<endl;
        }
        for (auto entry: g.Lb[v]){
            ofoutput<<"Lb "<<v<<' '<<entry.node<<' '<<entry.reliability<<' '<<entry.erd<<' ';
            for(auto d: entry.dcr){
                ofoutput<<d[0]<<'~'<<d[1]<<' ';
            }
            ofoutput<<endl;
        }
        // }
        // else{
        //     for (auto entry: g.Lf[v])
        //         ofoutput<<"Lf "<<v<<' '<<entry.node<<' '<<entry.reliability<<' '<<entry.erd<<' '<<entry.dcr[0][0]<<'~'<<entry.dcr[0][1]<<endl;
        //     for (auto entry: g.Lb[v])
        //         ofoutput<<"Lb "<<v<<' '<<entry.node<<' '<<entry.reliability<<' '<<entry.erd<<' '<<entry.dcr[0][0]<<'~'<<entry.dcr[0][1]<<endl;
        // }
        
    }
    ofoutput.close();
}



int cal(int a, int b = 10, int c = 20){
    return a+b+c;
}