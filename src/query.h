// #include "graph.h"
double countIndex(Graph &g){
    double r = 0;
    for(int i=0; i<g.n; i++){
        r += g.Lb[i].size();
        r += g.Lf[i].size();
    }
    return r/g.n;
}

int situation(Graph &g, int &s, int &t){
    clock_t t1, t2;
    t1 = clock();
    int ans = 0;
    if (g.inVertexCover[s]&&g.inVertexCover[t])
        return 1;
    else if (g.inVertexCover[s]&&(!g.inVertexCover[t]))
        return 2;
    else if ((!g.inVertexCover[s])&&g.inVertexCover[t])
        return 3;
    else
        return 4;
    t2 = clock();
    cout<<"Calculate status: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;
    return ans;
}

double min_d(double a, double b){
    if(a<b) return a;
    else return b;
}

double max_d(double a, double b){
    if(a>b) return a;
    else return b;
}

double reliabilityQuery(Graph &g, int &s, int &t){
    if(s==t) return 1;
    if(g.Lf[s].size()==0||g.Lb[t].size()==0) return 0;
    int status = situation(g, s, t);
    // clock_t t1, t2;
    // t1 = clock();
    double r = 0;
    int pos = 0;
    bool flag = false;
    if(status == 1){
        for(auto lf_s: g.Lf[s]){
            if(lf_s.node == t){
                flag = true;
                r = max(r, lf_s.reliability);
            }
        }
        for(auto lb_t: g.Lb[t]){
            if(lb_t.node == s){
                flag = true;
                r = max(r, lb_t.reliability);
            }
        }
        if(!flag){
            for(auto lf_s: g.Lf[s]){
                for (;pos<g.Lb[t].size();pos++){
                    IndexEntry lb_t = g.Lb[t][pos];
                    if (lb_t.node > lf_s.node) break;
                    if(lb_t.node == lf_s.node){
                        // if (lf_s.reliability*lb_t.reliability>r)
                        //     cout<<"New reliability: through "<<lb_t.node<<" is "<<lf_s.reliability*lb_t.reliability<<" > "<<r<<endl;
                        r = max(r, lf_s.reliability*lb_t.reliability);
                        break;
                    }
                }
            }
        }
    }
    else if (status == 2){
        for(auto lb_t: g.Lb[t]){
            int v = lb_t.node;
            if(situation(g, s, v)!=1) continue;
            r = max(r, reliabilityQuery(g, s, v)*lb_t.reliability);
        }
    }
    else if (status == 3){
        for(auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            if(situation(g, u, t)!=1) continue;
            r = max(r, reliabilityQuery(g, u, t)*lf_s.reliability);
        }
    }
    else if (status == 4){
        for (auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            for(auto lb_t: g.Lb[t]){
                int v = lb_t.node;
                if(situation(g, u, v)!=1) continue;
                r = max(r, lf_s.reliability*reliabilityQuery(g, u, v)*lb_t.reliability);
            }
        }
    }
    else{
        cout<<"Status Error"<<endl;
        handle_error(0);
    }
    return r;
}

double ERD_min(Graph &g, int &s, int &t){
    if(s==t) return 0;
    if(g.Lf[s].size()==0||g.Lb[t].size()==0) return 0;
    int status = situation(g, s, t);
    int pos = 0;
    bool flag = false;
    // double r = reliabilityQuery(g, s, t);
    double e = MAX;
    if(status == 1){
        for(auto lf_s: g.Lf[s]){
            if(lf_s.node == t){
                flag = true;
                e = min_d(lf_s.erd,e);
            }
        }
        for(auto lb_t: g.Lb[t]){
            if(lb_t.node == s){
                flag = true;
                e = min_d(lb_t.erd,e);
            }
        }
        if(!flag){
            for(auto lf_s: g.Lf[s]){
                for (;pos<g.Lb[t].size();pos++){
                    IndexEntry lb_t = g.Lb[t][pos];
                    if (lb_t.node > lf_s.node) break;
                    if(lb_t.node == lf_s.node){
                        // if (lf_s.reliability*lb_t.reliability>r)
                        //     cout<<"New reliability: through "<<lb_t.node<<" is "<<lf_s.reliability*lb_t.reliability<<" > "<<r<<endl;
                        e = min_d(e, lf_s.erd+lb_t.erd);
                        break;
                    }
                }
            }
        }
    }
    else if (status == 2){
        for(auto lb_t: g.Lb[t]){
            int v = lb_t.node;
            if(situation(g, s, v)!=1) continue;
            e = min_d(lb_t.erd+ERD_min(g, s, v), e);
        }
    }
    else if (status == 3){
        for(auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            if(situation(g, u, t)!=1) continue;
            e = min_d(lf_s.erd+ERD_min(g, u, t), e);
        }
    }
    else if (status == 4){
        for (auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            for(auto lb_t: g.Lb[t]){
                int v = lb_t.node;
                // erd = min(erd, ERDQuery(g, u, v)+2);
                if(situation(g, u, v)!=1) continue;
                e = min_d(lf_s.erd+lb_t.erd+ERD_min(g, u, v), e);
            }
        }
    }
    else{
        cout<<"Status Error"<<endl;
        handle_error(0);
    }
    if(e==MAX) e=0;
    return e;
}

double ERD_max(Graph &g, int &s, int &t){
    if(s==t) return 0;
    if(g.Lf[s].size()==0||g.Lb[t].size()==0) return 0;
    int status = situation(g, s, t);
    int pos = 0;
    bool flag = false;
    // double numerator = 0, denominator = 1;
    // double r = reliabilityQuery(g, s, t);
    double e = 0;
    if(status == 1){
        for(auto lf_s: g.Lf[s]){
            if(lf_s.node == t){
                flag = true;
                e = max_d(lf_s.erd,e);
            }
        }
        for(auto lb_t: g.Lb[t]){
            if(lb_t.node == s){
                flag = true;
                e = max_d(lb_t.erd,e);
            }
        }
        if(!flag){
            for(auto lf_s: g.Lf[s]){
                for (;pos<g.Lb[t].size();pos++){
                    IndexEntry lb_t = g.Lb[t][pos];
                    if (lb_t.node > lf_s.node) break;
                    if(lb_t.node == lf_s.node){
                        // if (lf_s.reliability*lb_t.reliability>r)
                        //     cout<<"New reliability: through "<<lb_t.node<<" is "<<lf_s.reliability*lb_t.reliability<<" > "<<r<<endl;
                        e = max_d(e, lf_s.erd+lb_t.erd);
                        break;
                    }
                }
            }
        }
    }
    else if (status == 2){
        for(auto lb_t: g.Lb[t]){
            int v = lb_t.node;
            if(situation(g, s, v)!=1) continue;
            // numerator  += (ERDQuery(g, s, v)[1]+1) * (r/lb_t.reliability);
            // numerator  += (ERD_max(g, s, v)+lb_t.erd) * (lb_t.reliability*reliabilityQuery(g, s, v));
            e = max_d(ERD_max(g, s, v) + lb_t.erd, e);
        }
        // denominator = r;
        // if(r==0)e=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e = numerator / denominator;
    }
    else if (status == 3){
        for(auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            if(situation(g, u, t)!=1) continue;
            // numerator  += (ERDQuery(g, u, t)[1]+1) * (r/lf_s.reliability);
            // numerator  += (ERD_max(g, u, t)+lf_s.erd) * (lf_s.reliability*reliabilityQuery(g, u, t));
            // denominator *= 1 - lf_s.reliability*reliabilityQuery(g, u, t);
            e = max_d(ERD_max(g, u, t) + lf_s.erd, e);
        }
        // denominator = r;
        // if(r==0)e=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e = numerator / denominator;
    }
    else if (status == 4){
        for (auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            for(auto lb_t: g.Lb[t]){
                int v = lb_t.node;
                if(situation(g, u, v)!=1) continue;
                e = max_d(ERD_max(g, u, v)+lf_s.erd+lb_t.erd, e);
                // numerator  += (ERDQuery(g, u, v)[1]+2) * (r/lf_s.reliability/lb_t.reliability);
                // numerator  += (ERD_max(g, u, v)+lf_s.erd+lb_t.erd) * (lf_s.reliability*reliabilityQuery(g, u, v)*lb_t.reliability);
                // denominator *= 1 - lf_s.reliability*reliabilityQuery(g, u, v)*lb_t.reliability;
            }
        }
        // denominator = r;
        // if(r==0)e=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e = numerator / denominator;
    }
    else{
        cout<<"Status Error"<<endl;
        handle_error(0);
    }
    return e;
}

vector<double> ERDQuery(Graph &g, int &s, int &t){
    if(s==t) return {0, 0};
    if(g.Lf[s].size()==0||g.Lb[t].size()==0) return {0, 0};
    int status = situation(g, s, t);
    int pos = 0;
    bool flag = false;
    double numerator = 0, denominator = 1;
    // double r = reliabilityQuery(g, s, t);
    vector<double> ERD = {MAX, 0};
    double e0 = MAX, e1 = 0;
    if(status == 1){
        for(auto lf_s: g.Lf[s]){
            if(lf_s.node == t){
                flag = true;
                e0 = min_d(lf_s.erd,e0);
                e1 = max_d(lf_s.erd,e1);
            }
        }
        for(auto lb_t: g.Lb[t]){
            if(lb_t.node == s){
                flag = true;
                e0 = min_d(lb_t.erd,e0);
                e1 = max_d(lb_t.erd,e1);
            }
        }
        if(!flag){
            for(auto lf_s: g.Lf[s]){
                for (;pos<g.Lb[t].size();pos++){
                    IndexEntry lb_t = g.Lb[t][pos];
                    if (lb_t.node > lf_s.node) break;
                    if(lb_t.node == lf_s.node){
                        // if (lf_s.reliability*lb_t.reliability>r)
                        //     cout<<"New reliability: through "<<lb_t.node<<" is "<<lf_s.reliability*lb_t.reliability<<" > "<<r<<endl;
                        e0 = min_d(e0, lf_s.erd+lb_t.erd);
                        e1 = max_d(e1, lf_s.erd+lb_t.erd);
                        break;
                    }
                }
            }
        }
    }
    else if (status == 2){
        for(auto lb_t: g.Lb[t]){
            int v = lb_t.node;
            if(situation(g, s, v)!=1) continue;
            e0 = min_d(lb_t.erd+ERD_min(g, s, v), e0);
            e1 = max_d(lb_t.erd+ERD_max(g, s, v), e1);
            // numerator  += (ERDQuery(g, s, v)[1]+1) * (r/lb_t.reliability);
            // numerator  += (ERD_max(g, s, v)+lb_t.erd) * (lb_t.reliability*reliabilityQuery(g, s, v));
        }
        // denominator = r;
        // if(r==0)e1=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e1 = numerator / denominator;
    }
    else if (status == 3){
        for(auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            if(situation(g, u, t)!=1) continue;
            e0 = min_d(lf_s.erd+ERD_min(g, u, t), e0);
            e1 = max_d(lf_s.erd+ERD_max(g, u, t), e1);
            // numerator  += (ERDQuery(g, u, t)[1]+1) * (r/lf_s.reliability);
            // numerator  += (ERD_max(g, u, t)+lf_s.erd) * (lf_s.reliability*reliabilityQuery(g, u, t));
            // denominator *= 1 - lf_s.reliability*reliabilityQuery(g, u, t);
        }
        // denominator = r;
        // if(r==0)e1=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e1 = numerator / denominator;
    }
    else if (status == 4){
        for (auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            for(auto lb_t: g.Lb[t]){
                int v = lb_t.node;
                if(situation(g, u, v)!=1) continue;
                // erd = min(erd, ERDQuery(g, u, v)+2);
                e0 = min_d(lf_s.erd+lb_t.erd+ERD_min(g, u, v), e0);
                e1 = max_d(lf_s.erd+lb_t.erd+ERD_max(g, u, v), e1);
                // numerator  += (ERDQuery(g, u, v)[1]+2) * (r/lf_s.reliability/lb_t.reliability);
                // numerator  += (ERD_max(g, u, v)+lf_s.erd+lb_t.erd) * (lf_s.reliability*reliabilityQuery(g, u, v)*lb_t.reliability);
                // denominator *= 1 - lf_s.reliability*reliabilityQuery(g, u, v)*lb_t.reliability;
            }
        }
        // denominator = r;
        // if(r==0)e1=0;
        // denominator = reliabilityQuery(g, s, t);
        // else e1 = numerator / denominator;
    }
    else{
        cout<<"Status Error"<<endl;
        handle_error(0);
    }
    // cout<<"Finish ERD query."<<endl;
    if(e0==MAX) e0=0;
    ERD[0] = e0;
    ERD[1] = e1;
    return ERD;
}

double dcrPlus(int num, vector<vector<double>> &v){
    double add=0;
    for(int i=0; i<num; i++){
        add += v[i][1];
    }
    return add;
}

double DCRQuery(Graph &g, int &s, int &t, int d, int h){
    if(s==t) return 1;
    if(g.Lf[s].size()==0||g.Lb[t].size()==0) return 0;
    int status = situation(g, s, t);
    int pos = 0;
    double dcr = 0;
    double tmp = 0;
    bool flag = false;
    if(status == 1){
        for(auto lf_s: g.Lf[s]){
            if(lf_s.node == t){
                flag = true;
                // double tmp = 0;
                // for (int i = 0; i < lf_s.dcr.size(); i++)
                //     tmp+=lf_s.dcr[i][1];
                int m = lf_s.dcr.size();
                tmp = dcrPlus(m,lf_s.dcr);
                dcr = max(dcr, tmp);
            }
        }
        for(auto lb_t: g.Lb[t]){
            if(lb_t.node == s){
                flag = true;
                // double tmp = 0;
                // for (int i = 0; i < lb_t.dcr.size(); i++)
                //     tmp += lb_t.dcr[i][1];
                int m = lb_t.dcr.size();
                tmp = dcrPlus(m,lb_t.dcr);
                dcr = max(dcr, tmp);
            }
        }
        if(!flag){
            for(auto lf_s: g.Lf[s]){
                for (;pos<g.Lb[t].size();pos++){
                    IndexEntry lb_t = g.Lb[t][pos];
                    if (lb_t.node > lf_s.node) break;
                    if(lb_t.node == lf_s.node){
                        for(int i=2; i<d; i++){
                            int m = d-i;
                            dcr = max(dcr, dcrPlus(i, lf_s.dcr)*dcrPlus(m, lb_t.dcr));
                            // dcr = max(dcr, DCRQuery(g, s, lf_s.node, i)*DCRQuery(g, lb_t.node, t, d-i));
                        }
                        // dcr = max(dcr, DCRQuery(g, s, lf_s.node, i)*DCRQuery(g, lb_t.node, t, d-i));
                        break;
                    }
                }
            }
        }
    }
    else if (status == 2){
        for(auto lb_t: g.Lb[t]){
            int v = lb_t.node;
            if(situation(g,s,v)!=1)continue;
            // for(int i=0;i<h;i++){
            //     dcr = max(dcr, DCRQuery(g, s, v, d-i, h)*dcrPlus(i, lb_t.dcr));
            // }
            // dcr = max(dcr, DCRQuery(g, s, v, d-1, h)*lb_t.reliability);
            dcr = max(dcr, DCRQuery(g, s, v, d-1, h)*lb_t.reliability);
            // dcr = max(dcr, DCRQuery(g, s, v, d-2, h)*dcrPlus(2,lb_t.dcr));
        }
    }
    else if (status == 3){
        for(auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            if(situation(g,u,t)!=1)continue;
            // for(int i=0;i<h;i++){
            //     dcr = max(dcr, DCRQuery(g, u, t, d-i, h)*dcrPlus(i, lf_s.dcr));
            // }
            dcr = max(dcr, DCRQuery(g, u, t, d-1, h)*lf_s.reliability);
            // dcr = max(dcr, DCRQuery(g, u, t, d-2, h)*dcrPlus(2, lf_s.dcr));
        }
    }
    else if (status == 4){
        for (auto lf_s: g.Lf[s]){
            int u = lf_s.node;
            for(auto lb_t: g.Lb[t]){
                int v = lb_t.node;
                if(situation(g,u,v)!=1)continue;
                // for(int i=0;i<h;i++){
                //     for(int j=0;j<h;j++){
                //         dcr = max(dcr, dcrPlus(i,lf_s.dcr)*DCRQuery(g, u, v, d-i-j, h)*dcrPlus(j,lb_t.dcr));
                //     }
                // }
                dcr = max(dcr, lf_s.reliability*DCRQuery(g, u, v, d-2, h)*lb_t.reliability);
            }
        }
    }
    else{
        cout<<"Status Error"<<endl;
        handle_error(0);
    }
    // cout<<"Finish DCR query."<<endl;
    return dcr;
}



vector<double> query(Graph &g, int s, int t, int d, int h){
    vector<double> res = vector<double>(4,0);
    clock_t t1, t2;
    double r, c;
    vector<double> e;
    if (s>=g.n || t>=g.n)
        handle_error("Vertex does not exist.");
    t1 = clock();
    r = reliabilityQuery(g, s, t);
    t2 = clock();
    g.time_R += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Reliability query "<<s<<' '<<t<<": "<<r<<endl;

    t1 = clock();
    e = ERDQuery(g, s, t);
    // double e1 = ERD_min(g,s,t);
    // double e2 = ERD_max(g,s,t);
    t2 = clock();
    g.time_ERD += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Expected Reliable Distance query "<<s<<' '<<t<<": "<<e[0]<<','<<e[1]<<endl;

    t1 = clock();
    c = DCRQuery(g, s, t, d, h);
    t2 = clock();
    g.time_DCR += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Distance Constrained Reliability query "<<s<<' '<<t<<": "<<c<<endl;
    // cout<<"Reliability: "<<r<<"\nExpected Reliable Distance: "<<e[0]<<", "<<e[1]<<"\nDistance Constrained Reliability: "<<c<<endl;
    res[0] = r;
    res[1] = e[0];
    res[2] = e[1];
    res[3] = c;
    return res;
}

void indexQuery(Graph &g, int d, int h){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+to_string(h)+"_Index_results0.txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    // clock_t t0, t1;
    // t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = query(g, s, t, d, h);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<' '<<res[3]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        // cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD_min: "<<res[1]<<", ERD_max: "<<res[2]<<", DCR: "<<res[3]<<endl;
        // if(cnt==100)break;
    }
    // t1 = clock();
    cout<<"Average Index query time:"<<(g.time_R+g.time_ERD+g.time_DCR)/cnt<<"s"<<endl;
    cout<<"Reliability query: "<<g.time_R/cnt<<"s"<<endl;
    cout<<"ERD query: "<<g.time_ERD/cnt<<"s"<<endl;
    cout<<"DCR query: "<<g.time_DCR/cnt<<"s"<<endl;
    file.close();
    of.close();
}

//一对query进行Naive Monte Calo的过程
vector<double> nmc_query(Graph &g, int s, int t, int d, int size){
    vector<double> res = vector<double>(3, 0);
    deque<int> q;
    vector<int> reach;
    vector<bool> visited = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, -1);
    bool flag = false;
    // srand(time(0));
    for(int i = 0; i<size; i++){
        fill(dis.begin(), dis.end(), -1);
        fill(visited.begin(),visited.end(),false);
        flag = false;
        q.clear();
        q.push_back(s);
        dis[s] = 0;
        while(!q.empty()){
            int node = q.front();
            q.pop_front();
            for(auto it:g.outneighbors[node]){
                if(visited[it.first]||it.second==0)continue;
                double randDouble = rand()%MAX/(double)(MAX+1);
                if(randDouble<=it.second){
                    dis[it.first] = dis[node]+1;
                    if(it.first==t){
                        flag = true;
                        break;
                    }
                    q.push_back(it.first);
                    visited[it.first] = true;
                }
            }
            if(flag){
                reach.push_back(dis[t]);
                break;
            }
        }
    }
    //res[0] = reliability; res[1] = ERD; res[2] = DCR(d)
    res[0] = 1.0*reach.size()/size;
    int sum = 0, cnt = 0;
    for(int i = 0; i<reach.size(); i++){
        if(reach[i]<=d) cnt++;
        sum += reach[i];
    }
    if(!reach.size()) res[1] = 0;
    else res[1] = 1.0*sum/reach.size();
    res[2] = 1.0*cnt/size;
    // cout<<"Reach size: "<<reach.size()<<endl;
    // for(int i = 0; i < reach.size(); i++) cout<<reach[i]<<' ';
    // cout<<endl;
    return res;
}

void nmc(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"NMC_results_"+to_string(K)+".txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    clock_t t0, t1;
    t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = nmc_query(g, s, t, d, K);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
    }
    t1 = clock();
    cout<<"Average NMC time: "<<(double)(t1-t0)/CLOCKS_PER_SEC/cnt<<"s"<<endl;
    file.close();
    of.close();
}

void readIndex(Graph &g, int h, int dc){
    ifstream file;
    file.open(g.folder+g.graph_file+to_string(h)+"_Index.txt");
    if(!file)
        handle_error("open");
    int line=0;
    string str;
    int src, dest;
    double prob, erd;
    vector<vector<double>> dcr_in = vector<vector<double>>(dc, vector<double>(2, 0));
    vector<vector<double>> dcr_out = vector<vector<double>>(h, vector<double>(2, 0));
    vector<vector<double>> dcr;
    while(getline(file,str)){
        // cout<<to_string(line)+str[str.size()-1];
        vector<string> sub = vector<string>();
        int j = 0;
        for(int i=0; i<str.size(); i++){
            if(str[i]==' '){
                // cout<<j<<' '<<i<<endl;
                sub.push_back(str.substr(j, i-j));
                j = i+1;
            }
        }
        src = atoi(sub[1].c_str());
        dest = atoi(sub[2].c_str());
        prob = atof(sub[3].c_str());
        erd = atof(sub[4].c_str());
        if(g.inVertexCover[src]) dcr = dcr_in;
        else dcr = dcr_out;

        if(sub.size()-5==dc){
            g.inVertexCover[src] = true;
            dcr = dcr_in;
        }
        else if(sub.size()-5==h){
            dcr = dcr_out;
        }
        else handle_error("Index");
        for(int i=5; i<sub.size(); i++){
            for(int j=0; j<sub[i].size(); j++){
                if(sub[i][j]=='~'){
                    dcr[i-5][0]=atof(sub[i].substr(0,j).c_str());
                    dcr[i-5][1]=atof(sub[i].substr(j+1).c_str());
                    break;
                }
            }
        }
        IndexEntry entry(dest, prob, erd, dcr);
        if(sub[0]=="Lf"){
            g.Lf[src].push_back(entry);
        }
        else if(sub[0]=="Lb"){
            g.Lb[src].push_back(entry);
        }
        else
            handle_error("Index");
        line++;
        if(line%1000000==0){
            cout<<"Finish loading indexes "<<line<<endl;
        }
    }
    file.close();
    cout<<"Finish reading index..."<<endl;
}

int random_walk(Graph &g, int s, int d){
    int t = -1;
    vector<int> q;
    q.push_back(s);
    while (d--)
    {
        int tmp = rand() % g.outneighbors[q.back()].size();
        if (!tmp) break;
        int cnt = 0;
        for (auto it:g.outneighbors[q.back()]){
            cnt++;
            if(cnt==tmp){
                if (find(q.begin(),q.end(),it.first)!=q.end()){
                    d++;
                    break;
                }
                q.push_back(it.first);
            }
        }
    }
    if(!(d+1)) t = q.back();
    return t;
}

void generate_queries(Graph &g, int size){
    ofstream of(g.folder+g.graph_file+"queries.txt",ios::ate); 
    int cnt=0;
    while (size--){
        // double temp = rand() % MAX/(double)(MAX+1);
        int s = rand() % g.n;
        int d = rand() % 5 + 2;
        int t = random_walk(g, s, d);
        if(t==-1||t==s){
            size++;
            continue;
        }
        else{
            of<<s<<' '<<t<<endl;
            if(size%10==0)
                cout<<"Generate "<<(++cnt)*10<<" queries."<<endl;
        }
    }
    
    of.close();
}


// BSS方法查询query
vector<double> BSS_query(Graph &g, int s, int t, int d, int size){
    vector<double> res = vector<double>(3, 0);
    map<int, double> m = g.outneighbors[s];
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0) return res;
    double omega, Omega=0;
    int newsize;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[s][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[s][it.first] = 1;
                break;
            }
            j++;
        }
        newsize = ceil(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        vector<double> tmpres = nmc_query(g0, s, t, d, newsize);
        for(int k=0; k<3; k++){
            if(k!=1)
                res[k]+=tmpres[k]*omega;
            else{
                if(tmpres[k]==0) Omega -= omega;
                else res[k]+=tmpres[k]*omega;
            }
            // cout<<tmpres[k]<<' '<<omega<<endl;
            // if(k==1)cout<<"Partial ERD: "<<tmpres[k]<<", Stratum Prob: "<<omega<<endl;
        }
        // cout<<endl;
    }
    // cout<<"Numerator: "<<res[1]<<", denominator: "<<Omega<<endl;
    if(Omega!=0)res[1] = res[1]/Omega;
    // for(auto v:res)cout<<v<<' ';
    // cout<<endl;
    return res;
}

void BSS(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"BSS_results_"+to_string(K)+".txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    clock_t t0, t1;
    t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = BSS_query(g, s, t, d, K);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==10) break;
    }
    t1 = clock();
    cout<<"Average BSS time: "<<(double)(t1-t0)/CLOCKS_PER_SEC/cnt<<"s"<<endl;
    file.close();
    of.close();
}


// RSS方法查询query
vector<double> RSS_query(Graph &g, int s, int t, int d, int size, int b){
    static int threshold = 10;
    if(size<=threshold||g.RSSvisit[b]) return nmc_query(g, s, t, d, size);
    vector<double> res = vector<double>(3, 0);
    map<int, double> m = g.outneighbors[b];
    // map<int, double> m = g.inneighbors[t];
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0) return res;
    double omega, Omega=0;
    int newsize, base;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[b][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[b][it.first] = 0;
                base = it.first;
                break;
            }
            j++;
        }
        newsize = ceil(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        // cout<<newsize<<' '<<omega<<endl;
        // vector<double> tmpres = nmc_query(g0, s, t, d, newsize);
        vector<double> tmpres = RSS_query(g0, s, t, d, newsize, base);
        for(int k=0; k<3; k++){
            if(k!=1)
                res[k]+=tmpres[k]*omega;
            else{
                if(tmpres[k]==0) Omega -= omega;
                else res[k]+=tmpres[k]*omega;
            }
            // cout<<tmpres[k]<<' ';
            // cout<<tmpres[k]<<' '<<omega<<endl;
            // if(k==1)cout<<"Partial ERD: "<<tmpres[k]<<", Stratum Prob: "<<omega<<endl;
        }
        // cout<<endl;
    }
    // cout<<"Numerator: "<<res[1]<<", denominator: "<<Omega<<endl;
    if(Omega!=0)res[1] = res[1]/Omega;
    // for(auto v:res)cout<<v<<' ';
    // cout<<endl;
    return res;
}

void RSS(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"RSS_results_"+to_string(K)+".txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    clock_t t0, t1;
    t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = RSS_query(g, s, t, d, K, s);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>0)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==10)break;
    }
    t1 = clock();
    cout<<"Average RSS time: "<<(double)(t1-t0)/CLOCKS_PER_SEC/cnt<<"s"<<endl;
    file.close();
    of.close();
}



// LPS的准备工作：生成一个几何采样的值
int Geo(double p){
    // srand(unsigned(time(NULL)));
    int num = 0;
    double randDouble = rand()%MAX/(double)(MAX+1);
    double sum = 0;
    while (true)
    {
        num++;
        sum += 1.0*pow(1-p, num-1)*p;
        // cout<<sum<<' ';
        if (randDouble<sum)
        break;
    }
    num -= 1;
    // num = rand()%(int)(2);
    return num;
}

bool compareHv(vector<double> a, vector<double> b){
	return (int)a[1]<(int)b[1];
}


// LPS方法查询query
vector<double> LPS_query(Graph &g, int s, int t, int d, int size){
    vector<double> res = vector<double>(3, 0);
    deque<int> q;
    vector<int> cv = vector<int>(g.n, -1);
    vector<deque<vector<double>>> hv = vector<deque<vector<double>>>(g.n ,deque<vector<double>>());
    vector<bool> visited = vector<bool>(g.n, false);
    vector<bool> init = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, -1);
    vector<int> reach;
    bool flag = false;
    if(s==t){
        res[0] = 1;
        res[1] = 0;
        res[2] = 1;
        return res;
    }
    for(int i = 0; i < size; i++){
        fill(visited.begin(),visited.end(),false);
        fill(dis.begin(), dis.end(), -1);
        q.clear();
        q.push_back(s);
        visited[s] = true;
        flag = false;
        dis[s] = 0;
        while (!q.empty()){
            int node = q.front();
            q.pop_front();
            if (!init[node]){
                // cout<<"Initiating node "<<node<<endl;
                cv[node] = 0;
                hv[node].clear();
                for(auto nbr: g.outneighbors[node]){
                    int geo = Geo(nbr.second);
                    hv[node].push_back({(double)nbr.first, (double)(geo+cv[node]), nbr.second});
                    // cout<<"Create new pair "<<nbr.first<<' '<<geo<<endl;
                }
                init[node] = true;
                sort(hv[node].begin(), hv[node].end(), compareHv);
            }
            // cout<<"Now visiting node "<<node<<" in round "<<i<<", its cv now is: "<<cv[node]<<", out_degree is "<<g.outdegree[node]<<endl;
            // cout<<"Its neighbors(hv) now is: "<<endl;
            // for(auto it:hv[node]){
            //     cout<<it[0]<<' '<<it[1]<<endl;
            // }
            if(g.outdegree[node]==0) {cv[node] += 1; continue;}
            while ((int)hv[node].front()[1]==cv[node]){
                vector<double> it = hv[node].front();
                hv[node].pop_front();
                // cout<<"Pop an Hv pair: "<<it[0]<<' '<<it[1]<<endl;
                if(!visited[it[0]]){
                    dis[(int)it[0]] = dis[node]+1;
                    q.push_back((int)it[0]);
                    visited[(int)it[0]] = true;
                    // cout<<"Round "<<i<<", visiting "<<node<<", push "<<it[0]<<" in hv, its cv: "<<cv[it[0]]<<endl;
                    // cout<<"Push "<<it[0]<<" in hv, its cv: "<<cv[it[0]]<<", Round "<<i<<endl;
                }
                int newX = Geo(it[2]);
                hv[node].push_back({it[0], (double)(newX+cv[node]+1), it[2]});
                // cout<<"Round "<<i<<", change "<<node<<" of neighbor "<<it[0]<<" from "<<it[1]<<" to "<<(double)(newX+cv[node]+1)<<endl;
                sort(hv[node].begin(), hv[node].end(), compareHv);
                if(it[0]==t){
                    flag = true;
                    reach.push_back(dis[t]);
                    // cout<<"Found on round "<<i<<endl;
                    // break;
                }
            }
            cv[node] += 1;
            if(flag) break;
            // cout<<"Node is "<<node<<", now cv is "<<cv[node]<<endl;
        }
        
    }
    res[0] = 1.0*reach.size()/size;
    int sum = 0, cnt = 0;
    for(int i = 0; i<reach.size(); i++){
        if(reach[i]<=d) cnt++;
        sum += reach[i];
    }
    if(!reach.size()) res[1] = 0;
    else res[1] = 1.0*sum/reach.size();
    res[2] = 1.0*cnt/size;
    // for(auto it:res)
    //     cout<<it<<' ';
    // cout<<endl;
    // cout<<endl;
    // for(int i=0; i<g.n; i++){
    //     if(cv[i]>50)
    //         cout<<i<<' '<<cv[i]<<endl;
    // }
    return res;
}

void LPS(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    // ofstream of(g.folder+g.graph_file+"LP+_results_"+to_string(K)+".txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    clock_t t0, t1;
    t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = LPS_query(g, s, t, d, K);
        // of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>0)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==100)break;
    }
    t1 = clock();
    cout<<"Average LP+ time: "<<(double)(t1-t0)/CLOCKS_PER_SEC/cnt<<"s"<<endl;
    file.close();
    // of.close();
}























double nmcR(Graph &g, int &s, int &t, int size){
    double r;
    if(s==t) return 1;
    int cnt = 0;
    deque<int> q;
    vector<bool> visited = vector<bool>(g.n, false);
    // srand(time(0));
    for(int i = 0; i<size; i++){
        fill(visited.begin(),visited.end(),false);
        q.clear();
        q.push_back(s);
        visited[s] = true;
        while(!q.empty()){
            int node = q.front();
            q.pop_front();
            for(auto it:g.outneighbors[node]){
                if(visited[it.first]||it.second==0)continue;
                double randDouble = rand()%MAX/(double)(MAX+1);
                if(randDouble<=it.second){
                    if(it.first==t){
                        cnt++;
                        break;
                    }
                    q.push_back(it.first);
                    visited[it.first] = true;
                }
            }
        }
    }
    r = 1.0*cnt/size;
    return r;
}

double nmcERD(Graph &g, int &s, int &t, int size){
    double erd;
    if(s==t) return 0;
    int d = 0, cnt = 0;
    deque<int> q;
    vector<bool> visited = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, 0);
    // srand(time(0));
    for(int i = 0; i<size; i++){
        fill(visited.begin(),visited.end(),false);
        fill(dis.begin(), dis.end(), -1);
        q.clear();
        q.push_back(s);
        dis[s] = 0;
        visited[s] = true;
        while(!q.empty()){
            int node = q.front();
            q.pop_front();
            for(auto it:g.outneighbors[node]){
                if(visited[it.first]||it.second==0)continue;
                double randDouble = rand()%MAX/(double)(MAX+1);
                if(randDouble<=it.second){
                    dis[it.first] = dis[node]+1;
                    if(it.first==t){
                        d += dis[it.first];
                        cnt++;
                        break;
                    }
                    q.push_back(it.first);
                    visited[it.first] = true;
                }
            }
        }
    }
    if(cnt==0) return 0;
    erd = 1.0*d/cnt;
    return erd;
}

double nmcDCR(Graph &g, int &s, int &t, int dc, int size){
    double dcr;
    if(s==t) return 1;
    int cnt = 0;
    deque<int> q;
    vector<bool> visited = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, 0);
    // srand(time(0));
    for(int i = 0; i<size; i++){
        fill(visited.begin(),visited.end(),false);
        fill(dis.begin(), dis.end(), -1);
        q.clear();
        q.push_back(s);
        dis[s] = 0;
        visited[s] = true;
        while(!q.empty()){
            int node = q.front();
            q.pop_front();
            for(auto it:g.outneighbors[node]){
                if(visited[it.first]||it.second==0)continue;
                double randDouble = rand()%MAX/(double)(MAX+1);
                if(randDouble<=it.second){
                    dis[it.first] = dis[node]+1;
                    if(it.first==t){
                        if(dis[it.first]<=dc)cnt++;
                        break;
                    }
                    q.push_back(it.first);
                    visited[it.first] = true;
                }
            }
        }
    }
    dcr = 1.0*cnt/size;
    return dcr;
}

double BSSR(Graph &g, int &s, int &t, int size){
    double r = 0;
    map<int, double> m = g.outneighbors[s];
    int num = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||num==0) return r;
    if(s==t) return 1;
    double omega, Omega=0;
    int newsize;
    Graph g0 = g;
    // cout<<s<<' '<<t<<' '<<size<<endl;
    for(int i=0; i<num; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[s][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[s][it.first] = 1;
                break;
            }
            j++;
        }
        newsize = floor(size * omega);
        // cout<<newsize<<' '<<size<<' ';
        if(newsize==0) continue;
        Omega += omega;
        double tmpr = nmcR(g0, s, t, newsize);
        r += tmpr * omega;
        // cout<<omega<<' '<<tmpr<<' '<<r<<endl;
    }
    return r;
}

double BSSERD(Graph &g, int &s, int &t, int size){
    double erd = 0;
    map<int, double> m = g.outneighbors[s];
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0||s==t) return 0;
    double omega, Omega=0;
    int newsize;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[s][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[s][it.first] = 1;
                break;
            }
            j++;
        }
        newsize = floor(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        double tmpe = nmcERD(g0, s, t, newsize);
        if(tmpe==0) Omega -= omega;
        else erd += tmpe*omega;
        // cout<<endl;
    }
    // cout<<"Numerator: "<<res[1]<<", denominator: "<<Omega<<endl;
    if(Omega!=0) erd = erd/Omega;
    return erd;
}

double BSSDCR(Graph &g, int &s, int &t, int &d, int size){
    double dcr = 0;
    map<int, double> m = g.outneighbors[s];
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0) return dcr;
    if(s==t) return 1;
    double omega, Omega=0;
    int newsize;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[s][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[s][it.first] = 1;
                break;
            }
            j++;
        }
        newsize = floor(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        double tmpr = nmcDCR(g0, s, t, d, newsize);
        dcr += tmpr*omega;
    }
    return dcr;
}


vector<double> BSSQUERY(Graph &g, int s, int t, int d, int size){
    // cout<<"Size: "<<size<<endl;
    vector<double> res = vector<double>(3,0);
    clock_t t1, t2;
    double r, e, c;
    if (s>=g.n || t>=g.n)
        handle_error("Vertex does not exist.");
    t1 = clock();
    r = BSSR(g, s, t, size);
    t2 = clock();
    g.BSS_R += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Reliability query "<<s<<' '<<t<<": "<<r<<endl;

    t1 = clock();
    e = BSSERD(g, s, t, size);
    // double e1 = ERD_min(g,s,t);
    // double e2 = ERD_max(g,s,t);
    t2 = clock();
    g.BSS_ERD += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Expected Reliable Distance query "<<s<<' '<<t<<": "<<e[0]<<','<<e[1]<<endl;

    t1 = clock();
    c = BSSDCR(g, s, t, d, size);
    t2 = clock();
    g.BSS_DCR += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Distance Constrained Reliability query "<<s<<' '<<t<<": "<<c<<endl;
    // cout<<"Reliability: "<<r<<"\nExpected Reliable Distance: "<<e[0]<<", "<<e[1]<<"\nDistance Constrained Reliability: "<<c<<endl;
    res[0] = r;
    res[1] = e;
    res[2] = c;
    return res;
}


void BSS_(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"BSSresults.txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    // clock_t t0, t1;
    // t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = BSSQUERY(g, s, t, d, K);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==50)break;
    }
    // t1 = clock();
    cout<<"Average BSS query time:"<<(g.BSS_R+g.BSS_ERD+g.BSS_DCR)/cnt<<"s"<<endl;
    cout<<"Reliability query: "<<g.BSS_R/cnt<<"s"<<endl;
    cout<<"ERD query: "<<g.BSS_ERD/cnt<<"s"<<endl;
    cout<<"DCR query: "<<g.BSS_DCR/cnt<<"s"<<endl;
    file.close();
    of.close();
}





double RSSR(Graph &g, int &s, int &t, int size, int b){
    static int threshold = 10;
    if(size<=threshold||g.RSSvisit[b]) return nmcR(g, s, t, size);
    double r = 0;
    map<int, double> m = g.outneighbors[b];
    g.RSSvisit[b] = true;
    int num = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<num<<", select "<<b<<endl;
    if(size==0||num==0) return r;
    if(s==t) return 1;
    double omega, Omega=0;
    int newsize, base;
    Graph g0 = g;
    // cout<<s<<' '<<t<<' '<<size<<endl;
    for(int i=0; i<num; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[b][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[b][it.first] = 1;
                base = it.first;
                break;
            }
            j++;
        }
        newsize = ceil(size * omega);
        // cout<<newsize<<' '<<size<<' ';
        if(newsize==0) continue;
        Omega += omega;
        // double tmpr = nmcR(g0, s, t, newsize);
        double tmpr = RSSR(g0, s, t, newsize, base);
        r += tmpr * omega;
        // cout<<omega<<' '<<tmpr<<' '<<r<<endl;
    }
    return r;
}

double RSSERD(Graph &g, int &s, int &t, int size, int b){
    static int threshold = 10;
    if(size<=threshold||g.RSSvisit[b]) return nmcERD(g, s, t, size);
    double erd = 0;
    map<int, double> m = g.outneighbors[b];
    g.RSSvisit[b] = true;
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0||s==t) return 0;
    if(size<=threshold) return nmcERD(g, s, t, size);
    double omega, Omega=0;
    int newsize, base;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[b][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[b][it.first] = 1;
                base = it.first;
                break;
            }
            j++;
        }
        newsize = ceil(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        double tmpe = RSSERD(g0, s, t, newsize, base);
        if(tmpe==0) Omega -= omega;
        else erd += tmpe*omega;
        // cout<<endl;
    }
    // cout<<"Numerator: "<<res[1]<<", denominator: "<<Omega<<endl;
    if(Omega!=0) erd = erd/Omega;
    return erd;
}

double RSSDCR(Graph &g, int &s, int &t, int &d, int size, int b){
    static int threshold = 10;
    if(size<=threshold||g.RSSvisit[b]) return nmcDCR(g, s, t, d, size);
    double dcr = 0;
    map<int, double> m = g.outneighbors[b];
    g.RSSvisit[b] = true;
    int r = m.size();
    // cout<<"Query "<<s<<" and "<<t<<": stratums "<<r<<endl;
    if(size==0||r==0) return dcr;
    if(s==t) return 1;
    if(size<=threshold) return nmcDCR(g, s, t, d, size);
    double omega, Omega=0;
    int newsize, base;
    Graph g0 = g;
    for(int i=0; i<r; i++){
        g0 = g;
        omega = 1;
        int j = 0;
        for(auto it: m){
            if(j<i){
                omega *= 1-it.second;
                g0.outneighbors[b][it.first] = 0;
            }
            else if(j==i){
                omega *= it.second;
                g0.outneighbors[b][it.first] = 1;
                base = it.first;
                break;
            }
            j++;
        }
        newsize = floor(size * omega);
        if(newsize==0) continue;
        Omega += omega;
        double tmpr = RSSDCR(g0, s, t, d, newsize, base);
        dcr += tmpr*omega;
    }
    return dcr;
}


vector<double> RSSQUERY(Graph &g, int s, int t, int d, int size){
    // cout<<"Size: "<<size<<endl;
    vector<double> res = vector<double>(3,0);
    clock_t t1, t2;
    double r, e, c;
    if (s>=g.n || t>=g.n)
        handle_error("Vertex does not exist.");
    t1 = clock();
    r = RSSR(g, s, t, size, s);
    t2 = clock();
    g.RSS_R += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Reliability query "<<s<<' '<<t<<": "<<r<<endl;

    t1 = clock();
    e = RSSERD(g, s, t, size, s);
    // double e1 = ERD_min(g,s,t);
    // double e2 = ERD_max(g,s,t);
    t2 = clock();
    g.RSS_ERD += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Expected Reliable Distance query "<<s<<' '<<t<<": "<<e[0]<<','<<e[1]<<endl;

    t1 = clock();
    c = RSSDCR(g, s, t, d, size, s);
    t2 = clock();
    g.RSS_DCR += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Distance Constrained Reliability query "<<s<<' '<<t<<": "<<c<<endl;
    // cout<<"Reliability: "<<r<<"\nExpected Reliable Distance: "<<e[0]<<", "<<e[1]<<"\nDistance Constrained Reliability: "<<c<<endl;
    res[0] = r;
    res[1] = e;
    res[2] = c;
    return res;
}


void RSS_(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"RSSresults.txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    // clock_t t0, t1;
    // t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = RSSQUERY(g, s, t, d, K);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==20)break;
    }
    // t1 = clock();
    cout<<"Average RSS query time:"<<(g.RSS_R+g.RSS_ERD+g.RSS_DCR)/cnt<<"s"<<endl;
    cout<<"Reliability query: "<<g.RSS_R/cnt<<"s"<<endl;
    cout<<"ERD query: "<<g.RSS_ERD/cnt<<"s"<<endl;
    cout<<"DCR query: "<<g.RSS_DCR/cnt<<"s"<<endl;
    file.close();
    of.close();
}




double LPSR(Graph &g, int &s, int &t, int size, int b){
    double r = 0;
    deque<int> q;
    vector<int> cv = vector<int>(g.n, -1);
    vector<deque<vector<double>>> hv = vector<deque<vector<double>>>(g.n ,deque<vector<double>>());
    vector<bool> visited = vector<bool>(g.n, false);
    vector<bool> init = vector<bool>(g.n, false);
    vector<int> reach;
    bool flag = false;
    if(s==t){
        return 1;
    }
    for(int i = 0; i < size; i++){
        fill(visited.begin(),visited.end(),false);
        q.clear();
        q.push_back(s);
        visited[s] = true;
        flag = false;
        while (!q.empty()){
            int node = q.front();
            q.pop_front();
            if (!init[node]){
                cv[node] = 0;
                hv[node].clear();
                for(auto nbr: g.outneighbors[node]){
                    int geo = Geo(nbr.second);
                    hv[node].push_back({(double)nbr.first, (double)(geo+cv[node]), nbr.second});
                }
                init[node] = true;
                sort(hv[node].begin(), hv[node].end(), compareHv);
            }
            if(g.outdegree[node]==0) {cv[node] += 1; continue;}
            while ((int)hv[node].front()[1]==cv[node]){
                vector<double> it = hv[node].front();
                hv[node].pop_front();
                // cout<<"Pop an Hv pair: "<<it[0]<<' '<<it[1]<<endl;
                if(!visited[it[0]]){
                    q.push_back((int)it[0]);
                    visited[(int)it[0]] = true;
                }
                int newX = Geo(it[2]);
                hv[node].push_back({it[0], (double)(newX+cv[node]+1), it[2]});
                sort(hv[node].begin(), hv[node].end(), compareHv);
                if(it[0]==t){
                    flag = true;
                }
            }
            cv[node] += 1;
            if(flag) break;
        }
        
    }
    r = 1.0*reach.size()/size;
    return r;
}

double LPSERD(Graph &g, int &s, int &t, int size, int b){
    double erd = 0;
    deque<int> q;
    vector<int> cv = vector<int>(g.n, -1);
    vector<deque<vector<double>>> hv = vector<deque<vector<double>>>(g.n ,deque<vector<double>>());
    vector<bool> visited = vector<bool>(g.n, false);
    vector<bool> init = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, -1);
    vector<int> reach;
    bool flag = false;
    if(s==t){
        return 0;
    }
    for(int i = 0; i < size; i++){
        fill(visited.begin(),visited.end(),false);
        fill(dis.begin(), dis.end(), -1);
        q.clear();
        q.push_back(s);
        visited[s] = true;
        flag = false;
        dis[s] = 0;
        while (!q.empty()){
            int node = q.front();
            q.pop_front();
            if (!init[node]){
                // cout<<"Initiating node "<<node<<endl;
                cv[node] = 0;
                hv[node].clear();
                for(auto nbr: g.outneighbors[node]){
                    int geo = Geo(nbr.second);
                    hv[node].push_back({(double)nbr.first, (double)(geo+cv[node]), nbr.second});
                    // cout<<"Create new pair "<<nbr.first<<' '<<geo<<endl;
                }
                init[node] = true;
                sort(hv[node].begin(), hv[node].end(), compareHv);
            }
            // cout<<"Now visiting node "<<node<<" in round "<<i<<", its cv now is: "<<cv[node]<<", out_degree is "<<g.outdegree[node]<<endl;
            // cout<<"Its neighbors(hv) now is: "<<endl;
            // for(auto it:hv[node]){
            //     cout<<it[0]<<' '<<it[1]<<endl;
            // }
            if(g.outdegree[node]==0) {cv[node] += 1; continue;}
            while ((int)hv[node].front()[1]==cv[node]){
                vector<double> it = hv[node].front();
                hv[node].pop_front();
                // cout<<"Pop an Hv pair: "<<it[0]<<' '<<it[1]<<endl;
                if(!visited[it[0]]){
                    dis[(int)it[0]] = dis[node]+1;
                    q.push_back((int)it[0]);
                    visited[(int)it[0]] = true;
                    // cout<<"Round "<<i<<", visiting "<<node<<", push "<<it[0]<<" in hv, its cv: "<<cv[it[0]]<<endl;
                    // cout<<"Push "<<it[0]<<" in hv, its cv: "<<cv[it[0]]<<", Round "<<i<<endl;
                }
                int newX = Geo(it[2]);
                hv[node].push_back({it[0], (double)(newX+cv[node]+1), it[2]});
                // cout<<"Round "<<i<<", change "<<node<<" of neighbor "<<it[0]<<" from "<<it[1]<<" to "<<(double)(newX+cv[node]+1)<<endl;
                sort(hv[node].begin(), hv[node].end(), compareHv);
                if(it[0]==t){
                    flag = true;
                    reach.push_back(dis[t]);
                    // cout<<"Found on round "<<i<<endl;
                    // break;
                }
            }
            cv[node] += 1;
            if(flag) break;
            // cout<<"Node is "<<node<<", now cv is "<<cv[node]<<endl;
        }
        
    }
    int sum = 0;
    for(int i = 0; i<reach.size(); i++){
        sum += reach[i];
    }
    if(!reach.size()) erd = 0;
    else erd = 1.0*sum/reach.size();
    return erd;
}

double LPSDCR(Graph &g, int &s, int &t, int &d, int size, int b){
    double dcr = 0;
    deque<int> q;
    vector<int> cv = vector<int>(g.n, -1);
    vector<deque<vector<double>>> hv = vector<deque<vector<double>>>(g.n ,deque<vector<double>>());
    vector<bool> visited = vector<bool>(g.n, false);
    vector<bool> init = vector<bool>(g.n, false);
    vector<int> dis = vector<int>(g.n, -1);
    vector<int> reach;
    bool flag = false;
    if(s==t){
        return 1;
    }
    for(int i = 0; i < size; i++){
        fill(visited.begin(),visited.end(),false);
        fill(dis.begin(), dis.end(), -1);
        q.clear();
        q.push_back(s);
        visited[s] = true;
        flag = false;
        dis[s] = 0;
        while (!q.empty()){
            int node = q.front();
            q.pop_front();
            if (!init[node]){
                // cout<<"Initiating node "<<node<<endl;
                cv[node] = 0;
                hv[node].clear();
                for(auto nbr: g.outneighbors[node]){
                    int geo = Geo(nbr.second);
                    hv[node].push_back({(double)nbr.first, (double)(geo+cv[node]), nbr.second});
                    // cout<<"Create new pair "<<nbr.first<<' '<<geo<<endl;
                }
                init[node] = true;
                sort(hv[node].begin(), hv[node].end(), compareHv);
            }
            if(g.outdegree[node]==0) {cv[node] += 1; continue;}
            while ((int)hv[node].front()[1]==cv[node]){
                vector<double> it = hv[node].front();
                hv[node].pop_front();
                // cout<<"Pop an Hv pair: "<<it[0]<<' '<<it[1]<<endl;
                if(!visited[it[0]]){
                    dis[(int)it[0]] = dis[node]+1;
                    q.push_back((int)it[0]);
                    visited[(int)it[0]] = true;
                }
                int newX = Geo(it[2]);
                hv[node].push_back({it[0], (double)(newX+cv[node]+1), it[2]});
                // cout<<"Round "<<i<<", change "<<node<<" of neighbor "<<it[0]<<" from "<<it[1]<<" to "<<(double)(newX+cv[node]+1)<<endl;
                sort(hv[node].begin(), hv[node].end(), compareHv);
                if(it[0]==t){
                    flag = true;
                    reach.push_back(dis[t]);
                    // cout<<"Found on round "<<i<<endl;
                    // break;
                }
            }
            cv[node] += 1;
            if(flag) break;
            // cout<<"Node is "<<node<<", now cv is "<<cv[node]<<endl;
        }
        
    }
    int cnt = 0;
    for(int i = 0; i<reach.size(); i++){
        if(reach[i]<=d) cnt++;
    }
    dcr = 1.0*cnt/size;
    return dcr;
}


vector<double> LPSQUERY(Graph &g, int s, int t, int d, int size){
    // cout<<"Size: "<<size<<endl;
    vector<double> res = vector<double>(3,0);
    clock_t t1, t2;
    double r, e, c;
    if (s>=g.n || t>=g.n)
        handle_error("Vertex does not exist.");
    t1 = clock();
    r = LPSR(g, s, t, size, s);
    t2 = clock();
    g.RSS_R += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Reliability query "<<s<<' '<<t<<": "<<r<<endl;

    t1 = clock();
    e = LPSERD(g, s, t, size, s);
    // double e1 = ERD_min(g,s,t);
    // double e2 = ERD_max(g,s,t);
    t2 = clock();
    g.RSS_ERD += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Expected Reliable Distance query "<<s<<' '<<t<<": "<<e[0]<<','<<e[1]<<endl;

    t1 = clock();
    c = LPSDCR(g, s, t, d, size, s);
    t2 = clock();
    g.RSS_DCR += (double)(t2-t1)/CLOCKS_PER_SEC;
    // cout<<"Distance Constrained Reliability query "<<s<<' '<<t<<": "<<c<<endl;
    // cout<<"Reliability: "<<r<<"\nExpected Reliable Distance: "<<e[0]<<", "<<e[1]<<"\nDistance Constrained Reliability: "<<c<<endl;
    res[0] = r;
    res[1] = e;
    res[2] = c;
    return res;
}


void LPS_(Graph &g, int K, int d){
    ifstream file;
    file.open(g.folder+g.graph_file+"queries.txt");
    ofstream of(g.folder+g.graph_file+"LPSresults.txt",ios::ate);
    if(!file)
        handle_error("open");
    string str;
    int s, t;
    int cnt = 0;
    // clock_t t0, t1;
    // t0 = clock();
    while(getline(file,str)){
        int size=str.size();
        for(int i=0;i<size;i++){
            if(str[i]==' '){
                s=atoi(str.substr(0,i).c_str());
                t=atoi(str.substr(i+1,size).c_str());
                break;
            }
        }
        vector<double> res = LPSQUERY(g, s, t, d, K);
        of<<res[0]<<' '<<res[1]<<' '<<res[2]<<endl;
        cnt++;
        // if(cnt>2)break;
        // if(cnt%50==0)
        cout<<cnt<<" Query between "<<s<<" and "<<t<<", Reliability: "<<res[0]<<", ERD: "<<res[1]<<", DCR: "<<res[2]<<endl;
        if(cnt==20)break;
    }
    // t1 = clock();
    cout<<"Average RSS query time:"<<(g.RSS_R+g.RSS_ERD+g.RSS_DCR)/cnt<<"s"<<endl;
    cout<<"Reliability query: "<<g.RSS_R/cnt<<"s"<<endl;
    cout<<"ERD query: "<<g.RSS_ERD/cnt<<"s"<<endl;
    cout<<"DCR query: "<<g.RSS_DCR/cnt<<"s"<<endl;
    file.close();
    of.close();
}