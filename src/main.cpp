// #pragma once
// #include <iostream>
#include "graph.h"
#include "sample.h"
#include "query.h"
// #include "memoryusage.h"
// #include "MeasureM.h"
// #include "independent_samples.h"
// using namespace std;

int main(){

    clock_t start_time, end_time, mc_start, mc_end, index_start, index_end, vc_start, vc_end;   
    start_time=clock();
    Graph g("epi/");
    // generate_queries(g, 1000);
    // // g.statistics();
    // // g.degreeSort();
    // // g.edgeSort();
    // // for(int i=0;i<10;i++){
    // //     vector<int> vc=vertexCover(g,1);
    // // }
    int K = 1000, h = 1, dc = 6;

    // Graph g0 = test(g);


    // nmc(g, 10000, dc);
    // BSS_(g, 1000, dc);
    // RSS_(g, 1000, dc);
    // BSS(g, 1000, dc);
    // RSS(g, 1000, dc);

    // while (true){
    //     int s, t;
    //     cout<<"Input your query vertices: ";
    //     cin>>s>>t;
    //     if (s == -1 && t == -1) break;
    //     clock_t query_start, query_end;
    //     query_start = clock();
    //     vector<double> nmc = nmc_query(g, s, t, dc, 10000);
    //     cout<<nmc[0]<<' '<<nmc[1]<<' '<<nmc[2]<<endl;
    //     query_end = clock();
    //     cout<<"Query time: "<<(double)(query_end-query_start)/CLOCKS_PER_SEC<<"s"<<endl;
    // }

    vc_start = clock();
    vector<int> vc=vertexCover(g,h);
    vc_end = clock();
    cout<<"Constructing vertex cover time: "<<(double)(vc_end-vc_start)/CLOCKS_PER_SEC<<"s"<<endl;

    // 分别进行正向、反向采样，构造采样结果
    mc_start = clock();
    vector<map<int,vector<int>>> reach = MonteCalo(true, g, vc, K, h);
    vector<map<int,vector<int>>> revReach = MonteCalo(false, g, vc, K, h);
    mc_end = clock();
    cout<<"Monte Calo time: "<<(double)(mc_end-mc_start)/CLOCKS_PER_SEC/(2*K)<<"s"<<endl;

    index_start = clock();
    constructIndex(true, g, reach, K, h, dc);
    constructIndex(false, g, revReach, K, h, dc);
    index_end = clock();
    cout<<"Construct index time: "<<(double)(index_end-index_start)/CLOCKS_PER_SEC<<"s"<<endl;

    end_time = clock();
    cout<<"Time spent: "<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    indexQuery(g, dc, h);

// // // // // // // // //根据需求设置一些可能的输出函数用于进行性能分析
//     indexOutput(g, h);
    // readIndex(g, h, dc);
    // indexQuery(g, dc, h);



    // string ss = "-0.1234";
    // cout<<atof(ss.c_str())<<endl;
    // cout<<ss.substr(0,3)<<endl;
    // cout<<ss.substr(0,ss.size())<<endl;
    // cout<<ss.substr(2,ss.size())<<endl;
    // cout<<atoi(ss.c_str())<<endl;
    // cout<<to_string(atof(ss.c_str()))<<endl;
    // cout<<to_string(atoi(ss.c_str()))<<endl;


    return 0;
}