#include <iostream>
#include <vector>
#include <map>

using namespace std;

class IndexEntry{
    public:
        int node;
        double erd, reliability;
        // int k;  // Distance constraint; 
        vector<vector<double>> dcr;
        bool operator == (const int &x) //查找数值x是否与node相等
        {
            return (this->node == x);
        }
    IndexEntry(int NODE, double R, double E, vector<vector<double>> DCR){
        node = NODE;
        reliability = R;
        erd = E;
        dcr = DCR;
    }
};