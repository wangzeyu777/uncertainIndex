class Edge{
    public:
        int start, end;
        double prob;
        int d;
    Edge(int u, int v, double p){
        start=u; end=v; prob=p; d=0;
    }
    Edge(){}
};