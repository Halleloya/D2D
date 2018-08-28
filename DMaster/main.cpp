#include<iostream>
#include <fstream>
#include <queue>
#include <ctime>
#include   <stdlib.h>
using namespace std;

ifstream fin;
ofstream fout;

int DeviceNum, LinkNum, SterileNum; //后sterileNum个为不孕结点
bool edge[2000][2000]={{0}};
bool inQ[2000][2000]={{0}};  //判断边是否已经在 Q 中
bool T[2000][2000]={{0}};   //就是原始论文中的 T
int fv[2000][1000]={{0}};  //就是fv
int fsize[2000]={0};

int Iv[2000][1000]={{0}};
int f,Lsize;    // f
bool W[1000]={0},S[1000]={0};     //集合W  S
int d[1000];  // T 的度
int Imin =0;
int L[20000]={0};

int oldgen=0;
int newRoot;


queue<int> *D;

struct Qnode
{
    int x,y;
    Qnode(int i,int j){x=i;y=j;}
};

queue<Qnode> Q;

struct node
{
    node* parent, *p, *lca;
    Qnode *New, *old;
    int rank,depth,vi;

    node(int a)
    {vi=a; rank=depth=0; p=parent=lca=NULL; New=old=NULL;}
    ~node(){}

};

node * vertex[10000];

void Reduce(node*x);
void MakeSet(node *x);
node *FindSet(node *x);
node *LCA(node *x);
void Link(node* x,node* y);
void Union(node*x, node*y);
void Initialization();
void BasicMerge(node *x);
node *CycleMerge(node*x, node*y);
bool SatisfyEndCondition();
void RefineTreeLoad();
void bfs_tree();
void Tbfs(int v);
void MinSPTLoad();

int dfs(int cNode, int c, bool* inT);
void show(bool TT[][2000]);
void dfs_tree();

void MakeSet(node  *x)
{
    x->p = x;
    x->rank = 0;
    x->lca = x;
}

node *FindSet(node *x)
{
    if (x != x->p)  x->p = FindSet(x->p);
    return x->p;
}

node *LCA(node *x)
{
    node *y = FindSet(x);
    return y->lca;
}

void Link(node* x,node* y)
{
    node *c = x->lca;
    if (c->depth > y->lca->depth)
        c = y->lca;

    if (x->rank > y->rank)
    {
        y->p = x;
        x->lca= c;
    }
    else
    {
        x->p = y;
        y->lca= c;
        if (x->rank == y->rank)
            y->rank = y->rank +1;
    }
}

void Union(node*x, node*y)
{
    Link(FindSet(x), FindSet(y));
}


void Tbfs(int v)
{
    bool inT[2000]={0};
    int cNode;
    queue<int> qe;

for (int i=0; i<DeviceNum; i++) vertex[i]->depth=-1;


    qe.push(v);
    inT[v]=1;
    vertex[v]->depth=0;
    vertex[v]->parent= vertex[v];


    while(!qe.empty())
    {
        cNode = qe.front(); qe.pop();

        for (int j=0; j<DeviceNum; j++)
        {
            if ( T[cNode][j] && (inT[j] ==false ) && j!=cNode )
            {
                inT[j]=1;
                qe.push(j);
                vertex[j]->depth= vertex[cNode]->depth+1;
                vertex[j]->parent= vertex[cNode];
            }
        }
    }

}

void Initialization()
{

    for ( newRoot =0; newRoot<DeviceNum; newRoot++)   //select a vertex in S
        if(S[newRoot]) { break; }

    oldgen=newRoot;

    Tbfs(newRoot);

    for (int x=0 ; x<DeviceNum; x++)
    {
        vertex[x]->New = vertex[x]->old = NULL;
    }
//    select a vertex in S, and perform BFS from this
//    vertex to establish the depth and parent of each
//    vertex. Set each vertex¨s New and old pointers
//    to NULL.

    for (int x=0 ; x<DeviceNum; x++)
        if ( !W[x] )
            MakeSet(vertex[x]);

  while (!Q.empty())  Q.pop();

    for (int x=0;x<DeviceNum; x++)
        for (int y=x+1; y<DeviceNum; y++)
        {
              inQ[x][y] = inQ[y][x]=0;

              if(edge[x][y] && !W[x] && !W[y] )
              {
                  if (T[x][y])
                        {
                            Union(vertex[x],vertex[y]);
                        }
                  else
                  {
                      if (!inQ[x][y])
                      {
                          Q.push(*new Qnode(x,y));
                          inQ[x][y] = inQ[y][x]=1;
                      }
                  }
              }
        }

}

void RefineTreeLoad()
{
    for (int v=0; v<DeviceNum; v++)
    {
        W[v]=S[v]=0;

        if( fv[v][d[v]] > L[Imin]) { W[v]=1; }
        if( fv[v][d[v]]==f )  { S[v]=1; }
    }

    D= new queue<int> [Lsize];
    for (int i=0; i<DeviceNum; i++)
    {
        D[ Iv[i][d[i]] ].push(i);
    }

 //   1: W ={v |fv[v][d[v]] > L[Imin ]}
  //  2: S ={v |f v (d(v))=f}
  //  3: let D be an array consisting of |L| linked lists of vertices
  //  4: insert all vertices with load L[i] into list D[i] using the indexing lists I

    Initialization();

    node *w = NULL;
    while (w == NULL)
    {
        if (!Q.empty())
        {
            Qnode  xQ=Q.front();
            w = CycleMerge(vertex[xQ.x] , vertex[xQ.y]);

            inQ[xQ.x][xQ.y] = inQ[xQ.y][xQ.x] =0;
            Q.pop();
        }
        else
        {
            Imin = Imin +1;
            for (; !D[Imin].empty()  ; D[Imin].pop() )
            {
                int v=D[Imin].front();

                if ( W[v] )
                {
                    if (fv[v][d[v]+1]>=f )
                        {return; }
                    BasicMerge( vertex[v] );
                }
            }
        }
    }

    Reduce(w);

    for (int i=0; i<DeviceNum; i++)     //统计 T 的度
    {
        d[i]=0;
        for (int j=0 ; j<DeviceNum; j++)
            d[i] += T[i][j];
    }
   // d[newRoot]++;

    f = fv[0][d[0]];     //the load of tree T;
    for (int j=1 ; j<DeviceNum; j++)
    { if (f < fv[j][d[j]]) f=fv[j][d[j]]; }

    delete D;

}

void Reduce(node*x)
{
    if (x->old !=NULL && !T[x->New->y][x->New->x] ){
        T[x->New->x][x->New->y] = T[x->New->y] [x->New->x] =1;
        T[x->old->x][x->old->y] = T[x->old->y] [x->old->x] =0;

        delete x->old;
        x->old = NULL;

        Reduce( vertex[x->New->x] );
        Reduce( vertex[x->New->y] );
       // add edge x->new to the tree
       // remove edge x->old from the tree
     //   let u and v be the two ends of x->new
      //  call Reduce(u) and Reduce(v)
    }
}


node *CycleMerge(node*x, node*y)
{
    if (FindSet(x) == FindSet(y))
        return NULL;
    node *u = LCA(x);
    node *v = LCA(y);
    node *w = NULL;

    while ( u!=v )
    {
        if (u->depth >= v->depth)
        {
            if ( W[u->vi] )
            {
                BasicMerge(u);
                if ( S[u->vi] && w ==NULL )
                   { w = u; }
            }
            if ( W[u->parent->vi] )
            {
                u->parent->old= new Qnode(u->vi,u->parent->vi);
                u->parent->New= new Qnode(x->vi,y->vi);
                u=u->parent;
            }
            else
            {
                u= LCA(u->parent);
            }
        }
        else
        {
            if (W[v->vi])
            {
                BasicMerge(v);
                if ( S[v->vi] && w == NULL)
                    { w = v;}
            }
            if (W[v->parent->vi])
            {
                v->parent->old= new Qnode(v->vi,v->parent->vi);
                v->parent->New= new Qnode(x->vi,y->vi);
                v =v->parent;
            }
            else
            {
                v = LCA(v->parent);
            }
        }
    }

    if ( W[u->vi] )
    {
        BasicMerge(u);
        if ( S[u->vi] && w ==NULL )
            { w =u; }
    }

    return w;
}

void BasicMerge(node *x)
{
    W[x->vi]=0;    //remove x from W;
    MakeSet(x);
    for ( int y =0 ; y<DeviceNum; y++ )
    {
        if ( edge[x->vi][y] && !W[y] )
        {
            if ( T[x->vi][y] )
                { Union(x,vertex[y]); }
            else
                {
                    if ( !inQ[x->vi][y] )
                    {
                        Q.push( *new Qnode(x->vi,y));
                        inQ[x->vi][y] = inQ[y][x->vi]=1;
                    }
                }
        }
    }

}

void bfs_tree()
{
    bool inT[2000]={0};
    int cNode;
    queue<int> qe;

    int ran = rand()%DeviceNum;

    qe.push(ran);
    inT[ran]=1;
  //  cout<<"gen = "<<ran<<endl;

    while(!qe.empty())
    {
        cNode = qe.front(); qe.pop();

        for (int j=0; j<DeviceNum; j++)
        {
            if ( edge[cNode][j] && (inT[j] ==false ) && j!=cNode )
            {
                inT[j]=1;
                qe.push(j);
                T[cNode][j]=T[j][cNode]=1;
            }
        }
    }
}

int dfs(int cNode, int c, bool* inT)
{
    if (c==DeviceNum)
    {
        return 1;
    }

    int ran= rand()%DeviceNum;
    int ranNode = rand()%DeviceNum;
    bool inj[2000]={0};

    for (int j=0; j<DeviceNum; j++)
    {
        while( inj[ran] ) {ran= rand()%DeviceNum;}
        inj[ran] = 1;
//ran=j;
        if ( edge[cNode][ran] && (inT[ran] ==false ) && ran!=cNode )
        {
            inT[ran] = 1;
            T[cNode][ran]=T[ran][cNode]=1;

            while ( !inT[ranNode])
                {ranNode= rand()%DeviceNum;}
//ranNode =j;
            if ( ! dfs(ranNode,c+1,inT)){
                inT[ran] =0;
                T[cNode][ran]=T[ran][cNode]=0;
            }
            else return 1;
        }
    }

    return 0;
}

void dfs_tree()
{
    bool inT[2000]={0};
    int ran = rand()%DeviceNum;

    inT[ran]=1;
  //  cout<<"gen = "<<ran;
    dfs(ran,1,inT);
}

void MinSPTLoad()
{
    dfs_tree();   //1: construct an arbitrary spanning tree T of graph G

     //   2: sort in ascending order the set L in Fact 4
    int minfv= fv[0][0];
    int maxfv= fv[0][0];
    int Lt[20000]={0};
    for (int i=0; i<DeviceNum; i++)
        for (int j=0; j<fsize[i]; j++)
        {
            Lt[fv[i][j]]=1;
            if (minfv> fv[i][j]) minfv= fv[i][j];
            if (maxfv< fv[i][j]) maxfv= fv[i][j];
        }

    int k=0;
    for (int i=minfv; i<=maxfv; i++)
    {
        if ( Lt[i] ) {
            L[k]=i;
            Lt[i]=k++;
        }
    }
    Lsize = k;

//  3: build the indexing lists I v for all v such that I v (j) is the index of f v (j) in L
    for (int i=0; i<DeviceNum; i++)
        for (int j=0; j<fsize[i]; j++)
            Iv[i][j]=Lt[fv[i][j]];

    Imin=0;

    for (int i=0; i<DeviceNum; i++)    //统计结点的度
    {
        vertex[i]=new node(i);
        d[i] = 0;

        for (int j=0 ; j<DeviceNum; j++)
            d[i] += T[i][j];
    }
   // d[0]++;

    f = fv[0][d[0]];     // 生成 T 的负载;
    for (int j=1 ; j<DeviceNum; j++)
    { if (f < fv[j][d[j]]) f=fv[j][d[j]];  }

fout<<f;      //输出优化前树的初始负载

    while (SatisfyEndCondition()==false )
        { RefineTreeLoad(); }
}


bool SatisfyEndCondition()
{
    for (int i=0 ; i<DeviceNum ; i++)
        if ((f <=fv[i][d[i]+1] && fv[i][d[i]] <= L[Imin])
            || (f ==fv[i][d[i]]&&fv[i][d[i]-1] <= L[ Imin ]))
                return true;
    return false;
}


void show(bool TT[][2000])
{
    for(int i=0; i<DeviceNum; i++)
    {
        for(int j=0; j<DeviceNum; j++)
        {
            cout<<TT[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
}


int main(int argc,char *argv[])
{

    fin.open(argv[1]);
    fout.open(argv[2],ios::app);

    fin>>DeviceNum>>LinkNum>>SterileNum;

    int e;
/*
    for(int i=0; i<LinkNum; i++)
    {
        fin>>a>>b;
        edge[a][b] = 1;
        edge[b][a] = 1;
    }*/
    for(int i=0; i<DeviceNum; i++)
        for(int j=0; j<DeviceNum; j++)
        {
            if(i==j) edge[i][j] = 0;
            else if (i == j+1 || i==j-1)
                edge[i][j] = 0;
            else edge[i][j] = 1;
        }

    for(int j=0; j<DeviceNum; j++)
    {
        for(int k=0; k<=DeviceNum ;k++)  //第j号结点的度为k的负载函数
        {
            fin>>e;
            fv[j][k] = e;
            fsize[j]++;
        }
    }

    srand(unsigned(time(NULL)));
    MinSPTLoad();   ///////////////////// 进入优化主函数

fout<<'\t'<<f<<endl;    //输出优化后树的负载

    fin.close();
    fout.close();
    cout<<"Output successful!"<<endl;
    return 0;

}
