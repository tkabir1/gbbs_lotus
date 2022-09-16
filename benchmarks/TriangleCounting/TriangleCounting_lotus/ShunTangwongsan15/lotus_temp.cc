#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include<iostream>
#include<string.h>
#include <vector>
#include <stack>
#include <bitset>
#include <cassert>
#include "parlay/sequence.h"

using namespace std;
//Number of triangles	1612010
// Define the maximum number of vertices in the graph
#define N 4039
#define PERCEN 10
 
// Data structure to store a graph object
struct Node_Element
{
    int oldId;
    int new_vertexId;
    int degree;
};

struct Graph
{
    // An array of pointers to Node to represent an adjacency list
    struct Node* head[N];
    struct Node_Element nodes[N];
    int hub_vertex;
};

struct Lotus_Graph
{
    // An array of pointers to Node to represent an adjacency list
    struct Graph *hub;
    struct Graph *non_hub;
    bool **hub_array;
    vector <parlay::sequence<typename T, typename Allocator>> hubs[N];
    vector <int> non_hubs[N];
    int hub_count;
    
};

// Data structure to store adjacency list nodes of the graph

// Data structure to store adjacency list nodes of the graph
struct Node
{
    int dest;
    struct Node* next;
};
 
// Data structure to store a graph edge
struct Edge {
    int src, dest;
};
 
// Function to create an adjacency list from specified edges
struct Graph* createGraph(struct Edge edges[], int n)
{
    // allocate storage for the graph data structure
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
 
    // initialize head pointer for all vertices
    for (int i = 0; i < N; i++) {
        graph->head[i] = NULL;
        graph->nodes[i].oldId=i;
    }
 
    // add edges to the directed graph one by one
    cout<<"Number"<<n;
    for (int i = 0; i < n; i++)
    {
        // get the source and destination vertex
        //cout<<i<<"\n";
        int src = edges[i].src;
        //cout<<"In create"<<src<<" ends";
        int dest = edges[i].dest;
 
        // allocate a new node of adjacency list from src to dest
        //struct Node_Element* node = (struct Node_Element*)malloc(sizeof(struct Node_Element));
        //node->oldId=i;
        struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
        newNode->dest = dest;
 
        // point new node to the current head
        newNode->next = graph->head[src];
 
        // point head pointer to the new node
        graph->head[src] = newNode;

        struct Node* newNode2 = (struct Node*)malloc(sizeof(struct Node));
        newNode2->dest = src;
 
        // point new node to the current head
        newNode2->next = graph->head[dest];
 
        // point head pointer to the new node
        graph->head[dest] = newNode2;
    }
    cout<<"Creation done";
    return graph;
}
 
// Function to print adjacency list representation of a graph
void printGraph(struct Graph* graph)
{
    for (int i = 0; i < N; i++)
    {
        /*if (i==10){
            break;
        }*/
        // print current vertex and all its neighbors
        struct Node* ptr = graph->head[i];
        while (ptr != NULL)
        {
            printf("(%d —> %d)\t", i, ptr->dest);
            ptr = ptr->next;
        }
 
        printf("\n");
    }
}

void countEdgeGraph(struct Graph* graph)
{
    int edge=0;
    for (int i = 0; i < N; i++)
    {
        /*if (i==10){
            break;
        }*/
        // print current vertex and all its neighbors
        struct Node* ptr = graph->head[i];
        while (ptr != NULL)
        {
            //printf("(%d —> %d)\t", i, ptr->dest);
            ptr = ptr->next;
            edge++;
        }
 
        //printf("\n");
    }
    cout<<"\nEdge count"<<edge<<"\n";
}

int get_degree(struct Graph* graph,int target)
{
    int count=0;
    for (int i = 0; i < N; i++)
    {
        // print current vertex and all its neighbors
        struct Node* ptr = graph->head[i];
        //cout<<i;
        if (i==target){
            //cout<<i<<"\n";
            while (ptr != NULL)
            {
                //printf("(%d —> %d)\t", i, ptr->dest);
                //cout<<"Here";
                ptr = ptr->next;
                count++;
            }
            break;
        }
    }
    return count;
}
int comparator(const void *p, const void *q) 
{ 
    return (*(int*)q)-(*(int*)p);
} 

int findIndex(int arr[], int idx, int K)
{
 
    // Base Case
    if (idx < 0)
        return -1;
 
    // Return Statement
    if (arr[idx] == K) {
        return idx;
    }
 
    // Recursive Call
    return findIndex(arr, idx - 1, K);
}
struct Graph* sort_graph(struct Graph* graph){
    int degrees[N];
    for (int i = 0; i < N; i++) {
        degrees[i]=get_degree(graph,i);
        graph->nodes[i].degree=degrees[i];
        //cout<<i<<" with degree "<<get_degree(graph,i)<<"\n";
    }
    qsort((void*)degrees, N, sizeof(degrees[0]), comparator); 
    float r_f=(PERCEN*1.0/100.0)*N;
    int r=(int)r_f;
    cout<<"\n"<<r<<"r";
    int max_degree=degrees[0];
    int min_degree_taken=degrees[r];
    r=findIndex(degrees,N, min_degree_taken);
    cout<<"\nmax:"<<max_degree<<"min"<<min_degree_taken<<"\n";
    cout<<"Updated"<<r<<"\n";
    int assigned_top=0;
    int assigned_normal=r+1;
    graph->hub_vertex=r;
    for (int i=0;i<N;i++){
        if (graph->nodes[i].degree>=min_degree_taken && assigned_top<=r){
            graph->nodes[i].new_vertexId=assigned_top;
            assigned_top++;
        }
        else{
            graph->nodes[i].new_vertexId=assigned_normal;
            assigned_normal++;
        }
    }
    /*for (int i=0;i<N;i++)
    {
        cout<<"old degree"<<graph->nodes[i].oldId<<" and new id "<<graph->nodes[i].new_vertexId<<"\n";
    }*/
    return graph;
}

/*struct Graph* recreateGraph(struct Graph* grpah_in)
{
    // allocate storage for the graph data structure
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
 
    // initialize head pointer for all vertices
    for (int i = 0; i < N; i++) {
        graph->head[i] = NULL;
        graph->nodes[i].oldId=i;
    }
 
    // add edges to the directed graph one by one
    cout<<"Number"<<n;
    for (int i = 0; i < n; i++)
    {
        // get the source and destination vertex
        //cout<<i<<"\n";
        int src = edges[i].src;
        //cout<<"In create"<<src<<" ends";
        int dest = edges[i].dest;
 
        // allocate a new node of adjacency list from src to dest
        //struct Node_Element* node = (struct Node_Element*)malloc(sizeof(struct Node_Element));
        //node->oldId=i;
        struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
        newNode->dest = dest;
 
        // point new node to the current head
        newNode->next = graph->head[src];
 
        // point head pointer to the new node
        graph->head[src] = newNode;
    }
 
    return graph;
}*/

struct Graph* create_relabeling_array(struct Graph* graph){
    int target;
    //scanf("%d",&target);
    //cout<<"Enter node number";
    //cin>>target;
    //cout<<get_degree(graph,target);
    //printf("%d",get_degree(graph,target));
    struct Graph *graph_sort = sort_graph(graph);
    return graph;
}

struct Graph* Lotus_Preprocessing(struct Graph* graph){
    int hub_count=graph->hub_vertex;
    cout<<"hub_vertex"<<graph->hub_vertex;
    int hub_array[hub_count]; //need optimization for bit
    struct Graph *graph_hub = (struct Graph*)malloc(sizeof(struct Graph));
    struct Graph *graph_nonhub = (struct Graph*)malloc(sizeof(struct Graph));
    struct Graph* lotus_graph = (struct Graph*)malloc(sizeof(struct Graph));
    lotus_graph->hub_vertex=graph->hub_vertex;
    for (int i = 0; i < N; i++) {
        lotus_graph->head[i] = NULL;
    }
    vector <int> hubs[N];
    vector <int> non_hubs[N];
    int p_vertex=0;
    for (int i=0;i<N;i++){
        int vertex=graph->nodes[i].new_vertexId;
        struct Node* ptr = graph->head[i];
        while (ptr != NULL)
        {
            int src=vertex;
            int dest=graph->nodes[ptr->dest].new_vertexId;
            struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
            newNode->dest = dest;
 
            // point new node to the current head
            newNode->next = lotus_graph->head[src];
 
            // point head pointer to the new node
            lotus_graph->head[src] = newNode;
            ptr = ptr->next;
        }
    }

    /*int p_vertex=0;
    for (int i=0;i<N;i++){
        int vertex=graph->nodes[i].new_vertexId;
        if (vertex==p_vertex){
            cout<<" "<<vertex<<",";
            struct Node* ptr = graph->head[vertex];
            while (ptr != NULL)
            {
                printf("(%d —> %d)\t", i, ptr->dest);
                int new_dest=graph->nodes[ptr->dest].new_vertexId;
                cout<<"after change";
                printf("(%d —> %d)\t", vertex, new_dest);
                ptr = ptr->next;
            }
        }
        p_vertex++;
    }*/
    

    return lotus_graph;
}

void print_vector(std::vector <int> const &a) {
   std::cout << "The vector elements are : ";

   for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
}
void printarray( int **array, int SIZE ){
    int i;
    int j;
    cout<<"printing"<<"\n\n";
    for( j = 0; j < SIZE; j++ ){
        for( i = 0; i < SIZE; i ++){
            printf( "%d ", array[j][i] );
        }
        printf( "\n" );
    }
}

struct Graph* createSubGraph(vector<int> *array)
{
    // allocate storage for the graph data structure
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
 
    // initialize head pointer for all vertices
    for (int i = 0; i < N; i++) {
        graph->head[i] = NULL;
        //graph->nodes[i].oldId=i;
    }
 
    // add edges to the directed graph one by one
    //cout<<"Number"<<n;
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j < array[i].size(); j++){
            int src = i;
            //cout<<"In create"<<src<<" ends";
            int dest = array[i].at(j);
    
            // allocate a new node of adjacency list from src to dest
            //struct Node_Element* node = (struct Node_Element*)malloc(sizeof(struct Node_Element));
            //node->oldId=i;
            struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
            newNode->dest = dest;
    
            // point new node to the current head
            newNode->next = graph->head[src];
    
            // point head pointer to the new node
            graph->head[src] = newNode;
        }
    }
 
    return graph;
}

struct Lotus_Graph* Lotus_Preprocessing_neighbours(struct Graph* graph){
    int hub_count=graph->hub_vertex+1;
    cout<<"hub_vertex"<<graph->hub_vertex;
    //int hub_array[hub_count][hub_count]; //need optimization for bit
    bool** hub_array = (bool**)malloc(hub_count * sizeof(bool*));
    /*for (int i = 0; i < hub_count; i++){
        hub_array[i] = (bool*)malloc(hub_count* sizeof(bool));
    }*/
    for (int i = 0; i < hub_count; i++){
        hub_array[i] = (bool*)malloc(i* sizeof(bool));
    }
    struct Graph *graph_hub = (struct Graph*)malloc(sizeof(struct Graph));
    struct Graph *graph_nonhub = (struct Graph*)malloc(sizeof(struct Graph));
    struct Graph* lotus_graph = (struct Graph*)malloc(sizeof(struct Graph));
    struct Lotus_Graph* lotus_pre = (struct Lotus_Graph*)malloc(sizeof(struct Lotus_Graph));
    //vector <int> hubs[N];
    //vector <int> non_hubs[N];
    for (int i = 0; i < N; i++)
    {
        // print current vertex and all its neighbors
        struct Node* ptr = graph->head[i];
        for (;ptr != NULL;ptr = ptr->next)
        {
            //printf("(%d —> %d)\t", i, ptr->dest);
            if (i==ptr->dest){
                continue;
            }
            if (ptr->dest>i){
                continue;
            }
            if (ptr->dest<hub_count){
                lotus_pre->hubs[i].push_back(ptr->dest);
                if (i<hub_count){
                    hub_array[i][ptr->dest]=1;
                }
            }
            else{
                lotus_pre->non_hubs[i].push_back(ptr->dest);
            }
        }
        //print_vector(non_hubs[1]);
        //printf("\n");
    }
    graph_hub=createSubGraph(lotus_pre->hubs);
    graph_nonhub=createSubGraph(lotus_pre->non_hubs);
    /*cout<<"\nHub graph\n";
    printGraph(graph_hub);
    cout<<"\nnonHub graph\n";
    printGraph(graph_nonhub);*/
    lotus_pre->hub=graph_hub;
    lotus_pre->non_hub=graph_nonhub;
    lotus_pre->hub_array=hub_array;
    /*cout<<"\nhub_array\n";
    printarray(lotus_pre->hub_array, (hub_count));*/
    lotus_pre->hub_count=hub_count;

    return lotus_pre;
}
int hhnh_traingle(struct Lotus_Graph* graph){
    int triangle=0;
    int start=graph->hub_count;
    cout<<"\nhub_count"<<start;
    for (int i = 0; i < N; i++)
    {
        // print current vertex and all its neighbors
        int v=i;
        //cout<<"vertex"<<v<<"\n";
        struct Node* ptr = graph->hub->head[v];
        vector <int> n_v=graph->hubs[i];
        sort(n_v.begin(), n_v.end());
        //print_vector(n_v);
        for(int i=0; i < n_v.size(); i++){
            int first=n_v.at(i);
            for (int j=0;j<n_v.size();j++){
                int second=n_v.at(j);
                //cout<<"first"<<first<<"second"<<second<<"\n";
                if (second<first && graph->hub_array[first][second]==1){
                    triangle++;
                }
                //cout<<"triangle number"<<triangle<<"\n";
            }
        }
        //printf("\n");
    }
    cout<<"trangle"<<triangle;
    return triangle;
}
/*int hhnh_traingle(struct Lotus_Graph* graph){
    int triangle=0;
    for (int i = 0; i < N; i++)
    {
        // print current vertex and all its neighbors
        int v=i;
        struct Node* ptr = graph->hub->head[v];
        vector <int> n_v=graph->hubs[i];
        cout<<"\nvertex"<<v;
        while (ptr != NULL)
        {
            //printf("(%d —> %d)\t", i, ptr->dest);
            int neighbor=ptr->dest;
            cout<<"first"<<neighbor;
            vector <int> n_u=graph->hubs[neighbor];
            struct Node* ptr_nu = graph->hub->head[neighbor];
            while(ptr_nu!=NULL){
                int sec_dest=ptr_nu->dest;
                cout<<"second"<<sec_dest;
                if (sec_dest<neighbor){
                    if (graph->hub_array[neighbor][sec_dest]==1){
                        triangle=triangle+1;
                    }
                }
                ptr_nu = ptr_nu->next;
            }
            cout<<"\ntrangle"<<triangle<<" ";
            ptr = ptr->next;
        }
        //printf("\n");
    }
    return triangle;
}*/
int hnn_traingle(struct Lotus_Graph* graph){
    int triangle=0;
    int start=graph->hub_count;
    cout<<"\nhub_count"<<start;
    //cout<<"Inside HNN counting\n";
    //printGraph(graph->non_hub);
    for (int i = start; i < N; i++)
    {
        // print current vertex and all its neighbors
        int v=i;
        struct Node* ptr = graph->non_hub->head[v];
        vector <int> n_v=graph->hubs[i];
        //cout<<"\nvertex"<<v;
        //print_vector(n_v);
        while (ptr != NULL)
        {
            //printf("(%d —> %d)\t", i, ptr->dest);
            int neighbor=ptr->dest;
            vector <int> n_u=graph->hubs[neighbor];
            /*cout<<"\nvertex"<<i;
            print_vector(n_v);
            cout<<"\nneighboor"<<neighbor;
            print_vector(n_u);*/
            vector<int> v(n_v.size() + n_u.size());
            vector<int>::iterator it, st;
            /*cout<<"\nFirst vector";
            print_vector(n_v);
            cout<<"\nSecond vector";
            print_vector(n_u);*/
            sort(n_v.begin(), n_v.end());
            sort(n_u.begin(), n_u.end());
            it = set_intersection(n_v.begin(),n_v.end(),n_u.begin(),n_u.end(),v.begin());
            //cout << "\nCommon elements:\n";
            int count=0;
            for (st = v.begin(); st != it; ++st){
                //cout << *st << ", ";
                count++;
            }
            //cout << '\n';
            triangle=triangle+count;
            //cout<<"\ntrangle"<<triangle<<" ";
            ptr = ptr->next;
        }
        //printf("\n");
    }
    return triangle;
}

int nnn_traingle(struct Lotus_Graph* graph){
    int triangle=0;
    int start=graph->hub_count;
    /*cout<<"\nhub_count"<<start;
    cout<<"Inside HNN counting\n";
    printGraph(graph->non_hub);*/
    for (int i = start; i < N; i++)
    {
        // print current vertex and all its neighbors
        int v=i;
        struct Node* ptr = graph->non_hub->head[v];
        vector <int> n_v=graph->non_hubs[i];
        //cout<<"\nvertex begin"<<v;
        while (ptr != NULL)
        {
            //printf("(%d —> %d)\t", i, ptr->dest);
            int neighbor=ptr->dest;
            vector <int> n_u=graph->non_hubs[neighbor];
            /*cout<<"\nvertex"<<v;
            print_vector(n_v);
            cout<<"\nneighboor"<<neighbor;
            print_vector(n_u);*/
            vector<int> v(n_v.size() + n_u.size());
            vector<int>::iterator it, st;
            /*cout<<"\nFirst vector";
            print_vector(n_v);
            cout<<"\nSecond vector";
            print_vector(n_u);*/
            sort(n_v.begin(), n_v.end());
            sort(n_u.begin(), n_u.end());
            it = set_intersection(n_v.begin(),n_v.end(),n_u.begin(),n_u.end(),v.begin());
            //cout << "\nCommon elements:\n";
            int count=0;
            for (st = v.begin(); st != it; ++st){
                //cout << *st << ", ";
                count++;
            }
            //cout << '\n';
            triangle=triangle+count;
            //cout<<"\ntrangle"<<triangle<<" ";
            ptr = ptr->next;
        }
        //printf("\n");
    }
    return triangle;
}


int triangle_counting(struct Lotus_Graph* graph)
{
    int triangle=0;
    int hhnh=hhnh_traingle(graph);
    cout<<"\nhhnh"<<hhnh;
    int hnn=hnn_traingle(graph);
    cout<<"\nhnn"<<hnn;
    int nnn=nnn_traingle(graph);
    cout<<"\nnnn"<<nnn;
    triangle=triangle+hhnh+hnn+nnn;
    return triangle;
}

// Directed graph implementation in C
int main(void)
{
    // input array containing edges of the graph (as per the above diagram)
    // (x, y) pair in the array represents 
    ifstream input("facebook_combined.txt");
    string line;
    int count=0;

    //struct Edge edges[];
    while(getline( input, line ) ) {
        //cout<<line<<'\n';
        count++;
    }
    ifstream input2("facebook_combined.txt");
    ofstream myfile;
    myfile.open ("example.txt");
    myfile.close();

    struct Edge edges[count];
    count=0;
    while(getline( input2, line ) ) {
        //cout<<line<<'\n';
        string arr[2];
        int i = 0;
        std::stringstream ssin(line);
        while (ssin.good() && i < 2){
            ssin >> arr[i];
            ++i;
        }
        int src=atoi(arr[0].c_str());
        int dest=atoi(arr[1].c_str());
        if (src==0)
            myfile << "Writing this to a file.\n";
        if (src<N && dest<N){   
            edges[count].src=src;
            //cout<<"src"<<edges[count].src<<"print";
            edges[count].dest=dest;
            count++;
        }
    }
    cout<<"input done";
 
    // calculate the total number of edges
    int n = sizeof(edges)/sizeof(edges[0]);
 
    // construct a graph from the given edges
    struct Graph *graph = createGraph(edges, count);
    countEdgeGraph(graph);
 
    // Function to print adjacency list representation of a graph
    //printGraph(graph);

    graph=create_relabeling_array(graph);
    countEdgeGraph(graph);
    ///
    //printGraph(graph);
    cout<<"before: "<<get_degree(graph,2)<<"\n";
    graph=Lotus_Preprocessing(graph);
    countEdgeGraph(graph);
    cout<<"after "<<get_degree(graph,2);
    //printGraph(graph);
    struct Lotus_Graph *graph_lotus=Lotus_Preprocessing_neighbours(graph);
    cout<<"In main\n";
    //printGraph(graph_lotus->hub);
    countEdgeGraph(graph_lotus->hub);
    countEdgeGraph(graph_lotus->non_hub);
    //printarray(graph_lotus->hub_array, graph_lotus->hub_count);
    int triangle=triangle_counting(graph_lotus);
    cout<<"\nFinal count of triangle="<<triangle;
 
    return 0;
}
