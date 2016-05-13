# leetcode
========================
Searching and Sorting
========================
Search
Linear Search (Binary search is better if u can practice)
Take care of boundary Condition 0/1, >=<, N/N+1
for (i = 1; i <= N ; i++)
{
    if (array[i] == S)
    break;
}
Sorting
Bubble sort (Quick sort is better if u can practice)
for (i = 1; i < N; i++)
{
    for (j = i + 1; j <= N; j++)
    {
        if (array[i]>array[j])
       {
           int temp = array[i];
           array[i] = array[j];
           array[j] = temp;
       }
    }
}

========================
Queue
========================
Why Queue
Q is mandatory for BFS 
#include <stdio.h>
#include <malloc.h>
struct node
{
int x,y,z;
};
struct Queue
{
    int front, rear, capacity;
    struct node* array;
};
struct Queue* createQueueStruct(int capacity)
{
    struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
    queue->front = -1;
    queue->rear = -1;
    queue->capacity = capacity;
    queue->array = (struct node*)malloc(capacity*(sizeof(struct node)));
    return queue;
}
int isFull(struct Queue* queue)
{
    if (queue->rear == queue->capacity - 1)
        return 1;
    else
        return 0;
}



int isEmpty(struct Queue* queue)
{
    if (queue->front == -1)
        return 1;
    else
        return 0;
}
void enQueue(struct Queue* queue, struct node* item)
{
    if (isFullStruct(queue))
        return;
    ++queue->rear;
    queue->array[queue->rear].x = item->x;
    queue->array[queue->rear].y = item->y;
    queue->array[queue->rear].z = item->z;
    if (queue->front == -1)
        ++queue->front;
}
struct node deQueue(struct Queue* queue)
{
    struct node item = { 0 };
    if (isEmptyStruct(queue))
        return item;
    item.x = queue->array[queue->front].x;
    item.y = queue->array[queue->front].y;
    item.z = queue->array[queue->front].z;
    if (queue->front == queue->rear)
        queue->front = queue->rear = -1;
    else
        ++queue->front;
    return item;
}

========================
Graph Representation
========================
Adjacency Matrix
Graph[101][101] = {0}
//Always give size more to avoid any segmentation fault, Always initialize or reset for every test case
scanf("%d %d", &X,&Y);
Graph[X][Y] = 1




Adjacency List
// a structure to represent a weighted edge in graph
struct Edge
{
int src, dest, weight;
};

// a structure to represent a connected, undirected and weighted graph
struct Graph
{
// V-> Number of vertices, E-> Number of edges
int V, E;

// graph is represented as an array of edges. Since the graph is
// undirected, the edge from src to dest is also edge from dest
// to src. Both are counted as 1 edge here.
struct Edge* edge;
};

// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
graph->V = V;
graph->E = E;

graph->edge = (struct Edge*) malloc(graph->E * sizeof(struct Edge));

return graph;
}



struct Graph* graph = createGraph(V, E);


// add edge 0-1
graph->edge[0].src = 0;
graph->edge[0].dest = 1;
graph->edge[0].weight = 10;

========================
BFS/DFS
========================
Graph[][V]
Int visited[V],reach[V],count;
 for (int i = 0; i < V; ++i)
   {
       visited[i] = -1;
       reach[i] = -1;
    }
void bfs(int src)
{
    int temp, u;
struct Queue* queue = createQueue(V);
    visited[src] = 1;
    enQueue(queue, src);
    while (!isEmpty(queue))
    {
        temp = deQueue(queue);
        reach[count++] = temp;  // Store the visited node
        for (u = 0; u < V; ++u)
        {
            if (Graph[temp][u])
           {
               if (visited[u] == -1)
              {
                   visited[u] = 1;
                  enQueue(queue, u);
             }
           }
        }
    }
}
Graph[][V]
Int visited[V],reach[V],count;
 for (int i = 0; i < V; ++i)
   {
       visited[i] = -1;
       reach[i] = -1;
    }
void dfs(int src)
{
    int temp, u;
    visited[src] = 1;
    reach[count++] = src;
    for (u = 0; u < V; ++u)
    {
        if (Graph[src][u])
       {
           if (visited[u] == -1)
          {
              dfs(u);
         }
       }
    }
}

========================
Cycle in Directed Graph
========================
DFS can be used to detect cycle in a Graph
There is a cycle in a graph only if there is a back edge present in the graph.
A back edge is an edge that is from a node to itself (selfloop) or one of its ancestor in the tree produced by DFS
To detect a back edge, we can keep track of vertices currently in recursion stack of function for DFS traversal
If we reach a vertex that is already in the recursion stack, then there is a cycle in the tree
Graph[][V]
Int visited[V], recStack[V];
 for (int i = 0; i < V; ++i)
   {
       visited[i] = -1;
       recStack[i] = -1;
    }

// check for all trees in a disconnected forest
    for(int i = 0; i < V; i++)
        if (Cycle(i, recStack))
            return true;

Int cycle(int src, int* recStack)
{
    int temp, u;
    if(visited[src] != -1)
   {
        visited[src] = 1;
        recStack[src] = 1;
        for (u = 0; u < V; ++u)
       {
           if (Graph[src][u])
          {
              if (visited[u] == -1 && cycle(u,recStack)
             {
                 return 1;
            }
            else if(recStack[u] != -1)
            {
                return 1;
            }
       }
    }
    recStack[src] = 0;
    return 0;
}
}


========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
========================
