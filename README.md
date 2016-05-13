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
Cycle in UnDirected Graph
========================
We can use DFS to detect cycle in an undirected graph
We do a DFS traversal of the given graph. For every visited vertex ‘v’, if there is an adjacent ‘u’ such that u is already visited and u is not parent of v, then there is a cycle in graph.
Graph[][V]
Int visited[V]
 for (int i = 0; i < V; ++i)
   {
       visited[i] = -1;
    }

// check for all trees in a disconnected forest
    for(int i = 0; i < V; i++)
     if(visited[src] != -1)
        if (Cycle(i, -1))
            return true;
Int cycle(int src, int parent)
{
    int temp, u;
    visited[src] = 1;
    for (u = 0; u < V; ++u)
   {
       if (Graph[src][u])
       {
          if (visited3[u] == -1)
             {
                  if(cycle(u,src)
                 return 1;
            }
            else if(u != parent)
            {
                return 1;
            }
       }
    }
}



========================
Topological Sorting
========================
Topological sorting in useful in cases we need to draw dependency or ordering
Topological sorting for Directed Acyclic Graph (DAG) is a linear ordering of vertices such that for every directed edge uv, vertex u comes before v in the ordering
Topological Sorting for a graph is not possible if the graph is not a DAG.
For example, a topological sorting of the following graph is “5 4 2 3 1 0″
DFS can be used for topological sorting
In DFS, we start from a vertex, we first print it and then recursively call DFS for its adjacent vertices. In topological sorting, we use a temporary stack. We don’t print the vertex immediately, we first recursively call topological sorting for all its adjacent vertices, then push it to a stack. Finally, print contents of stack.
Graph[][V]
Int visited[V],stack[V],indegree[V],stack_count;
 for (int i = 0; i < V; ++i)
   {
       visited[i] = -1;
       indegree[i]=0;
       stack = -1;
    }

for (i = 0; i < E; i++)
{
    int x, y;
   scanf("%d", &x);
   scanf("%d", &y);
   Graph[x][y] = 1;
   indegree[y] = indegree[y]+1;
}

// check for all trees in a disconnected forest
    for(int i = 0; i < V; i++)
     if(visited[src] != -1)
        if(indegree[i] == 0)
            topo(i);

void topo(int v)
{
    int i;
    visited[v] = 1;
    for (i = 1; i <= V; i++)
    {
        if (Graph[v][i])
       {
           if (visited[i] != 1)
                topo(i);

         }
     }
      stack[stack_count++] = v;
}

for (i = stack_count; i > 0; i--)
{
     printf(" %d", stack[i-1]);
}

========================
BiPartite
========================
A Bipartite Graph is a graph whose vertices can be divided into two independent sets, U and V such that every edge (u, v) either connects a vertex from U to V or a vertex from V to U. 
We can also say that there is no edge that connects vertices of same set.
It is possible to color a cycle graph with even cycle using two colors
It is not possible to color a cycle graph with odd cycle using two colors.
BFS can be used for checking the BiPartite graph by using different color at each level
Graph[][V]
int isBiPartite(int src)
{
    int temp, u;
    int colorMatrix[V]
    struct Queue* queue = createQueue(V);
    color = 1;
    colorMatrix[src] = color;    
    enQueue(queue, src);
    while (!isEmpty(queue))
    {
        temp = deQueue(queue);
        color = 1 - colorMatrix[temp];        
        for (u = 0; u < V; ++u)
        {
            if (Graph[temp][u])
           {
               if (colorMatrix[u] == -1)
              {
                  colorMatrix[u] = color;
                  enQueue(queue, u);
             }
             else if(colorMatrix[u] == colorMatrix[temp])
              {
                   return 0;
              }
           }
        }
    }
   return 1;
}

========================
Articulation Points (or Cut Vertices) 
========================
A vertex in an undirected connected graph is an articulation point (or cut vertex) iff removing it (and edges through it) disconnects the graph. 
Articulation points represent vulnerabilities in a connected network – single points whose failure would split the network into 2 or more disconnected components.


========================
How to find all articulation points in a given graph
========================
A simple approach is to one by one remove all vertices and see if removal of a vertex causes disconnected graph. Following are steps of simple approach for connected graph.
1) For every vertex v, do following…..a) Remove v from graph..…b) See if the graph remains connected (We can either use BFS or DFS)…..c) Add v back to the graph


========================
Bridges
========================
An edge in an undirected connected graph is a bridge iff removing it disconnects the graph. 
For a disconnected undirected graph, definition is similar, a bridge is an edge removing which increases number of connected components.
Like Articulation Points, bridges represent vulnerabilities in a connected network and are useful for designing reliable networks. 
For example, in a wired computer network, an articulation point indicates the critical computers and a bridge indicates the critical wires or connections.

========================
How to find all bridges in a given graph
========================
A simple approach is to one by one remove all edges and see if removal of a edge causes disconnected graph. Following are steps of simple approach for connected graph.
1) For every edge (u, v), do following…..a) Remove (u, v) from graph..…b) See if the graph remains connected (We can either use BFS or DFS)…..c) Add (u, v) back to the graph.


========================
Check if a given graph is tree or not
========================
An undirected graph is tree if it has following properties.1) There is no cycle.2) The graph is connected.
For an undirected graph we can either use BFS or DFS to detect above two properties.


========================
When to use BFS
========================
Traversal
Bipartite
Shortest path in a unweighted graph
Chess problem (Minimum steps to reach from one point to another)
Snake and Ladder
Time taken to traverse all nodes at different level
Laughing Bomb



========================
Chess
========================
int N, M;
int chess[51][51] = { 0 };
int S, K;
int mv[8][2] = { { -1, -2 }, { -2, -1 }, { -2, 1 }, { -1, 2 }, { 1, -2 }, { 2, -1 }, { 2, 1 }, { 1, 2 } };
int min_count;

int knight_move(int R, int C)
{
    struct Queue* Q;
    struct node Node;
    int i;
    struct node temp = { 0 };
    int f = 0;

    Q = createQueueStruct(N*M);

    Node.x = R;
    Node.y = C;
    Node.z = 0;
    chess[Node.x][Node.y] = 1;
    enQueueStruct(Q, &Node);

    while (!isEmptyStruct(Q))
    {
        if (f)
       {
            break;
       }
       temp = deQueueStruct(Q);
        if (temp.x < 1 || temp.x > N || temp.y < 1 || temp.y > M)
       {
            continue;
       }
       if (temp.x == S && temp.y == K)
       {
           f = 1;
          break;
       }
       for (i = 0; i < 8; i++)
      {
           struct node  temp1;
          temp1.x = temp.x + mv[i][0];
          temp1.y = temp.y + mv[i][1];
          temp1.z = temp.z + 1;
          if (temp1.x < 1 || temp1.x > N || temp1.y < 1 || temp1.y > M)
          {
              continue;
          }
           if (chess[temp1.x][temp1.y] == 0)
          {
               enQueueStruct(Q, &temp1);
               chess[temp1.x][temp1.y] = 1;
           }
       }
    }
    if (f)
        return temp.z;
    else
        return -1;
}

========================
Laughing Bomb
========================
int N, M;
int O, P;
int city[101][101] = { 0 };
int people_count;
int time_count;
int mv[4][2] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };

int laughing_bomb(int O, int P)
{
    struct Queue* Q;
    struct node Node;
    int i;
    struct node temp = { 0 };
    int f = 0;

    Q = createQueueStruct(N*M);

    Node.x = O;
    Node.y = P;
    Node.z = 1;

    if (city[Node.x][Node.y] == 1)
        people_count--;
    city[Node.x][Node.y] = 2; //visited
    enQueueStruct(Q, &Node);
    while (!isEmptyStruct(Q))
    {
         temp = deQueueStruct(Q);
         if (people_count == 0)
         {
              break;
          }
          for (i = 0; i < 4; i++)
          {
                struct node  temp1;
                temp1.x = temp.x + mv[i][0];
                temp1.y = temp.y + mv[i][1];
                temp1.z = temp.z + 1;
                if (temp1.x < 1 || temp1.x > N || temp1.y < 1 || temp1.y > M)
               {
                    continue;
               }

               if (city[temp1.x][temp1.y] != 0 && city[temp1.x][temp1.y] != 2)
               {
                   enQueueStruct(Q, &temp1);
                   city[temp1.x][temp1.y] = 2;
                   people_count--;
                   time_count = temp1.z;
              }
         }
     }
     return time_count;
}

========================
When to Use DFS
========================
Traversal
Find Cycle(Directed and UnDirected)
Topological sort
Find Articulation Point
Find Bridge
Flood Fill
Traverse all nodes with backtracking
Toll Problem

========================
int N, t_cost[25], t_hire[25], min_cost = 1000000;
void dfs(int t_pos, int bpool3, int bpool2, int bpool1, int cur_cost)
{
    int  tot_bpool = bpool3 + bpool2 + bpool1;
    if (cur_cost > min_cost) return;   // condition important to avoid unnecessary cpu cycle
    if (t_pos == N - 1)   //base condition to check last toll gate
    {
        if (tot_bpool  < t_hire[t_pos]) cur_cost += t_cost[t_pos];
        if (cur_cost   < min_cost) min_cost = cur_cost;
            return;
    }
    dfs(t_pos + 1, bpool3, bpool2, bpool1, cur_cost + t_cost[t_pos]);  //toll pay option
    dfs(t_pos + 1, bpool3 + t_hire[t_pos], bpool2, bpool1, cur_cost + 2 * t_cost[t_pos]);  //toll hire option
    if (tot_bpool >= t_hire[t_pos])  //toll battle option
    {
        if (t_hire[t_pos] > bpool2 + bpool1)
        {
            bpool3 = tot_bpool - t_hire[t_pos];
            bpool1 = bpool2 = 0;
       }
       else if (t_hire[t_pos] > bpool1)
       {
           bpool2 = (bpool1 + bpool2) - t_hire[t_pos];
           bpool1 = 0;
       }
       dfs(t_pos + 1, 0, bpool3, bpool2, cur_cost); // note: pool3 is zero, pool3 becomes pool2 and pool2 as pool1
    }
}


========================
Flood Fill
========================
void floodFillUtil(int screen[][N], int x, int y, int prevC, int newC)
{
    // Base cases
    if (x < 0 || x >= M || y < 0 || y >= N)
        return;
    if (screen[x][y] != prevC)
        return;
 
    // Replace the color at (x, y)
    screen[x][y] = newC;
 
    // Recur for north, east, south and west
    floodFillUtil(screen, x+1, y, prevC, newC);
    floodFillUtil(screen, x-1, y, prevC, newC);
    floodFillUtil(screen, x, y+1, prevC, newC);
    floodFillUtil(screen, x, y-1, prevC, newC);
}


========================
Count Stars using Flood Fill
========================
int Const[26][26] = { 0 };
int N1;
int max_itr = 0;
int count1;

for (i = 0; i < N1; i++)
{
    for (j = 0; j < N1; j++)
    {
        fill(i, j);
        if (count1 > 0)
       {
           Num++;
           if (count1 > Max)
               Max = count1;
          count1 = 0;
}
}
}

void fill(int x, int y)
{
    if (x<0 || x>= N1 || y<0 || y >= N1)
    {
         return;
    }
    if (Const[x][y] == 0)
    {
        return;
    }

    Const[x][y] = 0;
    count1++;

    fill(x + 1, y);
    fill(x - 1, y);
    fill(x, y + 1);
    fill(x, y - 1);
}

========================
Single Source Shortest Path- Dijkstra’s(positive weight’s)
========================
// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
void dijkstra(int graph[V][V], int src)
{
     int dist[V];     // The output array.  dist[i] will hold the shortest
    // distance from src to i
    bool sptSet[V]; // sptSet[i] will true if vertex i is included in shortest
    // path tree or shortest distance from src to i is finalized
   // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < V; i++)
   dist[i] = INT_MAX, sptSet[i] = false;
    // Distance of source vertex from itself is always 0
    dist[src] = 0;
    // Find shortest path for all vertices
    for (int count = 0; count < V - 1; count++)
   {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in first iteration.
        int u = minDistance(dist, sptSet);
        // Mark the picked vertex as processed
        sptSet[u] = true;
        // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < V; v++)
        // Update dist[v] only if is not in sptSet, there is an edge from 
       // u to v, and total weight of path from src to  v through u is 
       // smaller than current value of dist[v]
       if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
       && dist[u] + graph[u][v] < dist[v])
       dist[v] = dist[u] + graph[u][v];
    }
    // print the constructed distance array
    printSolution(dist, V);
}
int minDistance(int dist[], bool sptSet[])
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
    if (sptSet[v] == false && dist[v] <= min)
    min = dist[v], min_index = v;

    return min_index;
}

========================
Dijkstra for Matrix traversal
========================
void getMinDistNode(int*ci, int*cj, int n)
{
    int minDist = 999999;
    int i, j;
    for (i = 0; i < n; i++)
   {
       for (j = 0; j < n; j++)
      {
          if (!sptSet[i][j] && minDist>dist[i][j])
          {
               minDist = dist[i][j];
               *ci = i;
              *cj = j;
          }
       }
   }
}
int getShortestPath(int n)
{
    int i, j;    int currI, currJ;
   for (i = 0; i < n; i++)
   {
       for (j = 0; j < n; j++)
      {
          sptSet[i][j] = false; dist[i][j] = 999999;
       }
   }
    dist[0][0] = 0;
 for (i = 0; i < n; i++)
   {
       for (j = 0; j < n; j++)
      {
           getMinDistNode(&currI, &currJ, n);
           sptSet[currI][currJ] = true;
            if (currI + 1 < n)
           {
                if (!sptSet[currI + 1][currJ] && dist[currI + 1][currJ] > dist[currI][currJ] + array[currI + 1][currJ])
               {
                    dist[currI + 1][currJ] = dist[currI][currJ] + array[currI + 1][currJ];
               }
          }
          if (currI - 1 >= 0)
         {
              if (!sptSet[currI - 1][currJ] && dist[currI - 1][currJ] > dist[currI][currJ] + array[currI - 1][currJ])
             {
                 dist[currI - 1][currJ] = dist[currI][currJ] + array[currI - 1][currJ];
             }
         }
         if (currJ + 1 < n)
         {
              if (!sptSet[currI][currJ + 1] && dist[currI][currJ + 1] > dist[currI][currJ] + array[currI][currJ + 1])
             {
                 dist[currI][currJ + 1] = dist[currI][currJ] + array[currI][currJ + 1];
             }
        }
        if (currJ - 1 >= 0)
        {
             if (!sptSet[currI][currJ - 1] && dist[currI][currJ - 1] > dist[currI][currJ] + array[currI][currJ - 1])
            {
                 dist[currI][currJ - 1] = dist[currI][currJ] + array[currI][currJ - 1];
             }
         }
       }
   }
   return dist[n - 1][n - 1];
}

========================
Single Source Shortest Path- Bellman-Ford(negative weights weight’s)
========================
// The main function that finds shortest distances from src to
// all other vertices using Bellman-Ford algorithm.  The function
// also detects negative weight cycle
void BellmanFord(struct Graph* graph, int src)
{
    int V = graph->V;
    int E = graph->E;
    int dist[V];
 
    // Step 1: Initialize distances from src to all other vertices
    // as INFINITE
    for (int i = 0; i < V; i++)
        dist[i]   = INT_MAX;
    dist[src] = 0;
 
    // Step 2: Relax all edges |V| - 1 times. A simple shortest 
    // path from src to any other vertex can have at-most |V| - 1 
    // edges
    for (int i = 1; i <= V-1; i++)
    {
        for (int j = 0; j < E; j++)
        {
            int u = graph->edge[j].src;
            int v = graph->edge[j].dest;
            int weight = graph->edge[j].weight;
            if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
                dist[v] = dist[u] + weight;
        }
    }
 
 
    printArr(dist, V);
 
    return;
}
    // Step 3: check for negative-weight cycles.  The above step 
    // guarantees shortest distances if graph doesn't contain 
    // negative weight cycle.  If we get a shorter path, then there
    // is a cycle.
    for (int i = 0; i < E; i++)
    {
        int u = graph->edge[i].src;
        int v = graph->edge[i].dest;
        int weight = graph->edge[i].weight;
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
            printf("Graph contains negative weight cycle");
    }

========================
All-PairShortest Path- Floyd Warshall
========================
// Solves the all-pairs shortest path problem using Floyd Warshall algorithm
void floydWarshell (int graph[][V])
{
    /* dist[][] will be the output matrix that will finally have the shortest 
      distances between every pair of vertices */
    int dist[V][V], i, j, k;
 
    /* Initialize the solution matrix same as input graph matrix. Or 
       we can say the initial values of shortest distances are based
       on shortest paths considering no intermediate vertex. */
    for (i = 0; i < V; i++)
        for (j = 0; j < V; j++)
            dist[i][j] = graph[i][j];
 
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of a iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of a iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */
    for (k = 0; k < V; k++)
    {
        // Pick all vertices as source one by one
        for (i = 0; i < V; i++)
        {
            // Pick all vertices as destination for the
            // above picked source
            for (j = 0; j < V; j++)
            {
                // If vertex k is on the shortest path from
                // i to j, then update the value of dist[i][j]
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }
 
    // Print the shortest distance matrix
    printSolution(dist);
}
 
    printArr(dist, V);
 
    return;
}

========================
Shortest Path – Directed Acyclic
========================
- Create a toplogical order of all vertices.- Do following for every vertex u in topological order.………..Do following for every adjacent vertex v of u………………if (dist[v] > dist[u] + weight(u, v))………………………dist[v] = dist[u] + weight(u, v) 
// The function to find shortest paths from given vertex. It uses recursive 
// topologicalSortUtil() to get topological sorting of given graph.
void Graph::shortestPath(int s)
{
    stack<int> Stack;
    int dist[V];
 
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
 
    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);
 
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++)
        dist[i] = INF;
    dist[s] = 0;
 
    
// Process vertices in topological order
    while (Stack.empty() == false)
    {
        // Get the next vertex from topological order
        int u = Stack.top();
        Stack.pop();
 
        // Update distances of all adjacent vertices
        list<AdjListNode>::iterator i;
        if (dist[u] != INF)
        {
          for (i = adj[u].begin(); i != adj[u].end(); ++i)
             if (dist[i->getV()] > dist[u] + i->getWeight())
                dist[i->getV()] = dist[u] + i->getWeight();
        }
    }
 
    // Print the calculated shortest distances
    for (int i = 0; i < V; i++)
        (dist[i] == INF)? cout << "INF ": cout << dist[i] << " ";
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
