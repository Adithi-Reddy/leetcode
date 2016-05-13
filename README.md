
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
Shortest path with exactly k edges in a directed and weighted graph
========================
// A naive recursive function to count walks from u to v with k edges
int shortestPath(int graph[][V], int u, int v, int k)
{
   // Base cases
   if (k == 0 && u == v)             return 0;
   if (k == 1 && graph[u][v] != INF) return graph[u][v];
   if (k <= 0)                       return INF;
 
   // Initialize result
   int res = INF;
 
   // Go to all adjacents of u and recur
   for (int i = 0; i < V; i++)
   {
       if (graph[u][i] != INF && u != i && v != i)
       {
           int rec_res = shortestPath(graph, i, v, k-1);
           if (rec_res != INF)
              res = min(res, graph[u][i] + rec_res);
       }
   }
   return res;
}


========================
Minimum Spanning Tree
========================
Minimum Spanning Tree (MST) problem: Given connected graph G with positive edge weights, find a min weight set of edges that connects all of the vertices.
Network design.– telephone, electrical, hydraulic, TV cable, computer, road

========================
Minimum Spanning Tree (Kruskal’s)
========================
1. Sort all the edges in non-decreasing order of their weight. 
2. Pick the smallest edge. Check if it forms a cycle with the spanning tree formed so far. If cycle is not formed, include this edge. Else, discard it. 
3. Repeat step#2 until there are (V-1) edges in the spanning tree.
4. // A structure to represent a subset for union-find
struct subset
{
    int parent;
    int rank;
};
// A utility function to find set of an element i
// (uses path compression technique)
int find(struct subset subsets[], int i)
{
    // find root and make root as parent of i (path compression)
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);
    return subsets[i].parent;
}
// A function that does union of two sets of x and y
// (uses union by rank)
void Union(struct subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);
    // Attach smaller rank tree under root of high rank tree
    // (Union by Rank)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;
    // If ranks are same, then make one as root and increment
    // its rank by one
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

void KruskalMST(struct Graph* graph)
{
    int V = graph->V;
    struct Edge result[V];  // Tnis will store the resultant MST
    int e = 0;  // An index variable, used for result[]
    int i = 0;  // An index variable, used for sorted edges
    // Step 1:  Sort all the edges in non-decreasing order of their weight
    // If we are not allowed to change the given graph, we can create a copy of
    // array of edges
    qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);
    // Allocate memory for creating V ssubsets
    struct subset *subsets =
        (struct subset*) malloc( V * sizeof(struct subset) );
     // Create V subsets with single elements
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v; subsets[v].rank = 0;
    }
    // Number of edges to be taken is equal to V-1
    while (e < V - 1)
    {
        // Step 2: Pick the smallest edge. And increment the index
        // for next iteration
        struct Edge next_edge = graph->edge[i++];
        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);
        // If including this edge does't cause cycle, include it
        // in result and increment the index of result for next edge
        if (x != y)
        {
            result[e++] = next_edge;
            Union(subsets, x, y);
        }
        // Else discard the next_edge
    }
    // print the contents of result[] to display the built MST
    return;
}


========================
Minimum Spanning Tree (Prim’s)
========================
void primMST(int graph[V][V])
{
     int parent[V]; // Array to store constructed MST
     int key[V];   // Key values used to pick minimum weight edge in cut
     bool mstSet[V];  // To represent set of vertices not yet included in MST
     // Initialize all keys as INFINITE
     for (int i = 0; i < V; i++)
        key[i] = INT_MAX, mstSet[i] = false;
     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST 
     // The MST will have V vertices
     for (int count = 0; count < V-1; count++)
     {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < V; v++)
           // graph[u][v] is non zero only for adjacent vertices of m
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if graph[u][v] is smaller than key[v]
          if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v])
             parent[v]  = u, key[v] = graph[u][v];
     }
     printMST(parent, V, graph);      // print the constructed MST
}
// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int minKey(int key[], bool mstSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < V; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
 
   return min_index;
}
 
// A utility function to print the constructed MST stored in parent[]
int printMST(int parent[], int n, int graph[V][V])
{
   printf("Edge   Weight\n");
   for (int i = 1; i < V; i++)
      printf("%d - %d    %d \n", parent[i], i, graph[i][parent[i]]);
}


========================
Dynamic Programming
========================
Dynamic Programming is an algorithmic paradigm that solves a given complex problem by breaking it into subproblems and stores the results of subproblems to avoid computing the same results again. 
Overlapping Sub Problem
In dynamic programming, computed solutions to subproblems are stored in a table so that these don’t have to recomputed. 
There are following two different ways to store the values so that these values can be reused. 
	Memoization (Top Down): 	Tabulation (Bottom Up): 
Optimal Substructure: 
A given problems has Optimal Substructure Property if optimal solution of the given problem can be obtained by using optimal solutions of its subproblems.


========================
Cut Rod
========================
/* Returns the best obtainable price for a rod of length n and
   price[] as prices of different pieces */
int cutRod(int price[], int n)
{
   int val[n+1];
   val[0] = 0;
   int i, j;
 
   // Build the table val[] in bottom up manner and return the last entry
   // from the table
   for (i = 1; i<=n; i++)
   {
       int max_val = INT_MIN;
       for (j = 0; j < i; j++)
         max_val = max(max_val, price[j] + val[i-j-1]);
       val[i] = max_val;
   }
 
   return val[n];
}


========================
Bus Fare Problem based on Cut Rod
========================
int minFare(int fare[], int n, int distance) {
    int knapsack[9999];
    knapsack[0] = 0;
    int i, j;
    for (j = 1; j <= distance; j++)
   {
        int min = 1000000;
        for (i = 0; i<n && i<distance; i++) 
       {
            int x = j - (i + 1);
            if (x >= 0 && (knapsack[x] + fare[i]) < min) 
           {
               min = knapsack[x] + fare[i];
          }
       }
       knapsack[j] = min;
    }
    return knapsack[distance];
}

========================
Stairs Problem
========================
max_array[0] = 0;
max_array[1] = array[1];
max_array[2] = array[1] + array[2];
for (i = 3; i <= N; i++)
{
    int x = max_array[i - 3] + array[i - 1] + array[i];
    int y = max_array[i - 2] + array[i];
   if (x>y)
       max_array[i] = x;
   else
      max_array[i] = y;
}
Answer = max_array[N];

========================
Greedy
========================
Greedy is an algorithmic paradigm that builds up a solution piece by piece, always choosing the next piece that offers the most obvious and immediate benefit. 
Greedy algorithms are used for optimization problems. An optimization problem can be solved using Greedy if the problem has the following property: 
At every step, we can make a choice that looks best at the moment, and we get the optimal solution of the complete problem.
Examples
Kruskal’s Minimum Spanning Tree (MST): 
 Prim’s Minimum Spanning Tree
Dijkstra’s Shortest Path

========================
Binary to Decimal
========================
int binary_decimal(int* binary, int size)
{
    int dec = binary[0];

    for (int i = 0; i < size-1; i++)
   {
       //Double and Add
       dec = (dec * 2) + binary[i + 1];
   }
    return dec;
}

int binary_decimal(int num)
{
    int decimal_val = 0, base =1, rem=0;
    while (num > 0)
   {
       rem = num % 10;
       decimal_val = decimal_val + rem * base;
       num = num / 10;
       base = base * 2;
   }
   return decimal_val;
}

========================
Power
========================
/* Function to calculate x raised to the power y */
int power(int x, unsigned int y)
{
    if (y == 0)
        return 1;
    else if (y % 2 == 0)
        return power(x, y / 2)*power(x, y / 2);
     else
        return x*power(x, y / 2)*power(x, y / 2);
}

========================
GCD
========================
int gcd_it(int a, int b)
{
    int temp;
    while (b > 0)
    {
        temp = a%b;
        a = b;
        b = temp;
    }
    return a;
}

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b%a, a);
}


========================
Fibonacci
========================
int fib(int n)
{
    if (n <= 1)
        return n;
    return fib(n - 1) + fib(n - 2);
}

int fib(int n)
{
    /* Declare an array to store Fibonacci numbers. */
    int f[MAX + 1];
    int i;

    /* 0th and 1st number of the series are 0 and 1*/
    f[0] = 0;
    f[1] = 1;

    for (i = 2; i <= n; i++)
    {
        /* Add the previous 2 numbers in the series and store it */
        f[i] = f[i - 1] + f[i - 2];
    }
    return f[n];
}

========================
Min-Max
========================
struct pair
{
    int min;
    int max;
};

struct pair getMinMax(int arr[], int low, int high)
{
    struct pair minmax, mml, mmr;
    int mid;

   /* If there is only one element */
   if (low == high)
   {
       minmax.max = arr[low];
       minmax.min = arr[low];
       return minmax;
   }


   /* If there are two elements */
   if (high == low + 1)
   {
        if (arr[low] > arr[high])
       {
           minmax.max = arr[low];
           minmax.min = arr[high];
       }
       else
      {
           minmax.max = arr[high];
           minmax.min = arr[low];
      }
      return minmax;
   }
    /* If there are more than 2 elements */
    mid = (low + high) / 2;
    mml = getMinMax(arr, low, mid);
    mmr = getMinMax(arr, mid + 1, high);

   /* compare minimums of two parts*/
   if (mml.min < mmr.min)
       minmax.min = mml.min;
   else
       minmax.min = mmr.min;

    /* compare maximums of two parts*/
   if (mml.max > mmr.max)
       minmax.max = mml.max;
   else
       minmax.max = mmr.max;

    return minmax;
}

========================
Sub-Matrix
========================
1) Construct a sum matrix S[R][C] for the given M[R][C].
     a)	Copy first row and first columns as it is from M[][] to S[][]
     b)	For other entries, use following expressions to construct S[][]
         If M[i][j] is 1 then
            S[i][j] = min(S[i][j-1], S[i-1][j], S[i-1][j-1]) + 1
         Else /*If M[i][j] is 0*/
            S[i][j] = 0
2) Find the maximum entry in S[R][C]
3) Using the value and coordinates of maximum entry in S[i], print 
   sub-matrix of M[][]

void printMaxSubSquare(bool M[R][C])
{
  int i,j;
  int S[R][C];
  int max_of_s, max_i, max_j; 
  
  /* Set first column of S[][]*/
  for(i = 0; i < R; i++)
     S[i][0] = M[i][0];
  
  /* Set first row of S[][]*/    
  for(j = 0; j < C; j++)
     S[0][j] = M[0][j];


      
  /* Construct other entries of S[][]*/
  for(i = 1; i < R; i++)
  {
    for(j = 1; j < C; j++)
    {
      if(M[i][j] == 1) 
        S[i][j] = min(S[i][j-1], S[i-1][j], S[i-1][j-1]) + 1;
      else
        S[i][j] = 0;
    }    
  } 
  /* Find the maximum entry, and indexes of maximum entry 
     in S[][] */
  max_of_s = S[0][0]; max_i = 0; max_j = 0;
  for(i = 0; i < R; i++)
  {
    for(j = 0; j < C; j++)
    {
      if(max_of_s < S[i][j])
      {
         max_of_s = S[i][j];         max_i = i;          max_j = j;
      }        
    }                 
  }     
  printf("\n Maximum size sub-matrix is: \n");
  for(i = max_i; i > max_i - max_of_s; i--)
  {
    for(j = max_j; j > max_j - max_of_s; j--)
    {
      printf("%d ", M[i][j]);
    }  
  }  
}   

========================
Taekwondo
========================
for (i = 0; i < G1; i++) //small group
{
    for (j = 0; j < G2 - G1 + 1; j++) //big group
    {
        float diff;
        if (group1[i] > group2[i + j])
            diff = group1[i] - group2[i + j];
        else
            diff = group2[i + j] - group1[i];
        value[j] = value[j] + diff;
        if (j>0 && value[j]>value[j - 1])
             value[j] = value[j - 1];
    }
}

Answer = value[G2 - G1];

========================
Chocolate
========================
int N;
long long int interns[10001];
long long int i, j;
long long int min = 2000000000;
scanf("%d", &N);
for (i = 0; i < N; i++)
{
    scanf("%d", &interns[i]);
    if (interns[i] < min)
        min = interns[i];
}
Answer = 2000000000;
for (i = 0; i < 5; i++)
{
    long long int temp;
    temp = solve(min - i, N, interns);
    if (temp < Answer)
        Answer = temp;
}
long long int solve(long long int a, long long int N, long long int* array)
{
    long long int i,ans=0;
    for (i = 0; i < N; i++)
    {
        long long int temp = array[i] - a;
       //ans = ans + ((temp / 5) + ((temp % 5) / 2) + ((temp % 5) % 2));
         if (temp >= 5)
        {
            ans += temp / 5;
            temp = temp % 5;
        }
        if (temp >= 2)
        {
            ans += temp / 2;
           temp = temp % 2;
       }
       ans += temp;
    }
    return ans;
}

========================
Coin Change (Dynamic)
========================
int count( int S[], int m, int n )
{
    // table[i] will be storing the number of solutions for
    // value i. We need n+1 rows as the table is consturcted
    // in bottom up manner using the base case (n = 0)
    int table[n+1];
 
    // Initialize all table values as 0
    memset(table, 0, sizeof(table));
 
    // Base case (If given value is 0)
    table[0] = 1;
 
    // Pick all coins one by one and update the table[] values
    // after the index greater than or equal to the value of the
    // picked coin
    for(int i=0; i<m; i++)
        for(int j=S[i]; j<=n; j++)
            table[j] += table[j-S[i]];
 
    return table[n];
}

// m is size of coins array (number of different coins)
int minCoins(int coins[], int m, int V)
{
    // table[i] will be storing the minimum number of coins
    // required for i value.  So table[V] will have result
    int table[V+1];
 
    // Base case (If given value V is 0)
    table[0] = 0;
 
    // Initialize all table values as Infinite
    for (int i=1; i<=V; i++)
        table[i] = INT_MAX;
 
    // Compute minimum coins required for all
    // values from 1 to V
    for (int i=1; i<=V; i++)
    {
        // Go through all coins smaller than i
        for (int j=0; j<m; j++)
          if (coins[j] <= i)
          {
              int sub_res = table[i-coins[j]];
              if (sub_res != INT_MAX && sub_res + 1 < table[i])
                  table[i] = sub_res + 1;
          }
    }
    return table[V];
}


========================
Knapsack (0-1)
========================
// Returns the maximum value that can be put in a knapsack of capacity W
int knapSack(int W, int wt[], int val[], int n)
{
   int i, w;
   int K[n+1][W+1];
 
   // Build table K[][] in bottom up manner
   for (i = 0; i <= n; i++)
   {
       for (w = 0; w <= W; w++)
       {
           if (i==0 || w==0)
               K[i][w] = 0;
           else if (wt[i-1] <= w)
                 K[i][w] = max(val[i-1] + K[i-1][w-wt[i-1]],  K[i-1][w]);
           else
                 K[i][w] = K[i-1][w];
       }
   }
 
   return K[n][W];
}


========================
Knapsack – Fractional (Greedy)
========================
// Main greedy function to solve problem
double fractionalKnapsack(int W, struct Item arr[], int n)
{
    //    sorting Item on basis of ration
    sort(arr, arr + n, cmp);
    int curWeight = 0;  // Current weight in knapsack
    double finalvalue = 0.0; // Result (value in Knapsack)
 
    // Looping through all Items
    for (int i = 0; i < n; i++)
    {
        // If adding Item won't overflow, add it completely
        if (curWeight + arr[i].weight <= W)
        {
            curWeight += arr[i].weight;
            finalvalue += arr[i].value;
        }
        // If we can't add current Item, add fractional part of it
        else
        {
            int remain = W - curWeight;
            finalvalue += arr[i].value * ((double) remain / arr[i].weight);
            break;
        }
    }
    // Returning final value
    return finalvalue;
}


========================
Celebrity
========================
bool HaveAcquiantance(int a, int b) { return MATRIX[a][b]; }
int Celebrity[4];
int findCelebrity(int size)
{
    /* Initialize all as celebrity */
    for (int i = 0; i < size; i++)
    Celebrity[i] = 1;

    for (int i = 0; i < size; i++)
    {
        if (Celebrity[i] == 1)
        {
            for (int j = i+1; j < size; j++)
            {
                 if (HaveAcquiantance(i, j))
                 {
                      Celebrity[i] = 0;
                      break;
                 }
                  else
                 {
                      Celebrity[j] = 0;
                 }
              }
          }
     }
    int potential_celeb=-1;
    for (int i = 0; i < size; i++)
    {
        if (Celebrity[i] == 1)
        {
             potential_celeb = i;
              break;
          }
     }
     if (potential_celeb == -1)
         return -1;

     for (int i = 0; i < size; i++)
    {
         if (i == potential_celeb)
             continue;
         if (HaveAcquiantance(potential_celeb, i))
            return -1;
         if (!HaveAcquiantance(i, potential_celeb))
             return -1;
    }
    return potential_celeb;
}

========================
Subset Sum
========================
// inputs // s            - set vector// t            - tuplet vector// s_size       - set size
// t_size       - tuplet size so far// sum          - sum so far// ite          - nodes count  // target_sum   - sum to be found
void subset_sum(int s[], int t[],int s_size, int t_size,int sum, int ite,int const target_sum)
{
    total_nodes++;
    if (target_sum == sum)
    {
        printSubset(t, t_size); // We found sum
        if (ite + 1 < s_size && sum - s[ite] + s[ite + 1] <= target_sum)
        {
            // Exclude previous added item and consider next candidate
            subset_sum(s, t, s_size, t_size - 1, sum - s[ite], ite + 1, target_sum);
        }
        return;
    }
    else
    {
        if (ite < s_size && sum + s[ite] <= target_sum)
        {
             // generate nodes along the breadth
            for (int i = ite; i < s_size; i++)
           { 
               t[t_size] = s[i];
               if (sum + s[i] <= target_sum)
              {
                  // consider next level node (along depth)
                  subset_sum(s, t, s_size, t_size + 1, sum + s[i], i + 1, target_sum);
              }
           }
       }
   }
}

========================
Tower of Hanoi
========================
void tower(int n, char sourcePole, char destinationPole, char auxiliaryPole)
{
    // Base case (termination condition)
    if (0 == n)
        return;

    // Move first n-1 disks from source pole
    // to auxiliary pole using destination as
    // temporary pole
    tower(n - 1, sourcePole, auxiliaryPole,destinationPole);

    // Move the remaining disk from source
    // pole to destination pole
    printf("Move the disk %d from %c to %c\n", n,sourcePole, destinationPole);

    // Move the n-1 disks from auxiliary (now source)
    // pole to destination pole using source pole as
    // temporary (auxiliary) pole
    tower(n - 1, auxiliaryPole, destinationPole, sourcePole);
}

========================
Union- Find Algorithm
========================
// A utility function to find the subset of an element i
int find(int parent[], int i)
{
    if (parent[i] == -1)
        return i;
    return find(parent, parent[i]);
}
 
// A utility function to do union of two subsets 
void Union(int parent[], int x, int y)
{
    int xset = find(parent, x);
    int yset = find(parent, y);
    parent[xset] = yset;
}
 

// The main function to check whether a given graph contains 
// cycle or not
int isCycle( struct Graph* graph )
{
    // Allocate memory for creating V subsets
    int *parent = (int*) malloc( graph->V * sizeof(int) );
 
    // Initialize all subsets as single element sets
    memset(parent, -1, sizeof(int) * graph->V);
 
    // Iterate through all edges of graph, find subset of both
    // vertices of every edge, if both subsets are same, then 
    // there is cycle in graph.
    for(int i = 0; i < graph->E; ++i)
    {
        int x = find(parent, graph->edge[i].src);
        int y = find(parent, graph->edge[i].dest);
 
        if (x == y)
            return 1;
 
        Union(parent, x, y);
    }
    return 0;
}

========================
Quick sort
========================
int partition(int* array, int p, int r)
{
    int pivot = array[p];
    int i = p;
    int j;
    for (j = p+1 ; j <= r; j++)
   {
        if (array[j] <= pivot)
       {
           i = i + 1;
           swap(&array[i],&array[j]);
      }
   }
   swap(&array[i], &array[p]);
   return i;
}

int quick_sort(int* array, int p, int r)
{
    if (p < r)
    {
        int q = partition(array, p, r);
        quick_sort(array, p, q - 1);
        quick_sort(array, q + 1, r);
    }
    return 1;
}

========================
Back Tracking
========================
int N;
int Ballon[10] = { 0 };
int Visited[10] = { 0 };
int data[10] = { 0 };
int count;
void dfs_backtracking()
{
    int i = 0;
    for (i = 0; i < N; i++)
    {
        if (!Visited[i])
       {
           Visited[i] = 1;
           data[count++] = Ballon[i];
           if (count == N)
          {
                Visited[i] = 0;
                for (int j = 0; j < N; j++)
                printf("%d ", data[j]);
                printf("\n");
                data[count--] = 0;
                return;
            }
            dfs_backtracking();
            Visited[i] = 0;
            data[count--] = 0;
         }
    }
}

========================
Balloon Shot
========================
int max_cost;
int N;
int Ballon[10] = { 0 };
int Visited[10] = { 0 };
int data[10] = { 0 };
int count;
int get_left(int src)
{
    while (src>0)
    {
        src--;
        if (!Visited[src])
            return src;
    }
    return -1;
}

int get_right(int src)
{
    while (src<N-1)
    {
         src++;
         if (!Visited[src])
             return src;
    } 
    return -1;
}


void ballon_shot(int cur_cost)
{
    int i = 0;
    for (i = 0; i < N; i++)
    {
         if (!Visited[i])
        {
             int left, right;
             Visited[i] = 1;
             data[count++] = Ballon[i];
             left = get_left(i);
             right = get_right(i);
             if (left == -1 && right == -1)
            {
                 if (cur_cost + Ballon[i] > max_cost)
                     max_cost = cur_cost + Ballon[i];
                 Visited[i] = 0;
                 for (int j = 0; j < N;j++)
                     printf("%d ", data[j]);
                 printf("\n%d\n", cur_cost + Ballon[i]);
                 data[count--] = 0;
                return;
            }
            else if (left == -1 && right != -1)
            {
                ballon_shot(cur_cost + Ballon[right]);
            }
            else if (left != -1 && right == -1)
            {
                ballon_shot(cur_cost + Ballon[left]);
            }
            else
            {
                 ballon_shot(cur_cost + Ballon[left] * Ballon[right]);
            }
            Visited[i] = 0;
            data[count--] = 0;
        }     
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
