#include <stdio.h>
#include <stdlib.h>

#define MAX_CITIES 1001

void DFS(int city, int numVertices, int graph[MAX_CITIES][MAX_CITIES], int visited[])
{
    visited[city] = 1;

    // 若 node 比 neighbor 小
    for (int i = 1; i <= numVertices; i++)
    {
        if (graph[city][i] && !visited[i])
        {
            graph[city][i] = 2;
            DFS(i, numVertices, graph, visited);
        }
    }
    // 若 node 比 neighbor 大
    for (int i = 1; i <= numVertices; i++)
    {
        if (graph[city][i] && !visited[i])
        {
            graph[city][i] = 2;
            DFS(i, numVertices, graph, visited);
        }
    }
}

int main(void)
{
    int neighbor, num_vertices, node;
    char temp;
    scanf("%d", &num_vertices);

    int graph[MAX_CITIES][MAX_CITIES] = {0};
    int visited[MAX_CITIES] = {0};

    for (int i = 0; i < num_vertices; i++)
    {
        int size = 0, node;
        do
        {
            if (scanf("%d", &neighbor) == EOF)
            {
                if (node == neighbor)
                {
                    graph[node][neighbor] = 2;
                }
                else
                {
                    graph[node][neighbor] = 1;
                }

                break;
            }
            if (size == 0)
            {
                node = neighbor;
                size++;
            }
            else
            {
                graph[node][neighbor] = 1;
            }
            temp = getchar();
            if (temp == '\n')
            {
                break;
            }
        } while (temp != '\n');
    }

    // 開始查找
    int minAdditionalRoutes = 0;
    int redundantRoutes = 0;
    int next_start = 1;

    // 看是否有多餘的路段，走過的路就進visited
    for (int i = 1; i <= num_vertices; i++)
    {
        visited[i] = 0;
    }

    for (int i = 1; i < num_vertices + 1; i++)
    {
        DFS(next_start, num_vertices, graph, visited);

        if (!visited[i])
        {
            minAdditionalRoutes++;
            next_start = i;
        }
    }

    for (int i = 1; i < num_vertices + 1; i++)
    {
        for (int j = 1; j < num_vertices + 1; j++)
        {
            if (graph[i][j] && graph[i][j] == 1 && graph[j][i] == 1)
            {
                redundantRoutes = 1;
            }
        }
    }

    printf("%d %d\n", minAdditionalRoutes, redundantRoutes);

    return 0;
}
