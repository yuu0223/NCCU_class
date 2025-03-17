#include <stdio.h>
#include <stdlib.h>

int matrix_x, matrix_y;

int dfs(int **matrix, int i, int j, int target)
{
    if (i < 0 || i >= matrix_x || j < 0 || j >= matrix_y || matrix[i][j] != target)
    {
        return 0;
    }

    matrix[i][j] = 2;

    return 1 + dfs(matrix, i - 1, j, target) + dfs(matrix, i + 1, j, target) + dfs(matrix, i, j - 1, target) + dfs(matrix, i, j + 1, target);
}

int find_max_area(int **matrix, int target)
{
    // search area
    int final_max = 0;
    for (int x = 0; x < matrix_x; x++)
    {
        for (int y = 0; y < matrix_y; y++)
        {
            if (matrix[x][y] == target)
            {
                int temp_max = dfs(matrix, x, y, target);

                if (temp_max > final_max)
                {
                    final_max = temp_max;
                }
            }
        }
    }
    return final_max;
}

int main(void)
{
    for (int i = 0; i < 2; i++)
    {
        if (i == 0)
        {
            scanf("%d", &matrix_x);
        }
        else
        {
            scanf("%d", &matrix_y);
        }
    }

    // matrix size
    int **matrix;
    matrix = (int **)malloc(matrix_x * sizeof(int *));
    for (int i = 0; i < matrix_x; i++)
    {
        matrix[i] = (int *)malloc(matrix_y * sizeof(int));
    }

    // store matrix
    for (int x = 0; x < matrix_x; x++)
    {
        for (int y = 0; y < matrix_y; y++)
        {
            scanf("%d", &matrix[x][y]);
        }
    }

    int max1 = find_max_area(matrix, 1);
    int max0 = find_max_area(matrix, 0);

    if (max1 > max0)
    {
        printf("%d\n", max1);
    }
    else
    {
        printf("%d\n", max0);
    }

    for (int i = 0; i < matrix_x; i++)
    {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}