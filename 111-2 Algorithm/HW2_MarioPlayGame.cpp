#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    int input;
    int count = 0;
    int arr[100000] = {};

    // 將input值存入array中
    while (scanf("%d", &input) != EOF)
    {
        arr[count] = input;
        count += 1;
    }

    // 先跑第一遍array中的元素並找出最大的總和
    int tmp_max = 0;
    int final_max = 0;

    for (int i = 0; i < count; i++)
    {
        int add = tmp_max + arr[i];

        if (add <= 0) // 回傳負值或是歸零都等於0
        {
            tmp_max = 0;
        }
        else
        {
            tmp_max = add;

            if (tmp_max > final_max)
            {
                final_max = tmp_max;
            }
        }
    }

    // 再來將array裡的元素加總起來減掉array的最小總和
    int arr_sum = 0;
    int min_sum = 0;
    int small = 0;
    int final_max_2 = 0;

    for (int j = 0; j < count; j++)
    {
        int num = arr[j];
        arr_sum += num;

        // 找出array裡最小的總和
        min_sum += num;
        if (min_sum > num)
        {
            min_sum = num;
        }

        // 回傳最終最小值
        if (small > min_sum)
        {
            small = min_sum;
        }
    }

    final_max_2 = arr_sum - small;

    if (final_max > final_max_2)
    {
        printf("%d", final_max);
    }
    else
    {
        printf("%d", final_max_2);
    }

    printf("\n");

    return 0;
}