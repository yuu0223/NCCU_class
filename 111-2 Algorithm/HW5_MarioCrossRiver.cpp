#include <stdio.h>
#include <stdlib.h>

struct node_info
{
    int start;  // 木板起始點
    int end;    // 木板結束點
    int travel; // 是否路過這個木板了 0為是 1為否
};

int next_step(int position, struct node_info wood[], int wood_num)
{
    int temp_max_end = 0;
    int max_position = 0;
    int final_max = 0;
    for (int i = 0; i < wood_num; i++)
    {
        if (wood[i].travel == 0)
        {
            if (wood[i].start < wood[position].end)
            {
                temp_max_end = wood[i].end;
                wood[i].travel = 1;
                if (temp_max_end > final_max)
                {
                    if (temp_max_end >= wood[wood_num - 1].end)
                    {
                        final_max = wood[wood_num - 1].end;
                        max_position = wood_num - 1;
                    }
                    else
                    {
                        final_max = temp_max_end;
                        max_position = i;
                    }
                }
            }
        }
    }

    return max_position;
}

int main(void)
{
    struct node_info wood[1001];

    int input;
    int i = 0;
    int num = 0;
    int wood_num;
    while (scanf("%d", &input) != EOF)
    {
        if (i == 0)
        {
            wood_num = input;
        }
        else if (i % 2 == 1)
        {
            wood[num].start = input;
        }
        else
        {
            wood[num].end = input;
            wood[num].travel = 0;
            num++;
        }
        i++;
    }

    // 開始跳木板
    int position = 0;
    int step = 0;
    int next_list[1001] = {}; // 選擇最短路徑 會到達的點

    wood[0].travel = 1;
    while (position != (wood_num - 1))
    {
        position = next_step(position, wood, wood_num); // 回傳下一個值在wood中的位置
        next_list[step] = position;
        step++;
    }

    printf("%d\n", step);

    // int n = 0;
    // while (n != wood_num)
    // {
    //     printf("n:%d %d\n", n, wood[n].travel);
    //     n++;
    // }

    return 0;
}

// 用來檢測有沒有input成功的
// for (int i = 0; i < wood_num; i++)
// {
//     printf("start:%d", wood[i].start);
//     printf(" ,end:%d", wood[i].end);
//     printf(" ,dist:%d", wood[i].distance);
//     printf(" ,travel:%d\n", wood[i].travel);
// }