#include <stdio.h>
#include <stdlib.h>

struct work_time
{
    int start;
    int end;
    int days;
};

void sort(struct work_time arr[], int size)
{
    for (int i = 1; i < size; i++)
    {
        struct work_time temp = arr[i];
        int s = i - 1;
        while (s >= 0 && arr[s].end > temp.end)
        {
            arr[s + 1] = arr[s];
            s--;
        }
        arr[s + 1] = temp;
    }
}

int main(void)
{
    struct work_time real_time[1001];

    int input;
    int client; // 顧客數量
    int i = 0;
    int num = 0;
    while (scanf("%d", &input) != EOF)
    {
        if (i == 0)
        {
            client = input;
        }
        else if (i % 2 == 1)
        {
            real_time[num].start = input;
        }
        else
        {
            real_time[num].end = input;
            real_time[num].days = real_time[num].end - real_time[num].start;
            num++;
        }
        i++;
    }

    int start_work = real_time[client].start;
    int end_work = real_time[client].end;

    // 先將工作進行排序(order by工作時間end)
    sort(real_time, client);

    // 將超出工作時間的工作刪除
    struct work_time new_sort_time[client];
    int count = 0;
    for (int i = 0; i < client; i++)
    {
        if (real_time[i].start >= start_work && real_time[i].end <= end_work)
        {
            new_sort_time[count] = real_time[i];
            count++;
        }
    }

    // 進行工作排程
    int final_max = 0;
    int day[count];
    for (int a = 0; a < count; a++)
    {
        day[a] = new_sort_time[a].days;
        for (int b = 0; b < a; b++)
        {
            if (new_sort_time[b].end <= new_sort_time[a].start)
            {
                day[a] = new_sort_time[a].days + day[b];
            }

            if (day[a] > final_max)
            {
                final_max = day[a];
            }
        }
    }

    printf("%d\n", final_max);
    return 0;
}
