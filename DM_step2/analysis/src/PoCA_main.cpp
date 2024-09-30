#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

#include "../header/PoCA.h"

using namespace std;
extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;

int main(int argc, char **argv)
{
    if (argc > 2 ||
        strcmp(argv[1], "--help") == 0)
    {
        cout << "Usage:\n    "
             << argv[0] << " (max_point_count)" << endl;
        return 0;
    }

    int max_point_count = INT32_MAX;
    if (argc == 2)
    {
        try
        {
            max_point_count = std::stoi(argv[1]);
        }
        catch (exception &err)
        {
            cerr << err.what() << endl;
            return 1;
        }
    }

    vector<int> count(SPACE_N);
    vector<double> sig(SPACE_N);
    vector<double> lambda(SPACE_N);

    struct timeval start;
    struct timeval end;
    unsigned long timer;

    gettimeofday(&start, NULL); // 计时开始

    int point_count = 0;
    int bad_point_count = 0;
    double E, x[4], y[4];
    while (cin >> x[0] >> y[0] >>
               x[1] >> y[1] >>
               x[2] >> y[2] >>
               x[3] >> y[3] >> E &&
           point_count < max_point_count)
    {
        ++point_count;
        if (!CollectPath({x[0], y[0], Z1}, // p1
                         {x[1], y[1], Z2}, // p2
                         {x[2], y[2], Z3}, // p3
                         {x[3], y[3], Z4}, // p4
                         E, count, sig))
        {
            ++bad_point_count;
        }
    }

    gettimeofday(&end, NULL); // 计时结束
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("time = %ld us\n", timer);

    if (point_count)
    {
        cerr << "Summary: good/total ="
             << point_count - bad_point_count << '/' << point_count << '\n';

        for (int i = 0; i < SPACE_N; ++i)
            lambda[i] = count[i] == 0 ? 0 : sig[i] / count[i];

        DumpStep(lambda, "PoCA Algorithm", "poca_cache", 0, 0); // 注意要加第五个参数 flag，否则无法生成二维图像

        vector<double> med_sig(SPACE_N);
        median_filter(sig, med_sig);
        DumpStep(med_sig, "PoCA Algorithm", "poca_median_cache", 0, 0); // 注意要加第五个参数 flag，否则无法生成二维图像
    }
    return 0;
}
