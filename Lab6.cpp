#include <iostream>
#include <omp.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/cache_aligned_allocator.h>
#include "BMPFileRW.h"
#include <iomanip>
using namespace std;
using namespace tbb;

struct Dot
{
	double R;
	int X, Y;
};

vector<double> t1, t2, t3, t4, t5, t6, t7, t8;

const char* images[] = {
	"C:\\Users\\popov\\Desktop\\1.bmp",
	"C:\\Users\\popov\\Desktop\\1280x1024.bmp",
	"C:\\Users\\popov\\Desktop\\1920x1080.bmp",
	"C:\\Users\\popov\\Desktop\\1920x2877.bmp"
};

#include "sorts.h"
#include "special_dots.h"
#include "special_dots_omp.h"
#include "special_dots_tbb.h"
#include "omp_tbb.h"

void((*functions[]))(RGBTRIPLE**&, RGBTRIPLE**&, int, int, int) = {
    special_dots,
    special_dots_concurrent_vector,

    special_dots_omp,
    special_dots_omp_concurrent_vector,

    special_dots_tbb,
    special_dots_tbb_concurrent_vector,

    special_dots_omp_tbb,
    special_dots_omp_tbb_concurrent_vector
};

const char* function_names[] = {
    "Спец. точки последовательно",
    "Спец. точки последовательно vector",

    "Спец. точки OMP",
    "Спец. точки OMP vector",

    "Спец. точки TBB",
    "Спец. точки TBB vector",

    "Спец. точки OMP-TBB",
    "Спец. точки OMP-TBB vector"
};

void bubbleSort(vector<double> mas, int cnt) {
    bool exchanges;
    do {
        exchanges = false;
        for (int i = 0; i < cnt - 1; i++) {
            if (mas[i] > mas[i + 1]) {
                double temp = mas[i];
                mas[i] = mas[i + 1];
                mas[i + 1] = temp;
                exchanges = true;
            }
        }
    } while (exchanges);
}

double AvgTrustedIntervalMed(vector<double> times, int cnt)
{
    double avg = 0;
    for (int i = 0; i < cnt; i++)
    {
        avg += times[i];
    }
    avg /= cnt;
    bubbleSort(times, cnt);
    double med = times[cnt / 2];
    double sd = 0, newAVg = 0;
    int newCnt = 0;
    for (int i = 0; i < cnt; i++)
    {
        sd += (times[i] - avg) * (times[i] - avg);
    }
    sd /= (cnt - 1.0);
    sd = sqrt(sd);
    for (int i = 0; i < cnt; i++)
    {
        if (med - sd <= times[i] && times[i] <= med + sd)
        {
            newAVg += times[i];
            newCnt++;
        }
    }
    if (newCnt == 0) newCnt = 1;
    return newAVg / newCnt;
}

int main() {
    setlocale(LC_ALL, "en_US.UTF-8");
	double times_sum1 = 0;
	double times_sum2 = 0;
	double times_sum3 = 0;
	double times_sum4 = 0;
	double times_sum5 = 0;
	double times_sum6 = 0;
	double times_sum7 = 0;
	double times_sum8 = 0;

	//ofstream full_time("full_time.txt");
	//ofstream time_functions("time_functions.txt");
	//ofstream omp_tbb("omp_tbb_test.txt");

    cout.precision(3);
    //time_functions.precision(3);
    //omp_tbb.precision(3);

	for (int func = 0; func < 8; func++) {
        cout << function_names[func] << ":\n";
        //full_time << function_names[func] << ":\n";
		for (int th = 2; th <= 4; th++) {
			cout << "П-" << th << ": ";
            //full_time << "П-" << th << ": ";
            omp_set_num_threads(th);
            tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, th);
			for (int image = 0; image < 4; image++) {
				RGBTRIPLE** rgb_input, ** rgb_output;
				BITMAPFILEHEADER header;
				BITMAPINFOHEADER bmiHeader;
                BMPRead(rgb_input, header, bmiHeader, images[image]);

                int imWidth = 0, imHeight = 0;
                imWidth = bmiHeader.biWidth;
                imHeight = bmiHeader.biHeight;

                rgb_output = new RGBTRIPLE * [imHeight];
                rgb_output[0] = new RGBTRIPLE[imWidth * imHeight];
                for (int i = 1; i < imHeight; i++)
                {
                    rgb_output[i] = &rgb_output[0][imWidth * i];
                }

                for (int i = 0; i < 10; i++)
                    functions[func](rgb_input, rgb_output, imHeight, imWidth, 5);
                times_sum8 = AvgTrustedIntervalMed(t8, 10);
                cout.precision(3);
                //full_time.precision(3);
                cout << setw(7) << times_sum8 << "  " << fixed;
                //full_time << times_sum8 << "  " << fixed;

				BMPWrite(rgb_output, imWidth, imHeight, "out.bmp");

                t1.clear();
                t2.clear();
                t3.clear();
                t4.clear();
                t5.clear();
                t6.clear();
                t7.clear();
                t8.clear();
			}
			cout << endl;
            //full_time << endl;
		}
        cout << endl;
        //full_time << endl;
	}

    for (int func = 0; func < 8; func++) {
        cout << function_names[func] << ":\n";
        //time_functions << function_names[func] << ":\n";
        for (int th = 2; th <= 4; th++) {
            cout << "П-" << th << ": ";
            //time_functions << "П-" << th << ": ";
            omp_set_num_threads(th);
            tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, th);
            for (int image = 0; image < 4; image++) {
                RGBTRIPLE** rgb_input, ** rgb_output;
                BITMAPFILEHEADER header;
                BITMAPINFOHEADER bmiHeader;
                BMPRead(rgb_input, header, bmiHeader, images[image]);

                int imWidth = 0, imHeight = 0;
                imWidth = bmiHeader.biWidth;
                imHeight = bmiHeader.biHeight;

                rgb_output = new RGBTRIPLE * [imHeight];
                rgb_output[0] = new RGBTRIPLE[imWidth * imHeight];
                for (int i = 1; i < imHeight; i++)
                {
                    rgb_output[i] = &rgb_output[0][imWidth * i];
                }

                for (int i = 0; i < 10; i++)
                    functions[func](rgb_input, rgb_output, imHeight, imWidth, 5);
                times_sum1 = AvgTrustedIntervalMed(t1, 10);
                times_sum2 = AvgTrustedIntervalMed(t2, 10);
                times_sum3 = AvgTrustedIntervalMed(t3, 10);
                times_sum4 = AvgTrustedIntervalMed(t4, 10);
                times_sum5 = AvgTrustedIntervalMed(t5, 10);
                times_sum6 = AvgTrustedIntervalMed(t6, 10);
                times_sum7 = AvgTrustedIntervalMed(t7, 10);

                cout << times_sum1 << "  " << fixed;
                cout << times_sum2 << "  " << fixed;
                cout << times_sum3 << "  " << fixed;
                cout << times_sum4 << "  " << fixed;
                cout << times_sum5 << "  " << fixed;
                cout << times_sum6 << "  " << fixed;
                cout << times_sum7 << fixed;
                cout << " | ";

                /*time_functions << times_sum1 << "  " << fixed;
                time_functions << times_sum2 << "  " << fixed;
                time_functions << times_sum3 << "  " << fixed;
                time_functions << times_sum4 << "  " << fixed;
                time_functions << times_sum5 << "  " << fixed;
                time_functions << times_sum6 << "  " << fixed;
                time_functions << times_sum7 << fixed;
                time_functions << " | ";*/

                BMPWrite(rgb_output, imWidth, imHeight, "out.bmp");

                t1.clear();
                t2.clear();
                t3.clear();
                t4.clear();
                t5.clear();
                t6.clear();
                t7.clear();
                t8.clear();
            }
            cout << endl;
            //time_functions << endl;
        }
        cout << endl;
        //time_functions << endl;
    }
}