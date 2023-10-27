#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <cstring>
#include <math.h>
#include <string>
#include <tbb/tbb.h>
#include "BMPFileRW.h"

using namespace tbb;
using namespace std;

vector<double> t1, t2, t3, t4, t5, t6, t7, t8;

typedef double(*TestFunctTemp1)(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize);

struct Dot
{
	double R;
	int X, Y;
};

#include "Sort.h"
#include "SpecialDotsCons.h"
#include "SpecialDotsOMP.h"
#include "SpecialDotsTBB.h"
#include "SpecialDotsOMPTBB.h"
#include "Tests.h"



void test_functions(void** Functions, string(&function_names)[8])
{
    RGBTRIPLE** rgb_in;
    BITMAPFILEHEADER header;
    BITMAPINFOHEADER bmiHeader;
    int imWidth = 0, imHeight = 0;

    int iters = 2;
    int nd = 0;
    double times[3][8][8][3];
    for (int i = 1; i < 4; i++)
    {
        char* cstr = inBMP(i);
        BMPRead(rgb_in, header, bmiHeader, cstr);
        imWidth = bmiHeader.biWidth;
        imHeight = bmiHeader.biHeight;

        for (int threads = 1; threads <= 4; threads++)
        {
            omp_set_num_threads(threads);
            global_control global_limit(global_control::max_allowed_parallelism, threads);

            for (int alg = 0; alg < 6; alg++)
            {
                    if (threads == 1)
                    {
                        if (alg == 0 || alg == 1) {
                            vector<double> tn = TestIter(Functions[alg], rgb_in, imHeight, imWidth, 5, iters, i, alg, threads);

                            for (int t = 0; t < 8; t++)
                            {
                                times[nd][alg][t][0] = tn[t];
                                times[nd][alg][t][1] = times[nd][alg][t][0];
                                times[nd][alg][t][2] = times[nd][alg][t][0];
                            }
                        }
                    }
                    else
                    {
                        if (alg != 0 && alg != 1)
                        {
                            vector<double> tn = TestIter(Functions[alg], rgb_in, imHeight, imWidth, 5, iters, i, alg, threads);

                            for (int t = 0; t < 8; t++)
                            {
                                times[nd][alg][t][threads - 2] = tn[t];
                            }
                        }
                    }
            }
        }
        delete[] cstr;
        nd++;
    }
    ofstream fout1("output1.txt");
    fout1.imbue(locale("Russian"));
    ofstream fout2("output2.txt");
    fout2.imbue(locale("Russian"));
    for (int ND = 0; ND < 3; ND++)
    {
        switch (ND)
        {
        case 0:
            cout << "\n----------640*480----------" << endl;
            break;
        case 1:
            cout << "\n----------1024*768----------" << endl;
            break;
        case 2:
            cout << "\n----------1920*1285----------" << endl;
            break;
        default:
            break;
        }
        for (int t = 0; t < 7; t++)
        {
        for (int alg = 0; alg < 6; alg++)
        {
            for (int threads = 1; threads <= 4; threads++)
            {
                cout << "����� " << threads << " --------------" << endl;
                
                    if (threads == 1)
                    {
                        if (alg == 0 || alg == 1) {
                            cout << function_names[alg] << "\t" << times[ND][alg][t][0] << " ms." << endl;
                            if (t != 7) {
                                fout1 << times[ND][alg][t][0] << endl;
                            }
                            else if (t == 7) {
                                fout2 << times[ND][alg][t][0] << endl;
                            }
                        }
                    }
                    else
                    {
                        if (alg != 0 && alg != 1)
                        {
                            cout << function_names[alg] << "\t" << times[ND][alg][t][threads - 2] << " ms." << endl;
                            if (t != 7) {
                                fout1 << times[ND][alg][t][threads - 2] << endl;
                            }
                            else if (t == 7) {
                               fout2 << times[ND][alg][t][threads - 2] << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    fout1.close();
    fout2.close();
}

int main() {
    setlocale(LC_ALL, "RUS");

    void** Functions = new void* [8] { testSpecialDotsCons, testSpecialDotsConsVector, testSpecialDotsOMP, testSpecialDotsOMPVector, testSpecialDotsTBB, testSpecialDotsTBBVector,
        testSpecialDotsOMPTBB, testSpecialDotsOMPTBBVector};

    string function_names[8]{ "����. ����� ���������������", "����. ����� ���������������(�������)", "����. ����� OMP", "����. ����� OMP(�������)",
                              "����. ����� TBB", "����. ����� TBB(�������)", "����. ����� OMP + TBB", "����. ����� OMP + TBB(�������)" };
    test_functions(Functions, function_names);

	/*RGBTRIPLE** rgb_in, ** rgb_out;
	BITMAPFILEHEADER header;
	BITMAPINFOHEADER bmiHeader;
	int width = 0, height = 0;
	BMPRead(rgb_in, header, bmiHeader, "images\\input_3.bmp");
	width = bmiHeader.biWidth;
	height = bmiHeader.biHeight;

	rgb_out = new RGBTRIPLE * [height];
	rgb_out[0] = new RGBTRIPLE[height * width];

	for (int i = 0; i < height; i++)
	{
		rgb_out[i] = &rgb_out[0][i * width];
	}


	std::cout << "Image params:" << width << "x" << height << std::endl;

	specialDotsCons(rgb_in, rgb_out, height, width, 4);

	BMPWrite(rgb_out, width, height, "output_3.bmp");
	std::cout << "Image saved\n";*/

	return 0;
}