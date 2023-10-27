#pragma once

using namespace std;

void testSpecialDotsCons(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsCons(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsConsVector(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsConsVector(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsOMP(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsOMP(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsOMPVector(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsOMPVector(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsTBB(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsTBB(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsTBBVector(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsTBBVector(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsOMPTBB(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsOMPTBB(rgb_in, rgb_out, height, width, threads);
}
void testSpecialDotsOMPTBBVector(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int& height, int& width, int& threads, int& ksize) {
    return specialDotsOMPTBBVector(rgb_in, rgb_out, height, width, threads);
}

char* inBMP(int i) {
    string str;
    str = "images\\input";
    char* cstr;
    switch (i)
    {
    case 1:
        str += "_1.bmp";
        break;
    case 2:
        str += "_2.bmp";
        break;
    case 3:
        str += "_3.bmp";
        break;
    default:
        break;
    }

    cstr = new char[str.length() + 1];
    strcpy(cstr, str.c_str());
    return cstr;
}
//void outBMP(char*& OUTPATHM2, char*& OUTPATHU, char*& OUTPATHR, char*& OUTPATHE, int i, int alg, int ksize) {
//    string* str = new string[4];
//    for (int k = 0; k < 4; k++) {
//        str[k] = "images\\output";
//    }
//
//    char** cstr;
//
//    for (int k = 0; k < 4; k++)
//    {
//        switch (k)
//        {
//        case 0:
//            str[k] += "_�������2�������";
//            break;
//        case 1:
//            str[k] += "_������������";
//            break;
//        case 2:
//            str[k] += "_����������������������";
//            break;
//        case 3:
//            str[k] += "_��������";
//            break;
//        default:
//            break;
//        }
//
//        switch (i)
//        {
//        case 1:
//            str[k] += "_1600X900";
//            break;
//        case 2:
//            str[k] += "_1920X1280";
//            break;
//        default:
//            break;
//        }
//
//        switch (alg)
//        {
//        case 0:
//            str[k] += "_����";
//            break;
//        case 1:
//            str[k] += "_���(OMPfor)";
//            break;
//        case 2:
//            str[k] += "_���(TBBfor)";
//            break;
//        default:
//            break;
//        }
//        switch (ksize)
//        {
//        case 5:
//            str[k] += "5X5";
//            break;
//        case 7:
//            str[k] += "7X7";
//            break;
//        case 9:
//            str[k] += "9X9";
//            break;
//        default:
//            break;
//        }
//
//        str[k] += ".bmp";
//    }
//
//    OUTPATHM2 = new char[str[0].length() + 1];
//    strcpy(OUTPATHM2, str[0].c_str());
//    OUTPATHU = new char[str[1].length() + 1];
//    strcpy(OUTPATHU, str[1].c_str());
//    OUTPATHR = new char[str[2].length() + 1];
//    strcpy(OUTPATHR, str[2].c_str());
//    OUTPATHE = new char[str[3].length() + 1];
//    strcpy(OUTPATHE, str[3].c_str());
//
//    delete[] str;
//}

double AvgTrustedInterval(double& avg, double*& times, int& cnt)
{
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
        if (avg - sd <= times[i] && times[i] <= avg + sd)
        {
            newAVg += times[i];
            newCnt++;
        }
    }
    if (newCnt == 0) newCnt = 1;
    return newAVg / newCnt;
}

vector<double> TestIter(void* Funct, RGBTRIPLE** rgb_input, int height, int width, int ksize, int iterations, int i, int alg, int threads)
{
    double* Times1 = new double[iterations];
    double* Times2 = new double[iterations];
    double* Times3 = new double[iterations];
    double* Times4 = new double[iterations];
    double* Times5 = new double[iterations];
    double* Times6 = new double[iterations];
    double* Times7 = new double[iterations];
    double* Times8 = new double[iterations];

    vector<double> avgTimeT;

    double avgTimeT1 = 0;
    double avgTimeT2 = 0;
    double avgTimeT3 = 0;
    double avgTimeT4 = 0;
    double avgTimeT5 = 0;
    double avgTimeT6 = 0;
    double avgTimeT7 = 0;
    double avgTimeT8 = 0;

    double t1_sum = 0;
    double t2_sum = 0;
    double t3_sum = 0;
    double t4_sum = 0;
    double t5_sum = 0;
    double t6_sum = 0;
    double t7_sum = 0;
    double t8_sum = 0;

    cout << endl;

    RGBTRIPLE** rgb_out = new RGBTRIPLE * [height];
    rgb_out[0] = new RGBTRIPLE[height * width];

    for (int j = 0; j < height; j++)
    {
        rgb_out[j] = &rgb_out[0][j * width];
    }

    for (int j = 0; j < iterations; j++)
    {
        ((*(TestFunctTemp1)Funct)(rgb_input, rgb_out, height, width, threads, ksize));
        Times1[j] = t1[0];
        Times2[j] = t2[0];
        Times3[j] = t3[0];
        Times4[j] = t4[0];
        Times5[j] = t5[0];
        Times6[j] = t6[0];
        Times7[j] = t7[0];
        Times8[j] = t8[0];
        t1_sum += t1[0];
        t2_sum += t2[0];
        t3_sum += t3[0];
        t4_sum += t4[0];
        t5_sum += t5[0];
        t6_sum += t6[0];
        t7_sum += t7[0];
        t8_sum += t8[0];
        t1.clear();
        t2.clear();
        t3.clear();
        t4.clear();
        t5.clear();
        t6.clear();
        t7.clear();
        t8.clear();
        cout << "+";
    }
    cout << endl;

    t1_sum /= iterations;
    cout << "AvgTime t1:" << t1_sum << endl;
    t2_sum /= iterations;
    cout << "AvgTime t2:" << t2_sum << endl;
    t3_sum /= iterations;
    cout << "AvgTime t3:" << t3_sum << endl;
    t4_sum /= iterations;
    cout << "AvgTime t4:" << t4_sum << endl;
    t5_sum /= iterations;
    cout << "AvgTime t5:" << t5_sum << endl;
    t6_sum /= iterations;
    cout << "AvgTime t6:" << t6_sum << endl;
    t7_sum /= iterations;
    cout << "AvgTime t7:" << t7_sum << endl;
    t8_sum /= iterations;
    cout << "AvgTime t8:" << t8_sum << endl;


    avgTimeT1 = AvgTrustedInterval(t1_sum, Times1, iterations);
    cout << "AvgTimeTrusted t1:" << avgTimeT1 << endl;
    avgTimeT2 = AvgTrustedInterval(t2_sum, Times2, iterations);
    cout << "AvgTimeTrusted t2:" << avgTimeT2 << endl;
    avgTimeT3 = AvgTrustedInterval(t3_sum, Times3, iterations);
    cout << "AvgTimeTrusted t3:" << avgTimeT3 << endl;
    avgTimeT4 = AvgTrustedInterval(t4_sum, Times4, iterations);
    cout << "AvgTimeTrusted t4:" << avgTimeT4 << endl;
    avgTimeT5 = AvgTrustedInterval(t5_sum, Times5, iterations);
    cout << "AvgTimeTrusted t5:" << avgTimeT5 << endl;
    avgTimeT6 = AvgTrustedInterval(t6_sum, Times6, iterations);
    cout << "AvgTimeTrusted t6:" << avgTimeT6 << endl;
    avgTimeT7 = AvgTrustedInterval(t7_sum, Times7, iterations);
    cout << "AvgTimeTrusted t7:" << avgTimeT7 << endl;
    avgTimeT8 = AvgTrustedInterval(t8_sum, Times8, iterations);
    cout << "AvgTimeTrusted t8:" << avgTimeT8 << endl;

    avgTimeT.push_back(t1_sum);
    avgTimeT.push_back(t2_sum);
    avgTimeT.push_back(t3_sum);
    avgTimeT.push_back(t4_sum);
    avgTimeT.push_back(t5_sum);
    avgTimeT.push_back(t6_sum);
    avgTimeT.push_back(t7_sum);
    avgTimeT.push_back(t8_sum);

    switch (i) {
    case 1: 
        BMPWrite(rgb_out, width, height, "output_1.bmp");
        break;
    case 2:
        BMPWrite(rgb_out, width, height, "output_2.bmp");
        break;
    case 3:
        BMPWrite(rgb_out, width, height, "output_3.bmp");
        break;
    default:
        break;
    }
    

    delete[] rgb_out[0];
    delete[] rgb_out;

    return avgTimeT;
}