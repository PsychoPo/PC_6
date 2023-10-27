#pragma once

void specialDotsTBB(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int imHeight, int imWidth, int threads, int ksize = 5)
{
	double** I = new double* [imHeight];
	double** IDX = new double* [imHeight];
	double** IDY = new double* [imHeight];
	double** IDXY = new double* [imHeight];
	double** IDX0 = new double* [imHeight];
	double** IDY0 = new double* [imHeight];
	double** R = new double* [imHeight];
	double** H = new double* [imHeight];

	for (int i = 0; i < imHeight; i++) {
		I[i] = new double[imWidth];
		IDX[i] = new double[imWidth];
		IDY[i] = new double[imWidth];
		IDXY[i] = new double[imWidth];
		IDX0[i] = new double[imWidth];
		IDY0[i] = new double[imWidth];
		R[i] = new double[imWidth];
		H[i] = new double[imWidth];
	}
	for (int y = 0; y < imHeight; y++) {
		for (int x = 0; x < imWidth; x++) {
			rgb_out[y][x].rgbtRed = rgb_in[y][x].rgbtRed;
			rgb_out[y][x].rgbtGreen = rgb_in[y][x].rgbtGreen;
			rgb_out[y][x].rgbtBlue = rgb_in[y][x].rgbtBlue;
		}
	}

	double start = omp_get_wtime();
	double start_time = omp_get_wtime();
	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					I[y][x] = rgb_in[y][x].rgbtRed * 0.299 + rgb_in[y][x].rgbtGreen * 0.587 + rgb_in[y][x].rgbtBlue * 0.114;
				}
			}
		});


	int size = ksize * ksize;
	int rh = ksize / 2;
	int rw = ksize / 2;
	int M = rh * 2 + 1;
	int N = rw * 2 + 1;
	double sigma = (double)ksize / 2;

	double** MATR_coef = new double* [M];
	for (int i = 0; i < M; i++)
		MATR_coef[i] = new double[N];

	double sum = tbb::parallel_reduce(
		tbb::blocked_range2d<int>(-rh, rh + 1, -rw, rw + 1), 0.0, [&](tbb::blocked_range2d<int> r, double total) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int yk = y + rh;
					int xk = x + rw;

					double CF = (1 / (2 * 3.14 * pow(sigma, 2))) * exp(-1 * ((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))));
					MATR_coef[yk][xk] = CF;
					total += MATR_coef[yk][xk];
				}
			}
	return total;
		}, std::plus<double>());

	tbb::parallel_for(
		tbb::blocked_range2d<int>(-rh, rh, -rw, rw), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int yk = y + rh;
					int xk = x + rw;

					MATR_coef[yk][xk] /= sum;
				}
			}
		});


	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_Value = 0;

					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;
						for (int dx = -rw; dx <= rw; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > imWidth - 1) kx = imWidth - 1;

							LinF_Value += double(I[ky][kx]) * MATR_coef[dy + rh][dx + rw];
						}
					}
					if (LinF_Value < 0) LinF_Value = 0;
					if (LinF_Value > 255) LinF_Value = 255;

					I[y][x] = LinF_Value;
				}
			}
		});


	t1.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	//rw = 1;
	//rh = 1;
	double DiffX[3][3] = { {1, 0, -1}, {1, 0, -1}, {1, 0, -1} };
	double DiffY[3][3] = { {1, 1, 1}, {0, 0, 0}, {-1, -1, -1} };

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_ValueDX = 0;
					double LinF_ValueDY = 0;

					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;
						for (int dx = -rw; dx <= rw; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > imWidth - 1) kx = imWidth - 1;

							LinF_ValueDX += double(I[ky][kx]) * DiffX[dy + rh][dx + rw];
							LinF_ValueDY += double(I[ky][kx]) * DiffY[dy + rh][dx + rw];
						}
					}

					if (LinF_ValueDX < 0) LinF_ValueDX = 0;
					if (LinF_ValueDX > 255) LinF_ValueDX = 255;
					if (LinF_ValueDY < 0) LinF_ValueDY = 0;
					if (LinF_ValueDY > 255) LinF_ValueDY = 255;

					IDX0[y][x] = LinF_ValueDX;
					IDY0[y][x] = LinF_ValueDY;
				}
			}
		});


	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_ValueDX = 0;
					double LinF_ValueDY = 0;
					double LinF_ValuedXDY = 0;
					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;

						for (int dx = -rw; dx <= rw; dx++) {
							int KX = x + dx;
							if (KX < 0) KX = 0;
							if (KX > imWidth - 1) KX = imWidth - 1;

							LinF_ValueDX += double(IDX0[ky][KX]) * DiffX[dy + rh][dx + rw];
							LinF_ValueDY += double(IDY0[ky][KX]) * DiffY[dy + rh][dx + rw];
							LinF_ValuedXDY += double(IDX0[ky][KX]) * DiffY[dy + rh][dx + rw];
						}
					}
					if (LinF_ValueDX < 0) LinF_ValueDX = 0;
					if (LinF_ValueDX > 255) LinF_ValueDX = 255;
					if (LinF_ValueDY < 0) LinF_ValueDY = 0;
					if (LinF_ValueDY > 255) LinF_ValueDY = 255;

					IDX[y][x] = LinF_ValueDX;
					IDY[y][x] = LinF_ValueDY;
					IDXY[y][x] = LinF_ValuedXDY;
				}
			}
		});


	t2.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	double k = 0.04;
	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double A = IDX[y][x];
					double B = IDY[y][x];
					double C = IDXY[y][x];
					R[y][x] = (A * B - C * C) - (k * ((A + B) * (A + B)));
				}
			}
		});


	vector<Dot> harrisPoints;
	Dot tmp;
	double min_v = DBL_MAX;
	double max_v = -DBL_MAX;

	max_v = tbb::parallel_reduce(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), 0.0, [&](tbb::blocked_range2d<int> r, double max) {
			for (int i = r.rows().begin(); i < r.rows().end(); i++) {
				for (int j = r.cols().begin(); j < r.cols().end(); j++) {

					if (R[i][j] > max_v) max_v = R[i][j];

				}
			}
			return max;
		}, std::plus<double>());
	min_v = tbb::parallel_reduce(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), 0.0, [&](tbb::blocked_range2d<int> r, double min) {
			for (int i = r.rows().begin(); i < r.rows().end(); i++) {
				for (int j = r.cols().begin(); j < r.cols().end(); j++) {

					if (R[i][j] > min) min = R[i][j];

				}
			}
			return min;
		}, std::plus<double>());

	double TS = (max_v - min_v) * k + min_v;

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					if (R[y][x] > TS) {
						H[y][x] = R[y][x];
						tmp.R = R[y][x];
						tmp.X = x;
						tmp.Y = y;
					}
					else {
						H[y][x] = 0;
					}
				}
			}
		});

	t3.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	for (int y = 0; y < imHeight; y++) {
		for (int x = 0; x < imWidth; x++) {
			if (H[y][x] > 0) {
				tmp.R = R[y][x];
				tmp.X = x;
				tmp.Y = y;
				harrisPoints.push_back(tmp);
			}
		}
	}

	t4.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	harrisPoints = shell(harrisPoints.size(), harrisPoints);

	t5.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	int i = 0;
	bool flag = 1;

	do {
		flag = 1;
		if (i > 0) {
			for (int k = 0; k < i; k++) {
				if (sqrt(pow(harrisPoints[i].X - harrisPoints[k].X, 2) + pow(harrisPoints[i].Y - harrisPoints[k].Y, 2)) < 20) {
					harrisPoints.erase(harrisPoints.begin() + i);
					flag = 0;
					break;
				}
			}
		}
		if (flag != 0)
			i++;
	} while (i < harrisPoints.size());

	t6.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	tbb::parallel_for(
		tbb::blocked_range<int>(0, harrisPoints.size()), [&](tbb::blocked_range<int> r) {
			for (int i = r.begin(); i < r.end(); i++) {
				for (int r = -2; r < 2; r++) {
					for (int c = -2; c < 2; c++) {

						int yk = harrisPoints[i].Y + r;
						int xk = harrisPoints[i].X + c;
						//if (yk > 0 && xk > 0) {
						if (yk >= 0 && xk >= 0 && yk < imHeight && xk < imWidth) {
							rgb_out[yk][xk].rgbtRed = 255; //���
							rgb_out[yk][xk].rgbtGreen = 0;
							rgb_out[yk][xk].rgbtBlue = 0;
						}
					}
				}
			}
		});
	t7.push_back((omp_get_wtime() - start_time) * 1000);
	t8.push_back((omp_get_wtime() - start) * 1000);

	for (int i = 0; i < imHeight; i++) {
		delete[] I[i];
		delete[] IDX[i];
		delete[] IDY[i];
		delete[] IDXY[i];
		delete[] IDX0[i];
		delete[] IDY0[i];
		delete[] R[i];
		delete[] H[i];
	}
	delete[] I;
	delete[] IDX;
	delete[] IDY;
	delete[] IDXY;
	delete[] IDX0;
	delete[] IDY0;
	delete[] R;
	delete[] H;
}

void specialDotsTBBVector(RGBTRIPLE**& rgb_in, RGBTRIPLE**& rgb_out, int imHeight, int imWidth, int threads, int ksize = 5)
{
	concurrent_vector <concurrent_vector <double>> I(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> IDX(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> IDY(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> IDXY(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> IDX0(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> IDY0(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> R(imHeight, concurrent_vector <double>(imWidth));
	concurrent_vector <concurrent_vector <double>> H(imHeight, concurrent_vector <double>(imWidth));

	for (int y = 0; y < imHeight; y++) {
		for (int x = 0; x < imWidth; x++) {
			rgb_out[y][x].rgbtRed = rgb_in[y][x].rgbtRed;
			rgb_out[y][x].rgbtGreen = rgb_in[y][x].rgbtGreen;
			rgb_out[y][x].rgbtBlue = rgb_in[y][x].rgbtBlue;
		}
	}

	double start = omp_get_wtime();
	double start_time = omp_get_wtime();
	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					I[y][x] = rgb_in[y][x].rgbtRed * 0.299 + rgb_in[y][x].rgbtGreen * 0.587 + rgb_in[y][x].rgbtBlue * 0.114;
				}
			}
		});


	int size = ksize * ksize;
	int rh = ksize / 2;
	int rw = ksize / 2;
	int M = rh * 2 + 1;
	int N = rw * 2 + 1;
	double sigma = (double)ksize / 2;

	concurrent_vector <concurrent_vector <double>> MATR_coef(M, concurrent_vector <double>(N));

	double sum = tbb::parallel_reduce(
		tbb::blocked_range2d<int>(-rh, rh + 1, -rw, rw + 1), 0.0, [&](tbb::blocked_range2d<int> r, double total) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int yk = y + rh;
					int xk = x + rw;

					double CF = (1 / (2 * 3.14 * pow(sigma, 2))) * exp(-1 * ((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))));
					MATR_coef[yk][xk] = CF;
					total += MATR_coef[yk][xk];
				}
			}
	return total;
		}, std::plus<double>());

	tbb::parallel_for(
		tbb::blocked_range2d<int>(-rh, rh, -rw, rw), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int yk = y + rh;
					int xk = x + rw;

					MATR_coef[yk][xk] /= sum;
				}
			}
		});


	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_Value = 0;

					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;
						for (int dx = -rw; dx <= rw; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > imWidth - 1) kx = imWidth - 1;

							LinF_Value += double(I[ky][kx]) * MATR_coef[dy + rh][dx + rw];
						}
					}
					if (LinF_Value < 0) LinF_Value = 0;
					if (LinF_Value > 255) LinF_Value = 255;

					I[y][x] = LinF_Value;
				}
			}
		});

	t1.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	//rw = 1;
	//rh = 1;
	double DiffX[3][3] = { {1, 0, -1}, {1, 0, -1}, {1, 0, -1} };
	double DiffY[3][3] = { {1, 1, 1}, {0, 0, 0}, {-1, -1, -1} };

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_ValueDX = 0;
					double LinF_ValueDY = 0;

					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;
						for (int dx = -rw; dx <= rw; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > imWidth - 1) kx = imWidth - 1;

							LinF_ValueDX += double(I[ky][kx]) * DiffX[dy + rh][dx + rw];
							LinF_ValueDY += double(I[ky][kx]) * DiffY[dy + rh][dx + rw];
						}
					}

					if (LinF_ValueDX < 0) LinF_ValueDX = 0;
					if (LinF_ValueDX > 255) LinF_ValueDX = 255;
					if (LinF_ValueDY < 0) LinF_ValueDY = 0;
					if (LinF_ValueDY > 255) LinF_ValueDY = 255;

					IDX0[y][x] = LinF_ValueDX;
					IDY0[y][x] = LinF_ValueDY;
				}
			}
		});


	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double LinF_ValueDX = 0;
					double LinF_ValueDY = 0;
					double LinF_ValuedXDY = 0;
					for (int dy = -rh; dy <= rh; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > imHeight - 1) ky = imHeight - 1;

						for (int dx = -rw; dx <= rw; dx++) {
							int KX = x + dx;
							if (KX < 0) KX = 0;
							if (KX > imWidth - 1) KX = imWidth - 1;

							LinF_ValueDX += double(IDX0[ky][KX]) * DiffX[dy + rh][dx + rw];
							LinF_ValueDY += double(IDY0[ky][KX]) * DiffY[dy + rh][dx + rw];
							LinF_ValuedXDY += double(IDX0[ky][KX]) * DiffY[dy + rh][dx + rw];
						}
					}
					if (LinF_ValueDX < 0) LinF_ValueDX = 0;
					if (LinF_ValueDX > 255) LinF_ValueDX = 255;
					if (LinF_ValueDY < 0) LinF_ValueDY = 0;
					if (LinF_ValueDY > 255) LinF_ValueDY = 255;

					IDX[y][x] = LinF_ValueDX;
					IDY[y][x] = LinF_ValueDY;
					IDXY[y][x] = LinF_ValuedXDY;
				}
			}
		});


	t2.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	double k = 0.04;
	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					double A = IDX[y][x];
					double B = IDY[y][x];
					double C = IDXY[y][x];
					R[y][x] = (A * B - C * C) - (k * ((A + B) * (A + B)));
				}
			}
		});


	vector<Dot> harrisPoints;
	Dot tmp;
	//double max = R[0][0];
	double min = DBL_MAX;
	double max = -DBL_MAX;
	double TS = (max - min) * k + min;
	max = tbb::parallel_reduce(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), 0.0, [&](tbb::blocked_range2d<int> r, double max) {
			for (int i = r.rows().begin(); i < r.rows().end(); i++) {
				for (int j = r.cols().begin(); j < r.cols().end(); j++) {
					
					if (R[i][j] > max) max = R[i][j];
					
				}
			}
			return max;
		}, std::plus<double>());

	//double max = -DBL_MAX;
	//double min = DBL_MAX;
	/*tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int i = r.rows().begin(); i < r.rows().end(); i++) {
				for (int j = r.cols().begin(); j < r.cols().end(); j++) {
					if (R[i][j] > max) max = R[i][j];
				}
			}
		});*/

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, imHeight, 0, imWidth), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					if (R[y][x] > TS) {
						H[y][x] = R[y][x];
						tmp.R = R[y][x];
						tmp.X = x;
						tmp.Y = y;
					}
					else {
						H[y][x] = 0;
					}
				}
			}
		});

	t3.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	for (int y = 0; y < imHeight; y++) {
		for (int x = 0; x < imWidth; x++) {
			if (H[y][x] > 0) {
				tmp.R = R[y][x];
				tmp.X = x;
				tmp.Y = y;
				harrisPoints.push_back(tmp);
			}
		}
	}
	t4.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	harrisPoints = shell(harrisPoints.size(), harrisPoints);

	t5.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	int i = 0;
	bool flag = 1;

	do {
		flag = 1;
		if (i > 0) {
			for (int k = 0; k < i; k++) {
				if (sqrt(pow(harrisPoints[i].X - harrisPoints[k].X, 2) + pow(harrisPoints[i].Y - harrisPoints[k].Y, 2)) < 20) {
					harrisPoints.erase(harrisPoints.begin() + i);
					flag = 0;
					break;
				}
			}
		}
		if (flag != 0)
			i++;
	} while (i < harrisPoints.size());

	t6.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	tbb::parallel_for(
		tbb::blocked_range<int>(0, harrisPoints.size()), [&](tbb::blocked_range<int> r) {
			for (int i = r.begin(); i < r.end(); i++) {
				for (int r = -2; r < 2; r++) {
					for (int c = -2; c < 2; c++) {

						int yk = harrisPoints[i].Y + r;
						int xk = harrisPoints[i].X + c;
						if (yk >= 0 && xk >= 0 && yk < imHeight && xk < imWidth) {
							rgb_out[yk][xk].rgbtRed = 255;
							rgb_out[yk][xk].rgbtGreen = 0;
							rgb_out[yk][xk].rgbtBlue = 0;
						}
					}
				}
			}
		});
	t7.push_back((omp_get_wtime() - start_time) * 1000);
	t8.push_back((omp_get_wtime() - start) * 1000);

	I.clear();
	IDX.clear();
	IDY.clear();
	IDXY.clear();
	IDX0.clear();
	IDY0.clear();
	R.clear();
	H.clear();
}