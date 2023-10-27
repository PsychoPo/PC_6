void special_dots_omp(RGBTRIPLE**& rgb, RGBTRIPLE**& rgb_out, int height, int width, int ksize = 5)
{
	double** I = new double* [height];
	double** IDX = new double* [height];
	double** IDY = new double* [height];
	double** IDXY = new double* [height];
	double** IDX0 = new double* [height];
	double** IDY0 = new double* [height];
	double** R = new double* [height];
	double** H = new double* [height];

	for (int i = 0; i < height; i++) {
		I[i] = new double[width];
		IDX[i] = new double[width];
		IDY[i] = new double[width];
		IDXY[i] = new double[width];
		IDX0[i] = new double[width];
		IDY0[i] = new double[width];
		R[i] = new double[width];
		H[i] = new double[width];
	}

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++) {
			rgb_out[y][x].rgbtRed = rgb[y][x].rgbtRed;
			rgb_out[y][x].rgbtGreen = rgb[y][x].rgbtGreen;
			rgb_out[y][x].rgbtBlue = rgb[y][x].rgbtBlue;
		}

	double start = omp_get_wtime();
	double start_time = omp_get_wtime();
#pragma omp parallel for
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			I[y][x] = rgb[y][x].rgbtRed * 0.299 + rgb[y][x].rgbtGreen * 0.587 + rgb[y][x].rgbtBlue * 0.114;

	int size = ksize * ksize;
	int rh = ksize / 2;
	int rw = ksize / 2;
	int M = rh * 2 + 1;
	int N = rw * 2 + 1;
	double sigma = (double)ksize / 2;

	double** MATR_coef = new double* [M];
	for (int i = 0; i < M; i++)
		MATR_coef[i] = new double[N];

	double sum = 0;

#pragma omp parallel for
	for (int y = -rh; y <= rh; y++) {
		for (int x = -rw; x <= rw; x++) {
			int yk = y + rh;
			int xk = x + rw;

			double CF = (1 / (2 * 3.14 * pow(sigma, 2))) * exp(-1 * ((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))));
			MATR_coef[yk][xk] = CF;
			sum += MATR_coef[yk][xk];
		}
	}

#pragma omp parallel for
	for (int y = -rh; y <= rh; y++) {
		for (int x = -rw; x <= rw; x++) {
			int yk = y + rh;
			int xk = x + rw;

			MATR_coef[yk][xk] /= sum;
		}
	}

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double LinF_Value = 0;

			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

					LinF_Value += double(I[ky][kx]) * MATR_coef[dy + rh][dx + rw];
				}
			}
			if (LinF_Value < 0) LinF_Value = 0;
			if (LinF_Value > 255) LinF_Value = 255;

			I[y][x] = LinF_Value;
		}
	}

	t1.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	//rw = 1;
	//rh = 1;
	double DiffX[3][3] = { {1, 0, -1}, {1, 0, -1}, {1, 0, -1} };
	double DiffY[3][3] = { {1, 1, 1}, {0, 0, 0}, {-1, -1, -1} };
	double LinF_ValueDX = 0;
	double LinF_ValueDY = 0;
	double LinF_ValuedXDY = 0;

#pragma omp parallel for firstprivate(LinF_ValueDX, LinF_ValueDY)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			LinF_ValueDX = 0;
			LinF_ValueDY = 0;

			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

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

#pragma omp parallel for firstprivate(LinF_ValueDX, LinF_ValueDY, LinF_ValuedXDY)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			LinF_ValueDX = 0;
			LinF_ValueDY = 0;
			LinF_ValuedXDY = 0;
			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;

				for (int dx = -rw; dx <= rw; dx++) {
					int KX = x + dx;
					if (KX < 0) KX = 0;
					if (KX > width - 1) KX = width - 1;

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

	t2.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	double k = 0.04;
#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double A = IDX[y][x];
			double B = IDY[y][x];
			double C = IDXY[y][x];
			R[y][x] = (A * B - C * C) - (k * ((A + B) * (A + B)));
		}
	}

	vector<Dot> harrisPoints;
	Dot tmp;
	double max = R[0][0];
#pragma omp parallel for
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			if (R[i][j] > max) max = R[i][j];

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (R[y][x] > max * 0.1) {
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

	t3.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
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

#pragma omp parallel for
	for (int i = 0; i < harrisPoints.size(); i++) {
		for (int r = -2; r < 2; r++) {
			for (int c = -2; c < 2; c++) {

				int yk = harrisPoints[i].Y + r;
				int xk = harrisPoints[i].X + c;
				if (yk >= 0 && xk >= 0 && yk < height && xk < width) {
					rgb_out[yk][xk].rgbtRed = 255;
					rgb_out[yk][xk].rgbtGreen = 0;
					rgb_out[yk][xk].rgbtBlue = 0;
				}
			}
		}
	}
	t7.push_back((omp_get_wtime() - start_time) * 1000);
	t8.push_back((omp_get_wtime() - start) * 1000);

	for (int i = 0; i < height; i++) {
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

void special_dots_omp_concurrent_vector(RGBTRIPLE**& rgb, RGBTRIPLE**& rgb_out, int height, int width, int ksize = 5)
{
	concurrent_vector <concurrent_vector <double>> I(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> IDX(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> IDY(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> IDXY(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> IDX0(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> IDY0(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> R(height, concurrent_vector <double>(width));
	concurrent_vector <concurrent_vector <double>> H(height, concurrent_vector <double>(width));

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++) {
			rgb_out[y][x].rgbtRed = rgb[y][x].rgbtRed;
			rgb_out[y][x].rgbtGreen = rgb[y][x].rgbtGreen;
			rgb_out[y][x].rgbtBlue = rgb[y][x].rgbtBlue;
		}

	double start = omp_get_wtime();
	double start_time = omp_get_wtime();
#pragma omp parallel for
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			I[y][x] = rgb[y][x].rgbtRed * 0.299 + rgb[y][x].rgbtGreen * 0.587 + rgb[y][x].rgbtBlue * 0.114;

	int size = ksize * ksize;
	int rh = ksize / 2;
	int rw = ksize / 2;
	int M = rh * 2 + 1;
	int N = rw * 2 + 1;
	double sigma = (double)ksize / 2;

	concurrent_vector <concurrent_vector <double>> MATR_coef(M, concurrent_vector <double>(N));

	double sum = 0;

#pragma omp parallel for reduction(+:sum)
	for (int y = -rh; y <= rh; y++) {
		for (int x = -rw; x <= rw; x++) {
			int yk = y + rh;
			int xk = x + rw;

			double CF = (1 / (2 * 3.14 * pow(sigma, 2))) * exp(-1 * ((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))));
			MATR_coef[yk][xk] = CF;
			sum += MATR_coef[yk][xk];
		}
	}

#pragma omp parallel for
	for (int y = -rh; y <= rh; y++) {
		for (int x = -rw; x <= rw; x++) {
			int yk = y + rh;
			int xk = x + rw;

			MATR_coef[yk][xk] /= sum;
		}
	}

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double LinF_Value = 0;

			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

					LinF_Value += double(I[ky][kx]) * MATR_coef[dy + rh][dx + rw];
				}
			}
			if (LinF_Value < 0) LinF_Value = 0;
			if (LinF_Value > 255) LinF_Value = 255;

			I[y][x] = LinF_Value;
		}
	}

	t1.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	//rw = 1;
	//rh = 1;
	double DiffX[3][3] = { {1, 0, -1}, {1, 0, -1}, {1, 0, -1} };
	double DiffY[3][3] = { {1, 1, 1}, {0, 0, 0}, {-1, -1, -1} };
	double LinF_ValueDX = 0;
	double LinF_ValueDY = 0;
	double LinF_ValuedXDY = 0;

#pragma omp parallel for firstprivate(LinF_ValueDX, LinF_ValueDY)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			LinF_ValueDX = 0;
			LinF_ValueDY = 0;

			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -rw; dx <= rw; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

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

#pragma omp parallel for firstprivate(LinF_ValueDX, LinF_ValueDY, LinF_ValuedXDY)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			LinF_ValueDX = 0;
			LinF_ValueDY = 0;
			LinF_ValuedXDY = 0;
			for (int dy = -rh; dy <= rh; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;

				for (int dx = -rw; dx <= rw; dx++) {
					int KX = x + dx;
					if (KX < 0) KX = 0;
					if (KX > width - 1) KX = width - 1;

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

	t2.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	double k = 0.04;
#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			double A = IDX[y][x];
			double B = IDY[y][x];
			double C = IDXY[y][x];
			R[y][x] = (A * B - C * C) - (k * ((A + B) * (A + B)));
		}
	}

	vector<Dot> harrisPoints;
	Dot tmp;
	double max = R[0][0];
#pragma omp parallel for
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			if (R[i][j] > max) max = R[i][j];

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (R[y][x] > max * 0.1) {
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

	t3.push_back((omp_get_wtime() - start_time) * 1000);
	start_time = omp_get_wtime();

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
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

#pragma omp parallel for
	for (int i = 0; i < harrisPoints.size(); i++) {
		for (int r = -2; r < 2; r++) {
			for (int c = -2; c < 2; c++) {

				int yk = harrisPoints[i].Y + r;
				int xk = harrisPoints[i].X + c;
				if (yk >= 0 && xk >= 0 && yk < height && xk < width) {
					rgb_out[yk][xk].rgbtRed = 255;
					rgb_out[yk][xk].rgbtGreen = 0;
					rgb_out[yk][xk].rgbtBlue = 0;
				}
			}
		}
	}
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
