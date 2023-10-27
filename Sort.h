#pragma once

int increment(long inc[], long size) {
	int p1, p2, p3, s;
	p1 = p2 = p3 = 1;
	s = -1;
	do {
		if (++s % 2) {
			inc[s] = 8 * p1 - 6 * p2 + 1;
		}
		else {
			inc[s] = 9 * p1 - 9 * p3 + 1;
			p2 *= 2;
			p3 *= 2;
		}
		p1 *= 2;
	} while (3 * inc[s] < size);
	return s > 0 ? --s : 0;
}

vector<Dot> shell(size_t size, vector<Dot>arr) {
	long inc, j, seq[40];
	int s;
	s = increment(seq, size);
	while (s >= 0) {
		inc = seq[s--];
		for (int i = inc; i < size; i++) {
			Dot temp = arr[i];
			for (j = i - inc; (j >= 0) && (arr[j].R > temp.R); j -= inc)
				arr[j + inc] = arr[j];
			arr[j + inc] = temp;
		}
	}
	return arr;
}