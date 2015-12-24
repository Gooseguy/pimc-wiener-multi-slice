#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include<chrono>
#include <algorithm>
#include <chrono>
#include <ppl.h>
using namespace std;
#define PI 3.14159265359

mt19937 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());


const int grid_width = 128;
const double grid_spacing = 6.0 / grid_width;
const int num_samples = 5e6;
const int burn_in_samples = 1e6;
const int correlation_skip = 5; //because metropolis-hastings samples are correlated, skip some samples
const int tgrid_width = 128;
const double time_range = 32;
double hbar = 1;
const double omega = 1;
const double lambda = 0.3;

double gen_norm(double mean = 0, double stddev = 1) {
	normal_distribution<double> dist(mean, stddev);
	return dist(mt_rand);
}
double gen_uniform(double min = 0, double max = 1) {
	uniform_real_distribution<double> dist(min, max);
	return dist(mt_rand);
}

int gen_uniform(int min, int max) {
	uniform_int_distribution<int> dist(min, max);
	return dist(mt_rand);
}

const string outputFilename = "C:\\Users\\Christian\\Desktop\\path_integrals.csv";

void main()
{
	const int NUM_EVALS = 1;
	const int N = 256;
	vector<double> samples(NUM_EVALS);

	{
		ofstream st(outputFilename, ios::out | ios::trunc);
		st << "tgrid_width,energy" << endl;
	}

	auto t = chrono::high_resolution_clock::now();


	for (int k = 0; k < NUM_EVALS; k++)
	{

		//for (double tgrid_width = 16; tgrid_width < 256; tgrid_width*=2) {
			const double dt = time_range / tgrid_width;

			double sumO = 0, sum = 0;

			//vector<double> zetas(N);

			int samplesCount = 0;
			concurrency::parallel_for((size_t)0, (size_t)grid_width, [&](int i) {

				vector<double> path(tgrid_width);
				vector<double> newpath(tgrid_width);
					{
						double prevEnergy = 0;
						double aTotal = 0;
						double xprime = (2 * i - grid_width) * grid_spacing;
						for (int t = 0; t < tgrid_width; t++)
							path[t] = xprime;
						for (int j = 0; j < num_samples; j++)
						{

							int chosenMode = gen_uniform(1, N + 1);
							double randWeight = gen_uniform(-0.2, 0.2);
							//zetas[chosenMode] += randWeight * chosenMode * PI;

							double SE = 0;
							//double l = 0;
							for (int t = 0; t < tgrid_width; t++)
								newpath[t] = path[t] + randWeight*sin(chosenMode * PI * (double)t / tgrid_width);

							//for (double z : zetas) SE += z*z;

							//SE /= time_range * sqrt(time_range);

							for (int t = 1; t < tgrid_width - 1; t++)
							{
								SE += 0.5 * (newpath[t + 1] - newpath[t])*(newpath[t + 1] - newpath[t]) / dt;
								SE += (0.5 * omega * omega* newpath[t] * newpath[t] + lambda * newpath[t] * newpath[t] * newpath[t] * newpath[t]) * dt;
							}

							if (prevEnergy == 0) prevEnergy = SE;

							double avgE = 0;// omega* omega*newpath[tgrid_width / 2] * newpath[tgrid_width / 2];
							for (int t = 1; t < tgrid_width - 1; t++)
							{
								double x12 = newpath[t];
								double delta = newpath[t + 1] - newpath[t];
								avgE += omega * omega*x12 * x12 + 3 * lambda * x12 *x12 *x12 *x12;
							}
							avgE /= tgrid_width-2;

							double lenSquared = 0;
							for (double z : path) lenSquared += (z - xprime)*(z - xprime);

							double a = std::min(exp((-SE + prevEnergy) / hbar), 1.);

							if (gen_uniform() < a)
							{
								path = newpath;
								prevEnergy = SE;
								if (j > burn_in_samples && j % correlation_skip == 0)
								{ //burn-in
									sumO += avgE;
									sum += 1;
								}
							}
							aTotal += a;
							if (j % (100000) == 1) cout << "a = " << aTotal / j << "\t\t" << sumO / sum << endl;

							//exp(-SE / hbar + 0.5 * lenSquared)
						}
#ifdef _WIN32
						system("cls");
#else
						system("clear");
#endif
						cout << "\t\t" << sumO / sum << endl;
						//cout << "Done: " << (int)(100 * ((double)k / NUM_EVALS + (double)i / grid_width / NUM_EVALS)) << "%\n";
					}
			});
			samples[k] = sumO / sum;

			{
				ofstream st(outputFilename, ios::out | ios::app);
				st << tgrid_width << "," << sumO / sum << endl;
			}
		
	}
	double avgE = 0;
	double stdE = 0;
	for (double k : samples)
		avgE += k;
	avgE /= NUM_EVALS;
	for (double k : samples)
		stdE += (k - avgE)*(k - avgE);
	stdE = sqrt(stdE/NUM_EVALS);

	for (double k : samples) cout << k << endl;
	cout << endl;
	cout << "Final energy:" << avgE << " +/- " << stdE << endl;
	cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count() << endl;
	while (1);

}