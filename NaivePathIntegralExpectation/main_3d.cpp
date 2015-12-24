#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include<chrono>
#include <algorithm>
#include <chrono>
#include <ppl.h>
#include <algorithm>
using namespace std;
#define PI 3.14159265359

mt19937 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());

//MAIN_3D: hydrogen atom

const int grid_width = 8;
const double grid_spacing = 5.0 / grid_width;
//const int num_samples = 3e5;
const int burn_in_samples = 1e4;
const int correlation_skip = 2; //because metropolis-hastings samples are correlated, skip some samples
const int tgrid_width = 128; //AKA Trotter number
const double time_range = 64;
double hbar = 1;
const double omega = 1;
const double lambda = 0;

const double epsilon = 0.001;

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

struct vec3 {
	double x, y, z;
	vec3() : x(0), y(0), z(0) {}
	vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

	double lenSquared() { return x*x + y*y + z*z; }
};


std::vector<double> sinTable(tgrid_width);

double computeSin(int x, int freq)
{
	return sinTable[(x*freq / 2) % tgrid_width];
}

void setUpSinTable()
{
	for (int i = 0; i < sinTable.size(); i++)
	{
		sinTable[i] = sin(2 * PI * i / sinTable.size());
	}
}

void printParameters()
{
	ofstream st(outputFilename, ios::out | ios::trunc);
	st << "grid_width: " << grid_width << endl;
	st << "grid_spacing: " << grid_spacing << endl;
	st << "correlation_skip: " << correlation_skip << endl;
	st << "burn_in_samples: " << burn_in_samples << endl;
}


void main()
{
	const int NUM_EVALS = 1;
	const int N = tgrid_width;
	vector<double> samples(NUM_EVALS);

	{
		ofstream st(outputFilename, ios::out | ios::trunc);
		st << "tgrid_width,energy,time_range,num_samples,acceptance rate,stddev,time taken" << endl;
	}

	auto t = chrono::high_resolution_clock::now();

	setUpSinTable();
	//for (int i = 0; i < tgrid_width; i++) cout << computeSin(i, 2) << endl;

	for (int k = 0; k < NUM_EVALS; k++)
	{
		//for (double time_range = 1.0; time_range < 128; time_range *= 2)
		{
			//for (double tgrid_width = 32; tgrid_width < 256; tgrid_width += 32) {
			for (int num_samples = 1e5; num_samples < 1e7; num_samples += 2e5)
			{
				auto runTime = chrono::high_resolution_clock::now();
				double aTotal = 0;
				int acnt = 0;
				const double dt = time_range / tgrid_width;

				double sumO = 0, sum = 0, sumSq = 0;

				//vector<double> zetas(N);

				int samplesCount = 0;
				concurrency::parallel_for((size_t)0, (size_t)grid_width, [&](int i) {
					/*for (int i = 0; i < grid_width;i++)
					{*/
					vector<vec3> path(tgrid_width);
					vector<vec3> newpath(tgrid_width);
					for (int t = 0; t < tgrid_width; t++) {
						path[t] = vec3(0.01, 0.01, 0.01); newpath[t] = path[t];
					}
					for (int y = 0; y < grid_width; y++)
					{
						for (int z = 0; z < grid_width; z++)
						{
							double prevEnergy = 0;
							vec3 xprime((2 * i - grid_width) * grid_spacing, (2 * y - grid_width) * grid_spacing, (2 * z - grid_width) * grid_spacing);
							for (int t = 0; t < tgrid_width; t++)
								path[t] = xprime;
							for (int j = 0; j < num_samples; j++)
							{
								for (int t = 0; t < tgrid_width; t++) newpath[t] = path[t];
								const double deltaWeight = 0.15;
								//int n = gen_uniform(0, tgrid_width - 1);
								//newpath[n].x += gen_uniform(-deltaWeight, deltaWeight);
								//newpath[n].y += gen_uniform(-deltaWeight, deltaWeight);
								//newpath[n].z += gen_uniform(-deltaWeight, deltaWeight);
								int chosenMode = gen_uniform(1, N + 1);
								vec3 randWeight = vec3(gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight));
								////zetas[chosenMode] += randWeight * chosenMode * PI;

								double SE = 0;
								for (int t = 0; t < tgrid_width; t++)
								{
									double s = computeSin(t, chosenMode);// sin(chosenMode * PI * t / tgrid_width);
									newpath[t].x = path[t].x + randWeight.x*s;
									newpath[t].y = path[t].y + randWeight.y*s;
									newpath[t].z = path[t].z + randWeight.z*s;
								}

								//for (double z : zetas) SE += z*z;

								//SE /= time_range * sqrt(time_range);

								const double minV = pow(tgrid_width, 2. / 3) / 10;

								for (int t = 1; t < tgrid_width - 1; t++)
								{
									vec3 dx = vec3(newpath[t + 1].x - newpath[t].x, newpath[t + 1].y - newpath[t].y, newpath[t + 1].z - newpath[t].z);
									SE += 0.5 * (dx.lenSquared()) / dt;
									double x2 = sqrt(newpath[t].lenSquared());
									//SE += 0.5 * omega * omega * newpath[t].lenSquared() * dt; // harmnonic oscillator
									//SE += std::max(-dt/(x2),-5.0); 
									double e = exp(-minV * x2);
									SE += -(1 - e) / x2;
								}

								if (prevEnergy == 0) prevEnergy = SE;

								//avgE += tgrid_width * time_range / 2;

								double a = std::min(exp((-SE + prevEnergy) / hbar), 1.);

								if (gen_uniform() < a) //Accept or not accept
								{
									path = newpath;
									prevEnergy = SE;
									if (j > burn_in_samples && j % correlation_skip == 0)
									{ //burn-in
										double avgE = 0;// omega* omega*newpath[tgrid_width / 2] * newpath[tgrid_width / 2];
										for (int t = 1; t < tgrid_width - 1; t++) //naive energy estimator
										{
											double x = sqrt(newpath[t].lenSquared());
											//avgE += omega * omega * newpath[t].lenSquared();
											//avgE -= 0.5/ (sqrt(newpath[t].lenSquared()));
											double e = exp(-minV * x);
											avgE += e * (-1.0 + (1 - 2 * x) / e - (-2 + minV) * x) / (2 * x * x);
										}
										avgE /= tgrid_width - 2;

										sumO += avgE;
										sum += 1;
										sumSq += avgE*avgE;

									}
								}
								acnt++;
								if (j % (80000) == 80000 - 1) cout << "a = " << (aTotal / acnt) << "\t" << sumO / (sum) << "\t\t sample " << j << endl;
								if (!std::isnan(a)) aTotal += a;
								//exp(-SE / hbar + 0.5 * lenSquared)
							}
							//#ifdef _WIN32
							//					system("cls");
							//#else
							//					system("clear");
							//#endif
							//cout << "\t\t" << sumO / sum << endl;
							//cout << "Done: " << (int)(100 * ((double)k / NUM_EVALS + (double)i / grid_width / NUM_EVALS)) << "%\n";
						}
					}
				});
				samples[k] = sumO / sum;

				//compute variance
				double stddev = sqrt((sumSq - (sumO * sumO) / (sum)) / (sum - 1));
				//change to stddev

				{
					ofstream st(outputFilename, ios::out | ios::app);
					st << tgrid_width 
						<< "," << sumO / sum 
						<< "," << time_range 
						<< "," << num_samples 
						<< "," << aTotal / acnt 
						<< "," << stddev 
						<< "," << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count()
						<< endl;
				}
			}
		}

	}
	double avgE = 0;
	double stdE = 0;
	for (double k : samples)
		avgE += k;
	avgE /= NUM_EVALS;
	for (double k : samples)
		stdE += (k - avgE)*(k - avgE);
	stdE = sqrt(stdE / NUM_EVALS);

	for (double k : samples) cout << k << endl;
	cout << endl;
	cout << "Final energy:" << avgE << " +/- " << stdE << endl;
	cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count() << endl;
	while (1);

}