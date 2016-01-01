#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include<chrono>
#include <algorithm>
#include <chrono>
#include <ppl.h>
#include <algorithm>
#include<array>

#include <glm\glm.hpp>
#include <glm\gtx\norm.hpp>
#include "Nucleus.h"
#include "SystemState.h"

using namespace std;
#define PI 3.14159265359

mt19937 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());

//MAIN_3D: hydrogen atom

//For ground state energy, no need to sample at same grid points
//const int grid_width = 10;
//const double grid_spacing = 8.0 / grid_width;
//const int num_samples = 3e5;
const double box_bounds = 1;
const int burn_in_samples = 1e5;
const int correlation_skip = 2; //because metropolis-hastings samples are correlated, skip some samples
const int tgrid_width = 300; //AKA Trotter number
const double time_range = 50;
double hbar = 1;
const double omega = 1;
const double lambda = 0;

const int NUM_SAMPLE_POINTS = 300; // multiplies number of samples

const bool BIASED_UPDATES = false;
const bool ADAPTIVE_STEPS = false;
const double ADAPTIVE_STEP_RANGE = 0.05;
const double ADAPTIVE_STEP_TARGET = 0.3;
const double ADAPTIVE_STEP_CHANGE = 0.01;

const double MINV_MULT = 3 * pow(time_range, -2./3); //the prefactor according to Muser and Berne 1997

const double epsilon = 0.001;

const int NUM_SAMPLES_MIN = 1e6;
const int NUM_SAMPLES_MAX = 1e6;
const int NUM_SAMPLES_INC = 2e5;


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
	//st << "grid_width: " << grid_width << endl;
	//st << "grid_spacing: " << grid_spacing << endl;
	st << "correlation_skip: " << correlation_skip << endl;
	st << "burn_in_samples: " << burn_in_samples << endl;
}


vector<Nucleus> nuclei;
void initializeNuclei()
{
	const double bondLength = 2.0;
	//hydrogen molecule
	nuclei.push_back(Nucleus(glm::dvec3(0, 0, 0)));
	//nuclei.push_back(Nucleus(glm::dvec3(-2 / 2, 0, 0)));
	//nuclei.push_back(Nucleus(glm::dvec3(2 / 2, 0, 0)));
}

double effectiveCoulomb(double x, double minV, double e) { //where e = exp(-minV*x) for simplification
	return -(1 - e) / x;
}
double effectiveCoulombD1(double x, double minV, double e) { //first derivative
	return (1 - e * (1 + minV * x)) / (x*x);
}
double effectiveCoulombD2(double x, double minV, double e) { //second derivative
	return e * (2.0 - 2.0 / e + minV * x * (2 + minV * x)) / (x*x*x);
}
double effectiveCoulombD3(double x, double minV, double e) {//third derivative
	return e * (-6.0 + 6.0 / e - minV * x * (6 + minV * x * (3 + minV * x))) / (x*x*x*x);
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
	initializeNuclei();
	//for (int i = 0; i < tgrid_width; i++) cout << computeSin(i, 2) << endl;

	for (int k = 0; k < NUM_EVALS; k++)
	{
		//for (double time_range = 1.0; time_range < 128; time_range *= 2)
		{
			//for (double tgrid_width = 32; tgrid_width < 256; tgrid_width += 32) {
			for (int num_samples = NUM_SAMPLES_MIN; num_samples <= NUM_SAMPLES_MAX; num_samples += NUM_SAMPLES_INC)
			{
				auto runTime = chrono::high_resolution_clock::now();
				double aTotal = 0;
				double aPrev = 0;
				int acnt = 0;
				const double dt = time_range / tgrid_width;

				double sumO = 0, sum = 0, sumSq = 0;

				//vector<double> zetas(N);

				int samplesCount = 0;
				concurrency::parallel_for((size_t)0, (size_t)NUM_SAMPLE_POINTS, [&](int i) {
					double deltaWeight = 0.1;
					/*for (int i = 0; i < grid_width;i++)
					{*/
					vector<SystemState> path(tgrid_width);
					vector<SystemState> newpath(tgrid_width);
					for (int t = 0; t < tgrid_width; t++) {
						for (glm::dvec3& pos : path[t].positions)
							pos = glm::dvec3(0.01, 0.01, 0.01);
						newpath[t] = path[t];
					}
					double prevEnergy = 0;
					//random xprime
					SystemState xprime;
					for (glm::dvec3& pos : xprime.positions)
						pos = glm::dvec3(gen_uniform(-box_bounds, box_bounds), gen_uniform(-box_bounds, box_bounds), gen_uniform(-box_bounds, box_bounds));
					for (int t = 0; t < tgrid_width; t++)
						path[t] = xprime;
					for (int j = 0; j < num_samples; j++)
					{
						for (int t = 0; t < tgrid_width; t++) newpath[t] = path[t];
						if (ADAPTIVE_STEPS && j < burn_in_samples) // only use adaptive steps during burn-in period
						{
							if (aTotal / acnt > ADAPTIVE_STEP_TARGET + ADAPTIVE_STEP_RANGE) deltaWeight += ADAPTIVE_STEP_CHANGE;
							if (aTotal / acnt < ADAPTIVE_STEP_TARGET + ADAPTIVE_STEP_RANGE) deltaWeight -= ADAPTIVE_STEP_CHANGE;
							if (deltaWeight < 0) deltaWeight *= -1;
							//cout << deltaWeight<<endl;
						}/*
						int n = gen_uniform(0, tgrid_width - 1);

						for (glm::dvec3& pos : newpath[n].positions)
							pos += glm::dvec3(gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight));*/
						for (int np = 0; np < NUM_PARTICLES; np++)
						{
							int chosenMode = gen_uniform(1, N + 1);
							glm::dvec3 randWeight = glm::dvec3(gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight), gen_uniform(-deltaWeight, deltaWeight));

							for (int t = 0; t < tgrid_width; t++)
							{
								double s = sin(chosenMode * PI * t / tgrid_width);
								if (BIASED_UPDATES) s /= chosenMode;
								newpath[t].positions[np] = path[t].positions[np] + randWeight*s;
							}
						}
						double SE = 0;

						const double minV = pow(tgrid_width, 2. / 3) * MINV_MULT;

						for (int np = 0; np < NUM_PARTICLES; np++)
						{
							for (int t = 0; t < tgrid_width - 1; t++)
							{
								//Calculate KE
								glm::dvec3 dx = newpath[t + 1].positions[np] - newpath[t].positions[np];
								SE += 0.5 * (glm::length2(dx)) / dt;
								//Calculate PE
								for (Nucleus n : nuclei) {
									double x = glm::length(newpath[t].positions[np] - n.pos);
									double e = exp(-minV * x);
									SE += -n.charge * (1 - e) / x * dt;
								}
							}
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
								double avgE = 0;
								for (int np = 0; np < NUM_PARTICLES; np++)
								{
									for (int t = 0; t < tgrid_width - 1; t++) //naive energy estimator
									{
										double V0 = 0, V1 = 0, V2 = 0, V3 = 0;//derivatives of the potential, with relevant weights
										double dx2 = glm::length2((newpath[t + 1].positions[np] - newpath[t].positions[np]));

										//calculate potential and derivatives of potential for each charge source

										for (Nucleus n : nuclei) {
											double x = glm::length(newpath[t].positions[np] - n.pos);
											//avgE += omega * omega * newpath[t].lenSquared();
											//avgE -= 0.5/ (sqrt(newpath[t].lenSquared()));

											double e = exp(-minV * x);
											//avgE += e * (-1.0 + (1 - 2 * x) / e - (-2 + minV) * x) / (2 * x * x); //first order energy estimator
											//double inc = -0.5 / x - 0.5 * e * minV + 0.5 * e / x;//first order energy estimator
											//if (!std::isnan(inc)) avgE += inc;

											//energy estimators from Grujic 2006 

											double
												i0 = effectiveCoulomb(x, minV, e),
												i1 = x * effectiveCoulombD1(x, minV, e) / 2,
												i2 = (dt / 6 + dx2 / 12) * effectiveCoulombD2(x, minV, e),
												i3 = (x * dt / 24 + x * dx2 / 48) * effectiveCoulombD3(x, minV, e);

											V0 += i0;
											V1 += i1;

											/*if (!std::isnan(i0)) V0 += i0;
											if (!std::isnan(i1)) V1 += i1;*/
											/*if (!std::isnan(i2)) V2 += i2;
											if (!std::isnan(i3)) V3 += i3;*/

											//double inc = e * (
											//	-24.0 * x * x * (-1.0 + 1.0 / e + minV * x)
											//	+ dt * (4.0 - 4.0 / e + 4.0 * minV * x
											//	+ 2.0 * minV * minV * x * x
											//	- 2.0 * minV * minV * minV * x * x * x)
											//	+ dx2 * (2.0 - 2.0 / e + 2.0 * minV * x + minV * minV * x * x - minV * minV * minV * x * x * x)) / (48.0 * x * x * x);//second order energy estimator according to Grujic 2006
											//if (!std::isnan(inc)) avgE += inc;
										}
										for (int ni = 0; ni < NUM_PARTICLES; ni++)
										{
											if (ni == np) continue; //no self-interaction
											double x = glm::length(newpath[t].positions[np] - newpath[t].positions[ni]);
											double e = exp(-minV * x);
											//energy estimators from Grujic 2006 
											V0 += effectiveCoulomb(x, minV, e);
											V1 += x * effectiveCoulombD1(x, minV, e) / 2;
											V2 += (dt / 6 + dx2 / 12) * effectiveCoulombD2(x, minV, e);
											V3 += (x * dt / 24 + x * dx2 / 48) * effectiveCoulombD3(x, minV, e);
										}
										avgE += V0;
										avgE += V1;
										avgE += V2;
										avgE += V3;
									}


								}
								avgE /= tgrid_width;

								sumO += avgE;
								sum++;
								sumSq += avgE*avgE;

							}
						}
						acnt++;
						if (j == num_samples - 1) cout << "a = " << (aTotal / acnt) << "\t" << sumO / (sum) << "\t\t sample " << j << endl;
						if (!std::isnan(a)) aTotal += a;
						aPrev = a;
						//exp(-SE / hbar + 0.5 * lenSquared)
					}
					//#ifdef _WIN32
					//					system("cls");
					//#else
					//					system("clear");
					//#endif
					//cout << "\t\t" << sumO / sum << endl;
					//cout << "Done: " << (int)(100 * ((double)k / NUM_EVALS + (double)i / grid_width / NUM_EVALS)) << "%\n";
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