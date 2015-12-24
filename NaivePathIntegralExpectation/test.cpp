#include <stdio.h>
#include <math.h>
#define N 400
#define nbrOfParticles 2
#define pathSize 200 //N nbrOfParticles
#define epsilon 0.1 f
#define epsiloninv 10.0 f
#define alpha 1.0 f

void xorshift(unsigned int *x, unsigned int *y, unsigned int *z, unsigned int *w)
{
	//This performs the necessary XORs to obtain a new random unsigned integer w.
	unsigned int t;
	t = *x ^ (*x << 11);
	*x = *y; *y = *z; *z = *w;
	*w = (*w ^ (*w >> 19) ^ (t ^ (t >> 8)));
}

float randFloat(unsigned int *x, unsigned int *y, unsigned int *z, unsigned int *w)
{
	xorshift(x, y, z, w);
	return (*w)*2.328306437080797e-10;
}

float oper(float *x)
{
	return 3.0f(x[0]  x[0] + x[N] x[N]);
}

float potential(float x)
{
	//This returns the value of the potential for a s p e c i f i e d point in time (
	index in[0, N1]) .
		return 1.5 f(x[0]  x[0] + x[N] x[N]);
}

float kineticDiff(float leftX, float midX, float rightX, float newMidX)
{
		float leftChangeOld = midX  leftX;
	float rightChangeOld = midX  rightX;
	float leftChangeNew = newMidX  leftX;
	float rightChangeNew = newMidX  rightX;
	return 0.5 f(leftChangeNewleftChangeNew + rightChangeNewrightChangeNew
		leftChangeOld  leftChangeOld  rightChangeOldrightChangeOld)  epsiloninv
		epsiloninv;
}




void metropolis(float path, float opMean, int  accepts, unsigned int  seeds,
	unsigned int nbrOfLoopings)
	f
	int modPathPoint = 0;
int modPoint;
int rightIndex;
int l e f t I n d e x;
float diffE;
float modx;
float oldX;
unsigned int loopLimit = nbrOfLoopings;
unsigned int threeLoopRuns = loopLimit / pathSize;

float partialOp = 0;
float opMeanBuffer = 0;

//This s e t s the seeds for the Xorshift PRNG.
unsigned int x = seeds[0];
unsigned int y = seeds[1];
unsigned int z = seeds[2];
unsigned int w = seeds[3];

for (int i = 0; i < N; i++)
	partialOp += oper(path + i);

opMean[0] = 0;
accepts[0] = 0;
for (int i = 0; i < threeLoopRuns; i++)
{
	for (int p a r t i c l e = 0; p a r t i c l e < nbrOfParticles; p a r t i c l e++)
	{
		for (modPoint = 0; modPoint < N; modPoint++)
		{
			modPathPoint = modPoint + N p a r t i c l e;

			//Arrange with the periodic boundary conditions .
			rightIndex = modPathPoint + 1;
			leftIndex = modPathPoint  1;
			if (modPoint == 0)
			{
				leftIndex = particleN + N  1;
				g else i f(modPoint == N  1)
					f
					rightIndex = particleN;
				g

					oldX = path[modPathPoint];
				modx = oldX + alpha(2.0 f  randFloat(&x, &y, &z, &w)  1.0 f);

				diffE = kineticDiff(path[l e f t I n d e x], path[modPathPoint], path[
					rightIndex], modx);




					diffE = potential(path + modPoint);
					path[modPathPoint] = modx;
					diffE += potential(path + modPoint);
					path[modPathPoint] = oldX;

					//Determine whether or not to accept the change in the path .
					i f(exp(epsilon  diffE) > randFloat(&x, &y, &z, &w)) f
						partialOp = oper(path + modPoint);
					path[modPathPoint] = modx;
					partialOp += oper(path + modPoint);
					accepts[0]++;
					g

						opMeanBuffer += partialOp;
					g
						g
						opMean[0] += opMeanBuffer;
					opMeanBuffer = 0;
					g

						opMean[0] += opMeanBuffer;
					opMean[0] /= (float)(threeLoopRuns nbrOfParticles NN);

					// Store the current s t a t e of the Xorshift PRNG for use in the next kernel run

					seeds[0] = x;
					seeds[1] = y;
					seeds[2] = z;
					seeds[3] = w;
			}
			//	132
			//	133 int main()
			//	134 f
			//	135     // I n i t i a l i z e  v a r i a b l e s
			//	136     float path[pathSize];
			//137     float opMean;
			//138     int accepts;
			//139    unsigned int seeds[4] = f123456789, 362436069, 521288629, 88675123g;
			//140    unsigned int nbrOfLoopings;
			//141     float sum = 0;
			//142     float squareSum = 0;
			//143     int runs = 10;
			//144     for (int i = 0; i < pathSize; i++)
			//	145            path[i] = 0;
			//146     nbrOfLoopings = 10000000;
			//147     //Run kernel once for burn in
			//	148     p r i n t f("Doing burn in . . . n n");
			//149     metropolis(path, &opMean, &accepts, seeds, &nbrOfLoopings);
			//150     //Run metropolis () runs times to produce runs values for the energy , each
			//	c a l l for metropolis() s t a r t s where the l a s t one ended(path and seeds are
			//	saved) .
			//	151     for (int i = 1; i <= runs; i++)
			//	152     f
			//	153            p r i n t f("Run %2d :n n", i);
			//154            metropolis(path, &opMean, &accepts, seeds, &nbrOfLoopings);
			//155            p r i n t f("Mean of oper : %f nn", opMean);
			//156            p r i n t f("Acceptance rate : %.2 f%%nnnn", 100.0 f((float)accepts) / ((float)
			//	((nbrOfLoopings / pathSize)  pathSize)));
			//157           sum += opMean;
			//158            squareSum += opMeanopMean;
			//159     g
			//	160     // Print r e s u l t s from runs c a l l s to metropolis () .

			//	E SOURCE CODE                                         101



			//	161         float mean = sum / ((float)runs);
			//162         float stdDev = sqrt(squareSum / ((float)runs)  meanmean);
			//163         float stdErr = stdDev / sqrt(runs);
			//164         p r i n t f("Mean : %f nn", mean);
			//165         p r i n t f("Standard error : %f nn", stdErr);
			//166         p r i n t f("65%% CI : [%f , %f ]n n", mean  stdErr, mean + stdErr);
			//167         p r i n t f("95%% CI : [%f , %f ]n n", mean  2 stdErr, mean + 2 stdErr);
			//168         return 0;
			//169 g