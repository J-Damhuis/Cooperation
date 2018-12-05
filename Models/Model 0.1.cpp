#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>

using namespace std;

const double b = 4.0;		//Benefits
const double c = 2.0;		//Costs
const double mu = 0.01;		//Mutation rate
const double sigma = 0.01;	//Stdev of change after mutation
const double Pc1 = 0.9;		//Initial mean Pc of population 1
const double Pc2 = 0.1;		//Initial mean Pc of population 2
const double Pcsigma = 0.3;	//Initial stdev of Pc
const int n = 1000;			//Population size
const int ni = 10;			//Number of interactions in individual's life
const int ng = 100;			//Number of generations

mt19937_64 rng;

struct Individual
{
	double strategy;		//0 always uses D; 1 always uses C
	double fitness;			//Cumulative payoff values
};

int main()
{
	//Seed rng
	rng.seed(1);

	//Open output file
	ofstream ofs("Model 0.1.csv");
	if (!ofs.is_open()) {
		cout << "Error: Unable to open output file\n";
		return 1;
	}

	//Initialise population
	vector<Individual> vecPopulation1(n), vecPopulation2(n);
	uniform_real_distribution<double> choosePc(0.0, 1.0);
	normal_distribution<double> chooseInitialPc1(Pc1, Pcsigma);
	normal_distribution<double> chooseInitialPc2(Pc2, Pcsigma);
	double mean1 = 0.0, mean2 = 0.0;
	for (int i = 0; i < n; ++i) {
		vecPopulation1[i].fitness = 0.0;
		vecPopulation2[i].fitness = 0.0;
		vecPopulation1[i].strategy = chooseInitialPc1(rng);
		vecPopulation2[i].strategy = chooseInitialPc2(rng);
		if (vecPopulation1[i].strategy < 0.0) {
			vecPopulation1[i].strategy = 0.0;
		}
		else if (vecPopulation1[i].strategy > 1.0) {
			vecPopulation1[i].strategy = 1.0;
		}
		if (vecPopulation2[i].strategy < 0.0) {
			vecPopulation2[i].strategy = 0.0;
		}
		else if (vecPopulation2[i].strategy > 1.0) {
			vecPopulation2[i].strategy = 1.0;
		}
		mean1 += vecPopulation1[i].strategy;
		mean2 += vecPopulation2[i].strategy;
	}
	mean1 /= n;
	mean2 /= n;
	double stdev1 = 0.0, stdev2 = 0.0;
	for (int i = 0; i < n; ++i) {
		stdev1 += pow(vecPopulation1[i].strategy - mean1, 2);
		stdev2 += pow(vecPopulation2[i].strategy - mean2, 2);
	}
	stdev1 = sqrt(stdev1 / n);
	stdev2 = sqrt(stdev2 / n);
	cout << "0\t" << mean1 << "\t" << stdev1 << "\t" << mean2 << "\t" << stdev2 << "\n";
	ofs << "0\t" << mean1 << "\t" << stdev1 << "\t" << mean2 << "\t" << stdev2 << "\n";

	//Simulate
	for (int g = 1; g <= ng; ++g) {

		//Interactions
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < ni; ++k) {
				uniform_int_distribution<int> pickPartner(0, n - 1);
				int j = pickPartner(rng);
				double r = choosePc(rng);
				double s = choosePc(rng);
				if (vecPopulation1[i].strategy > r && vecPopulation1[j].strategy > s) {		//Both cooperate
					vecPopulation1[i].fitness += b - c / 2;
				}
				else if (vecPopulation1[i].strategy > r && vecPopulation1[j].strategy < s) {	//Only focal individual cooperates
					vecPopulation1[i].fitness += b - c;
				}
				else if (vecPopulation1[i].strategy < r && vecPopulation1[j].strategy > s) {	//Only focal individual defects
					vecPopulation1[i].fitness += b;
				}
				else if (vecPopulation1[i].strategy < r && vecPopulation1[j].strategy < s) {	//Both defect
					vecPopulation1[i].fitness += 0;
				}

				if (vecPopulation2[i].strategy > r && vecPopulation2[j].strategy > s) {		//Both cooperate
					vecPopulation2[i].fitness += b - c / 2;
				}
				else if (vecPopulation2[i].strategy > r && vecPopulation2[j].strategy < s) {	//Only focal individual cooperates
					vecPopulation2[i].fitness += b - c;
				}
				else if (vecPopulation2[i].strategy < r && vecPopulation2[j].strategy > s) {	//Only focal individual defects
					vecPopulation2[i].fitness += b;
				}
				else if (vecPopulation2[i].strategy < r && vecPopulation2[j].strategy < s) {	//Both defect
					vecPopulation2[i].fitness += 0;
				}
			}
		}

		//Determine offspring
		vector<double> vecWeights1(n), vecWeights2(n);
		for (int i = 0; i < n; ++i) {
			vecWeights1[i] = vecPopulation1[i].fitness;
			vecWeights2[i] = vecPopulation2[i].fitness;
		}
		vector<Individual> vecPopulationNew1(n), vecPopulationNew2(n);
		discrete_distribution<int> chooseParent1(vecWeights1.begin(), vecWeights1.end());
		discrete_distribution<int> chooseParent2(vecWeights2.begin(), vecWeights2.end());
		vector<double> vecMutation = {mu, 1 - mu};
		discrete_distribution<int> chooseMutation(vecMutation.begin(), vecMutation.end());
		for (int i = 0; i < n; ++i) {
			vecPopulationNew1[i] = vecPopulation1[chooseParent1(rng)];
			vecPopulationNew2[i] = vecPopulation2[chooseParent2(rng)];
			vecPopulationNew1[i].fitness = 0.0;
			vecPopulationNew2[i].fitness = 0.0;
			if (chooseMutation(rng) == 0) {
				normal_distribution<double> defineMutation(0.0, sigma);
				vecPopulationNew1[i].strategy += defineMutation(rng);
				vecPopulationNew2[i].strategy += defineMutation(rng);
				if (vecPopulationNew1[i].strategy < 0.0) {
					vecPopulationNew1[i].strategy = 0.0;
				}
				else if (vecPopulationNew1[i].strategy > 1.0) {
					vecPopulationNew1[i].strategy = 1.0;
				}
				if (vecPopulationNew2[i].strategy < 0.0) {
					vecPopulationNew2[i].strategy = 0.0;
				}
				else if (vecPopulationNew2[i].strategy > 1.0) {
					vecPopulationNew2[i].strategy = 1.0;
				}
			}
		}
		vecPopulation1 = vecPopulationNew1;
		vecPopulation2 = vecPopulationNew2;

		//Calculate mean Pc
		mean1 = 0.0;
		mean2 = 0.0;
		for (int i = 0; i < n; ++i) {
			mean1 += vecPopulation1[i].strategy;
			mean2 += vecPopulation2[i].strategy;
		}
		mean1 /= n;
		mean2 /= n;

		//Calculate stdev of Pc
		stdev1 = 0.0;
		stdev2 = 0.0;
		for (int i = 0; i < n; ++i) {
			stdev1 += pow(vecPopulation1[i].strategy - mean1, 2);
			stdev2 += pow(vecPopulation2[i].strategy - mean2, 2);
		}
		stdev1 = sqrt(stdev1 / n);
		stdev2 = sqrt(stdev2 / n);

		//Output
		cout << g << "\t" << mean1 << "\t" << stdev1 << "\t" << mean2 << "\t" << stdev2 << "\n";
		ofs << g << "\t" << mean1 << "\t" << stdev1 << "\t" << mean2 << "\t" << stdev2 << "\n";
	}
	ofs.close();

	return 0;
}