#include <iostream>
#include <fstream>
#include <random>
#include <vector>

using namespace std;

const double b = 4.0;		//Benefits
const double c = 2.0;		//Costs
const double mu = 0.01;		//Mutation rate
const int n = 1000;			//Population size
const double Pc1 = 0.95;	//Initial fraction of population 1 with C
const double Pc2 = 0.05;	//Initial fraction of population 2 with C
const int ni = 10;			//Number of interactions in individual's life
const int ng = 100;			//Number of generation

mt19937_64 rng;

struct Individual
{
	int strategy;			//0 always uses D; 1 always uses C
	double fitness;			//Cumulative payoff values
};

int main()
{
	//Seed rng
	rng.seed(1);

	//Open output file
	ofstream ofs("Model 0.0.csv");
	if (!ofs.is_open()) {
		cout << "Error: Unable to open output file\n";
		return 1;
	}

	//Initialise population
	vector<Individual> vecPopulation1(n), vecPopulation2(n);
	double k1 = 0.0, k2 = 0.0;
	for (int i = 0; i < n; ++i) {
		vecPopulation1[i].fitness = 0.0;
		vecPopulation2[i].fitness = 0.0;
		vecPopulation1[i].strategy = i < n * Pc1 ? 1 : 0;
		vecPopulation2[i].strategy = i < n * Pc2 ? 1 : 0;
		if (vecPopulation1[i].strategy == 1) {
			k1 += 1.0;
		}
		if (vecPopulation2[i].strategy == 1) {
			k2 += 1.0;
		}
	}
	cout << "0\t" << k1 / n << "\t" << k2 / n << "\n";
	ofs << "0\t" << k1 / n << "\t" << k2 / n << "\n";

	//Simulate
	for (int g = 1; g <= ng; ++g) {

		//Interactions
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < ni; ++k) {
				uniform_int_distribution<int> pickPartner(0, n - 1);
				int j = pickPartner(rng);
				if (vecPopulation1[i].strategy == 1 && vecPopulation1[j].strategy == 1) {		//Both cooperate
					vecPopulation1[i].fitness += b - c / 2;
				}
				else if (vecPopulation1[i].strategy == 1 && vecPopulation1[j].strategy == 0) {	//Only focal individual cooperates
					vecPopulation1[i].fitness += b - c;
				}
				else if (vecPopulation1[i].strategy == 0 && vecPopulation1[j].strategy == 1) {	//Only focal individual defects
					vecPopulation1[i].fitness += b;
				}
				else if (vecPopulation1[i].strategy == 0 && vecPopulation1[j].strategy == 0) {	//Both defect
					vecPopulation1[i].fitness += 0;
				}

				if (vecPopulation2[i].strategy == 1 && vecPopulation2[j].strategy == 1) {		//Both cooperate
					vecPopulation2[i].fitness += b - c / 2;
				}
				else if (vecPopulation2[i].strategy == 1 && vecPopulation2[j].strategy == 0) {	//Only focal individual cooperates
					vecPopulation2[i].fitness += b - c;
				}
				else if (vecPopulation2[i].strategy == 0 && vecPopulation2[j].strategy == 1) {	//Only focal individual defects
					vecPopulation2[i].fitness += b;
				}
				else if (vecPopulation2[i].strategy == 0 && vecPopulation2[j].strategy == 0) {	//Both defect
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
				vecPopulationNew1[i].strategy = vecPopulationNew1[i].strategy == 1 ? 0 : 1;
				vecPopulationNew2[i].strategy = vecPopulationNew2[i].strategy == 1 ? 0 : 1;
			}
		}
		vecPopulation1 = vecPopulationNew1;
		vecPopulation2 = vecPopulationNew2;

		//Calculate Pc
		k1 = k2 = 0.0;
		for (int i = 0; i < n; ++i) {
			if (vecPopulation1[i].strategy == 1) {
				k1 += 1.0;
			}
			if (vecPopulation2[i].strategy == 1) {
				k2 += 1.0;
			}
		}

		//Output
		cout << g << "\t" << k1 / n << "\t" << k2 / n << "\n";
		ofs << g << "\t" << k1 / n << "\t" << k2 / n << "\n";
	}
	
	return 0;
}