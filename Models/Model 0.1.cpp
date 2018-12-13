#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include <exception>

using namespace std;

const double b = 4.0;							//Benefits
const double c = 2.0;							//Costs
const double mu = 0.01;							//Mutation rate
const double sigma = 0.01;						//Stdev of change after mutation
const int n = 1000;								//Population size
const int nPopulations = 3;						//Number of populations
const int nInteractions = 10;					//Number of interactions in individual's life
const int nGenerations = 10000;					//Number of generations
const int nGenSav = 100;						//Save every n generations
const vector<double> Pc = {0.95, 0.05, 0.67};	//Initial mean tendency of populations to cooperate


mt19937_64 rng;

struct Individual
{
	double strategy;		//0 always uses D; 1 always uses C
	double fitness;			//Cumulative payoff values
};

int main()
{
	try {
		//Seed rng
		rng.seed(1);

		//Open output file
		ofstream ofs("Model 0.1.csv");
		if (!ofs.is_open()) {
			throw logic_error("Unable to open output file\n");
		}

		//More checks
		if (nPopulations != Pc.size()) {
			throw logic_error("Number of populations and number of initial mean Pc values are different\n");
		}

		//Initialise population
		vector<vector<Individual> > Populations(nPopulations, vector<Individual>(n));
		cout << "0";
		ofs << "0";
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			double mean = 0.0;
			for (int i = 0; i < n; ++i) {
				Populations[nPop][i].fitness = 0.0;
				Populations[nPop][i].strategy = Pc[nPop];
				mean += Populations[nPop][i].strategy;
			}
			mean /= n;
			double stdev = 0.0;
			for (int i = 0; i < n; ++i) {
				stdev += pow(Populations[nPop][i].strategy - mean, 2);
			}
			stdev = sqrt(stdev / n);
			cout << "\t" << mean << "\t" << stdev;
			ofs << "\t" << mean << "\t" << stdev;
		}
		cout << "\n";
		ofs << "\n";

		//Simulate
		uniform_real_distribution<double> chooseFraction(0.0, 1.0);
		uniform_int_distribution<int> pickPartner(0, n - 1);
		for (int g = 1; g <= nGenerations; ++g) {
			if (g % nGenSav == 0) {
				cout << g;
				ofs << g;
			}
			for (int nPop = 0; nPop < nPopulations; ++nPop) {

				//Interactions
				for (int i = 0; i < n; ++i) {
					for (int k = 0; k < nInteractions; ++k) {
						int j = pickPartner(rng);
						double r = chooseFraction(rng);
						double s = chooseFraction(rng);
						if (Populations[nPop][i].strategy > r && Populations[nPop][j].strategy > s) {		//Both cooperate
							Populations[nPop][i].fitness += b - c / 2;
						}
						else if (Populations[nPop][i].strategy > r && Populations[nPop][j].strategy < s) {	//Only focal individual cooperates
							Populations[nPop][i].fitness += b - c;
						}
						else if (Populations[nPop][i].strategy < r && Populations[nPop][j].strategy > s) {	//Only focal individual defects
							Populations[nPop][i].fitness += b;
						}
						else if (Populations[nPop][i].strategy < r && Populations[nPop][j].strategy < s) {	//Both defect
							Populations[nPop][i].fitness += 0;
						}
					}
				}

				//Determine offspring
				vector<double> vecWeights(n);
				for (int i = 0; i < n; ++i) {
					vecWeights[i] = Populations[nPop][i].fitness;
				}
				vector<Individual> PopulationNew(n);
				discrete_distribution<int> chooseParent(vecWeights.begin(), vecWeights.end());
				for (int i = 0; i < n; ++i) {
					PopulationNew[i] = Populations[nPop][chooseParent(rng)];
					PopulationNew[i].fitness = 0.0;
					if (chooseFraction(rng) < mu) {
						normal_distribution<double> defineMutation(0.0, sigma);
						PopulationNew[i].strategy += defineMutation(rng);
						if (PopulationNew[i].strategy < 0.0) {
							PopulationNew[i].strategy = 0.0;
						}
						else if (PopulationNew[i].strategy > 1.0) {
							PopulationNew[i].strategy = 1.0;
						}
					}
				}
				Populations[nPop] = PopulationNew;

				if (g % nGenSav == 0) {
					//Calculate mean Pc
					double mean = 0.0;
					for (int i = 0; i < n; ++i) {
						mean += Populations[nPop][i].strategy;
					}
					mean /= n;

					//Calculate stdev of Pc
					double stdev = 0.0;
					for (int i = 0; i < n; ++i) {
						stdev += pow(Populations[nPop][i].strategy - mean, 2);
					}
					stdev = sqrt(stdev / n);

					//Output
					cout << "\t" << mean << "\t" << stdev;
					ofs << "\t" << mean << "\t" << stdev;
				}
			}
			if (g % nGenSav == 0) {
				cout << "\n";
				ofs << "\n";
			}
		}
		ofs.close();
	}
	catch (exception &error) {
		cerr << "Error: " << error.what();
		return 1;
	}

	return 0;
}