#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <exception>

using namespace std;

const double b = 4.0;							//Benefits
const double c = 2.0;							//Costs
const double mu = 0.01;							//Mutation rate
const int n = 1000;								//Population size
const int nPopulations = 3;						//Number of populations
const int nInteractions = 10;					//Number of interactions in individual's life
const int nGenerations = 100;					//Number of generation
const vector<double> Pc = {0.95, 0.05, 0.67};	//Initial fractions of populations that cooperates

mt19937_64 rng;

struct Individual
{
	int strategy;			//0 always uses D; 1 always uses C
	double fitness;			//Cumulative payoff values
};

int main()
{
	try {
		//Seed rng
		rng.seed(1);

		//Open output file
		ofstream ofs("Model 0.0.csv");
		if (!ofs.is_open()) {
			throw logic_error("Unable to open output file\n");
		}

		//More checks
		if (nPopulations != Pc.size()) {
			throw logic_error("Number of populations and number of initial mean Pc values are different\n");
		}

		//Initialise populations and output initial state
		vector<vector<Individual> > Populations(nPopulations, vector<Individual>(n));
		cout << "Generation";
		ofs << "Generation";
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			cout << "\tMeanP0";
			ofs << "\tMeanP0";
		}
		cout << "\n0";
		ofs << "\n0";
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			double mean = 0.0;
			for (int i = 0; i < n; ++i) {
				Populations[nPop][i].fitness = 0.0;
				Populations[nPop][i].strategy = i < n * Pc[nPop] ? 1 : 0;
				mean += Populations[nPop][i].strategy;
			}
			cout << "\t" << mean / n;
			ofs << "\t" << mean / n;
		}
		cout << "\n";
		ofs << "\n";

		//Simulate
		uniform_int_distribution<int> pickPartner(0, n - 1);
		uniform_real_distribution<double> chooseFraction(0.0, 1.0);
		for (int g = 1; g <= nGenerations; ++g) {
			cout << g;
			ofs << g;
			for (int nPop = 0; nPop < nPopulations; ++nPop) {

				//Interactions
				for (int i = 0; i < n; ++i) {
					for (int k = 0; k < nInteractions; ++k) {
						int j = pickPartner(rng);
						if (Populations[nPop][i].strategy == 1 && Populations[nPop][j].strategy == 1) {			//Both cooperate
							Populations[nPop][i].fitness += b - c / 2;
						}
						else if (Populations[nPop][i].strategy == 1 && Populations[nPop][j].strategy == 0) {	//Only focal individual cooperates
							Populations[nPop][i].fitness += b - c;
						}
						else if (Populations[nPop][i].strategy == 0 && Populations[nPop][j].strategy == 1) {	//Only focal individual defects
							Populations[nPop][i].fitness += b;
						}
						else if (Populations[nPop][i].strategy == 0 && Populations[nPop][j].strategy == 0) {	//Both defect
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
						PopulationNew[i].strategy = PopulationNew[i].strategy == 1 ? 0 : 1;
					}
				}
				Populations[nPop] = PopulationNew;

				//Calculate Pc
				double mean = 0.0;
				for (int i = 0; i < n; ++i) {
					if (Populations[nPop][i].strategy == 1) {
						mean += 1.0;
					}
				}

				//Output
				cout << "\t" << mean / n;
				ofs << "\t" << mean / n;
			}
			cout << "\n";
			ofs << "\n";
		}
	}
	catch (exception &error) {
		cerr << "Error: " << error.what();
		return 1;
	}
	
	return 0;
}