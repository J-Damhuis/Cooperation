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
const double price = 0.1;						//Costs for obtaining information
const int n = 1000;								//Population size
const int nPopulations = 3;						//Number of populations
const int nInteractions = 10;					//Number of interactions in individual's life
const int nGenerations = 100;					//Number of generations
const int nGenSav = 1;							//Save every n generations
const vector<double> Pc = { 0.95, 0.05, 0.67 };	//Initial mean tendency of populations to cooperate
const vector<double> Pi = { 0.5, 0.5, 0.5 };	//Initial fraction of population which obtains information

mt19937_64 rng;

struct Individual
{
	double strategy;		//0 always uses D; 1 always uses C
	double info;			//0 never obtains information; 1 always does so
	double fitness;			//Cumulative payoff values
};

int main()
{
	try {
		//Seed rng
		rng.seed(1);

		//Open output file
		ofstream ofs("Model 1.0.csv");
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
		double pmean;
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			double mean = 0.0, info = 0.0;
			for (int i = 0; i < n; ++i) {
				Populations[nPop][i].fitness = 0.0;
				Populations[nPop][i].strategy = Pc[nPop];
				Populations[nPop][i].info = i < n * Pi[nPop] ? 1 : 0;
				mean += Populations[nPop][i].strategy;
				info += Populations[nPop][i].info;
			}
			mean /= n;
			pmean = mean;
			info /= n;
			double stdev = 0.0;
			for (int i = 0; i < n; ++i) {
				stdev += pow(Populations[nPop][i].strategy - mean, 2);
			}
			stdev = sqrt(stdev / n);
			cout << "\t" << pmean << "\t" << mean << "\t" << stdev << "\t" << info;
			ofs << "\t" << pmean << "\t" << mean << "\t" << stdev << "\t" << info;
		}
		cout << "\n";
		ofs << "\n";

		//Simulate
		for (int g = 1; g <= nGenerations; ++g) {
			double coop = 0.0;
			if (g % nGenSav == 0) {
				cout << g;
				ofs << g;
			}
			for (int nPop = 0; nPop < nPopulations; ++nPop) {

				//Interactions
				for (int i = 0; i < n; ++i) {
					for (int k = 0; k < nInteractions; ++k) {
						uniform_real_distribution<double> choosePc(0.0, 1.0);
						uniform_int_distribution<int> pickPartner(0, n - 1);
						int j = pickPartner(rng);
						double r = choosePc(rng);
						double s = choosePc(rng);
						if (Populations[nPop][i].info == 1) {														//If focal individual obtains info
							Populations[nPop][i].fitness -= price;
							if (Populations[nPop][j].info == 0) {													//If partner does not obtain info
								if (Populations[nPop][j].strategy > r && Populations[nPop][j].strategy > s) {		//Only focal individual defects
									Populations[nPop][i].fitness += b;
								}
								else if (Populations[nPop][j].strategy > r && Populations[nPop][j].strategy < s) {	//Both defect
									Populations[nPop][i].fitness += 0.0;
								}
								else if (Populations[nPop][j].strategy < r && Populations[nPop][j].strategy > s) {	//Both cooperate
									coop += 1.0;
									Populations[nPop][i].fitness += b - c / 2.0;
								}
								else if (Populations[nPop][j].strategy < r && Populations[nPop][j].strategy < s) {	//Only focal individual cooperates
									Populations[nPop][i].fitness += b - c;
									coop += 1.0;
								}
							}
							else if (Populations[nPop][j].info == 1) {												//If partner does obtain info
								if (pmean > r && Populations[nPop][j].strategy > s) {								//Both cooperate
									Populations[nPop][i].fitness += b - c / 2.0;
									coop += 1.0;
								}
								else if (pmean > r && Populations[nPop][j].strategy < s) {							//Only focal individual cooperates
									Populations[nPop][i].fitness += b - c;
									coop += 1.0;
								}
								else if (pmean < r && Populations[nPop][j].strategy > s) {							//Only focal individual defects
									Populations[nPop][i].fitness += b;
								}
								else if (pmean < r && Populations[nPop][j].strategy < s) {							//Both defect
									Populations[nPop][i].fitness += 0.0;
								}
							}
						}
						else {																					//If focal individual does not obtain info
							if (Populations[nPop][i].strategy > r && Populations[nPop][j].strategy > s) {		//Both cooperate
								Populations[nPop][i].fitness += b - c / 2.0;
								coop += 1.0;
							}
							else if (Populations[nPop][i].strategy > r && Populations[nPop][j].strategy < s) {	//Only focal individual cooperates
								Populations[nPop][i].fitness += b - c;
								coop += 1.0;
							}
							else if (Populations[nPop][i].strategy < r && Populations[nPop][j].strategy > s) {	//Only focal individual defects
								Populations[nPop][i].fitness += b;
							}
							else if (Populations[nPop][i].strategy < r && Populations[nPop][j].strategy < s) {	//Both defect
								Populations[nPop][i].fitness += 0.0;
							}
						}
					}
				}
				pmean = coop / (nInteractions * n);

				//Determine offspring
				vector<double> vecWeights(n);
				for (int i = 0; i < n; ++i) {
					vecWeights[i] = Populations[nPop][i].fitness;
				}
				vector<Individual> PopulationNew(n);
				discrete_distribution<int> chooseParent(vecWeights.begin(), vecWeights.end());
				vector<double> vecMutation = { mu, 1 - mu };
				discrete_distribution<int> chooseMutation(vecMutation.begin(), vecMutation.end());
				for (int i = 0; i < n; ++i) {
					PopulationNew[i] = Populations[nPop][chooseParent(rng)];
					PopulationNew[i].fitness = 0.0;
					if (chooseMutation(rng) == 0) {
						normal_distribution<double> defineMutation(0.0, sigma);
						PopulationNew[i].strategy += defineMutation(rng);
						if (PopulationNew[i].strategy < 0.0) {
							PopulationNew[i].strategy = 0.0;
						}
						else if (PopulationNew[i].strategy > 1.0) {
							PopulationNew[i].strategy = 1.0;
						}
					}
					if (chooseMutation(rng) == 0) {
						PopulationNew[i].info = PopulationNew[i].info == 1 ? 0 : 1;
					}
				}
				Populations[nPop] = PopulationNew;

				if (g % nGenSav == 0) {
					//Calculate means
					double mean = 0.0, info = 0.0;
					for (int i = 0; i < n; ++i) {
						mean += Populations[nPop][i].strategy;
						info += Populations[nPop][i].info;
					}
					mean /= n;
					info /= n;

					//Calculate stdev of Pc
					double stdev = 0.0;
					for (int i = 0; i < n; ++i) {
						stdev += pow(Populations[nPop][i].strategy - mean, 2);
					}
					stdev = sqrt(stdev / n);

					//Output
					cout << "\t" << pmean << "\t" << mean << "\t" << stdev << "\t" << info;
					ofs << "\t" << pmean << "\t" << mean << "\t" << stdev << "\t" << info;
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