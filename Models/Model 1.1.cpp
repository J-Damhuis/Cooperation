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
const double price = 0.001;						//Costs for obtaining information
const int n = 1000;								//Population size
const int nPopulations = 1;						//Number of populations
const int nInteractions = 10;					//Number of interactions in individual's life
const int nGenerations = 500;					//Number of generations
const int nGenSav = 1;							//Save every n generations
const vector<double> Pc = { 0.1 };				//Initial mean tendency of populations to cooperate
const vector<double> Pi = { 0.0 };				//Initial fraction of population which obtains information

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
		ofstream ofs("Model 1.1 nul.csv");
		if (!ofs.is_open()) {
			throw logic_error("Unable to open output file\n");
		}
		ofstream spread1("spread_model1.1_0.95.csv"), spread2("spread_model1.1_0.05.csv"), spread3("spread_model1.1_0.67.csv");
		vector<ofstream> spread;
		spread.push_back(move(spread1));
		spread.push_back(move(spread2));
		spread.push_back(move(spread3));
		if (!spread[0].is_open() || !spread[1].is_open() || !spread[2].is_open()) {
			throw logic_error("Unable to open one of the spread file \n");
		}

		//More checks
		if (nPopulations != Pc.size()) {
			throw logic_error("Number of populations and number of initial mean Pc values are different\n");
		}
		if (nPopulations != Pi.size()) {
			throw logic_error("Number of populations and number of initial mean Pi values are different\n");
		}

		//Initialise population
		vector<vector<Individual> > Populations(nPopulations, vector<Individual>(n));
		cout << "Generation";
		ofs << "Generation";
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			cout << "\tPmean\tMeanP0\tStdevP0\tMeanInfo";
			ofs << "\tPmean\tMeanP0\tStdevP0\tMeanInfo";
		}
		cout << "\n0";
		ofs << "\n0";
		spread[0] << "0";
		spread[1] << "0";
		spread[2] << "0";
		vector<double> pmean = Pc;
		for (int nPop = 0; nPop < nPopulations; ++nPop) {
			double mean = 0.0, info = 0.0;
			int unresponsive = 0;
			for (int i = 0; i < n; ++i) {
				Populations[nPop][i].fitness = 0.0;
				Populations[nPop][i].strategy = Pc[nPop];
				Populations[nPop][i].info = i < n * Pi[nPop] ? 1 : 0;
				if (Populations[nPop][i].info == 0) {
					mean += Populations[nPop][i].strategy;
					++unresponsive;
				}
				info += Populations[nPop][i].info;
				spread[nPop] << "\t" << Populations[nPop][i].strategy;
			}
			mean /= unresponsive;
			info /= n;
			double stdev = 0.0;
			for (int i = 0; i < n; ++i) {
				if (Populations[nPop][i].info == 0) {
					stdev += pow(Populations[nPop][i].strategy - mean, 2);
				}
			}
			stdev = sqrt(stdev / unresponsive);
			cout << "\t" << pmean[nPop] << "\t" << mean << "\t" << stdev << "\t" << info;
			ofs << "\t" << pmean[nPop] << "\t" << mean << "\t" << stdev << "\t" << info;
		}
		cout << "\n";
		ofs << "\n";
		spread[0] << "\n";
		spread[1] << "\n";
		spread[2] << "\n";

		//Simulate
		uniform_real_distribution<double> chooseFraction(0.0, 1.0);
		uniform_int_distribution<int> pickPartner(0, n - 1);
		for (int g = 1; g <= nGenerations; ++g) {
			if (g % nGenSav == 0) {
				cout << g;
				ofs << g;
				spread[0] << g;
				spread[1] << g;
				spread[2] << g;
			}
			for (int nPop = 0; nPop < nPopulations; ++nPop) {

				//Interactions
				double coop = 0.0;
				for (int i = 0; i < n; ++i) {
					for (int k = 0; k < nInteractions; ++k) {
						int j = pickPartner(rng);
						double r = chooseFraction(rng);
						double s = chooseFraction(rng);
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
								if (pmean[nPop] > r && pmean[nPop] > s) {											//Both defect
									Populations[nPop][i].fitness += 0.0;
								}
								else if (pmean[nPop] > r && pmean[nPop] < s) {										//Only focal individual defects
									Populations[nPop][i].fitness += b;
								}
								else if (pmean[nPop] < r && pmean[nPop] > s) {										//Only focal individual cooperates
									Populations[nPop][i].fitness += b - c;
									coop += 1.0;
								}
								else if (pmean[nPop] < r && pmean[nPop] < s) {										//Both cooperate
									Populations[nPop][i].fitness += b - c / 2.0;
									coop += 1.0;
								}
							}
						}
						else if (Populations[nPop][i].info == 0) {						 							//If focal individual does not obtain info
							if (Populations[nPop][j].info == 0) {													//If partner does not obtain info
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
							else if (Populations[nPop][j].info == 1) {												//If partner obtains info
								if (Populations[nPop][i].strategy > r && Populations[nPop][i].strategy > s) {		//Only focal individual cooperates
									Populations[nPop][i].fitness += b - c;
									coop += 1.0;
								}
								else if (Populations[nPop][i].strategy > r && Populations[nPop][i].strategy < s) {	//Both cooperate
									Populations[nPop][i].fitness += b - c / 2.0;
									coop += 1.0;
								}
								else if (Populations[nPop][i].strategy < r && Populations[nPop][i].strategy > s) {	//Both defect
									Populations[nPop][i].fitness += 0.0;
								}
								else if (Populations[nPop][i].strategy < r && Populations[nPop][i].strategy < s) {	//Only focal individual defects
									Populations[nPop][i].fitness += b;
								}
							}
						}
					}
				}
				pmean[nPop] = coop / (nInteractions * n);

				//Determine offspring
				vector<double> vecWeights(n);
				for (int i = 0; i < n; ++i) {
					vecWeights[i] = Populations[nPop][i].fitness < 0.0 ? 0.0 : Populations[nPop][i].fitness;
				}
				vector<Individual> PopulationNew(n);
				//discrete_distribution<int> chooseParent(vecWeights.begin(), vecWeights.end());
				for (int i = 0; i < n; ++i) {
					PopulationNew[i] = Populations[nPop][pickPartner(rng)];
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
					if (chooseFraction(rng) < mu) {
						PopulationNew[i].info = PopulationNew[i].info == 1 ? 0 : 1;
					}
				}
				Populations[nPop] = PopulationNew;

				if (g % nGenSav == 0) {
					//Calculate means
					double mean = 0.0, info = 0.0;
					int unresponsive = 0;
					for (int i = 0; i < n; ++i) {
						if (Populations[nPop][i].info == 0) {
							mean += Populations[nPop][i].strategy;
							++unresponsive;
						}
						info += Populations[nPop][i].info;
					}
					mean /= unresponsive;
					info /= n;

					//Calculate stdev of Pc
					double stdev = 0.0;
					for (int i = 0; i < n; ++i) {
						if (Populations[nPop][i].info == 0) {
							stdev += pow(Populations[nPop][i].strategy - mean, 2);
						}
					}
					stdev = sqrt(stdev / unresponsive);

					//Output
					cout << "\t" << pmean[nPop] << "\t" << mean << "\t" << stdev << "\t" << info;
					ofs << "\t" << pmean[nPop] << "\t" << mean << "\t" << stdev << "\t" << info;
					for (int i = 0; i < n; ++i) {
						spread[nPop] << "\t" << Populations[nPop][i].strategy;
					}
					spread[nPop] << "\n";
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