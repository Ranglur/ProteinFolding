#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono> 


//Interaction matrix functions and utility:
//------------------------------------------


std::vector<std::vector<double>> generateNewInteractionMatrix(); // Genreates a new interaction matrix
std::vector<std::vector<double>> readInteractionMatrixFromFile(const std::string& filename); // Reads in an interaction matrix form file
void writeInteractionMatrixToFile(std::vector<std::vector<double>> interactionMatrix, std::string filename); // Writes interaction matrix to file
bool isValidInteractionMatrix(std::vector<std::vector<double>> interactionMatrix); // Validates interaction matrix

double averageVector(std::vector<double> vec, int startingIndex  = 0, int finalIndex = -1);

//Logging data struct for MC
//-----------------------------------------



struct MonteCarloData
{
	int dataSize;
	std::vector<double> endToEndDistance;
	std::vector<double> RoG;
	std::vector<double> energy;
	std::vector<int> monteCarloSweeps;

	MonteCarloData();
	MonteCarloData(int n);
	void writeToFile(const std::string& fileName);
	MonteCarloData& operator+=(const MonteCarloData& rhs); // appending onto the last
};


MonteCarloData operator+(const MonteCarloData& lhs, const MonteCarloData& rhs);


struct PhaseChangeData
{
	int dataSize;
	std::vector<double> avgEndToEndDistance;
	std::vector<double> avgRoG;
	std::vector<double> avgEnergy;
	std::vector<double> temperature;

	PhaseChangeData();
	PhaseChangeData(int n);
	void writeToFile(const std::string& fileName);
};



//ProteinSystem:
//------------------------------------------

class ProteinSystem
{
private:
	int N = 0; // Stores number of monomers in the chain
	double T = 0; // Temperature of the system
	const double k_b = 1;// = 1.380e-23; // Boltzmans constant
	

	std::vector<int> xPositions; // Stores the x-positions of all monomers in chain 
	std::vector<int> yPositions; // Stores the y-positions of all monomers in chain 
	std::vector<int> zPositions; // Stores the z-positions of all monomers in chain (unused and set to 0 if system is in 2D mode)

	const bool is3D; // Boolean determining if system is in 3D mode or not. 
	double energy = 0; // Stores the current system energy
	std::vector<int> aminoAcidType; // Stores the amino acid type for each monomer in chain. 
	gsl_rng* r; //Stores the random number generator for 


	std::vector<std::vector<int>> determineLegalMoves(int index); // Determines legal moves for monomer with index index
	std::vector<int> determineNearestNeigboursIndexI(int i, int posX, int posY, int posZ); // determines nearest neigbours of monomer i. Written with extra pos input, as it is also used by Monte Carlo sweep 
	void generateNearestNeigbours(); // initializes nearest neigbour list
	
	

public:
	std::vector<std::vector<int>> nearestNeigbours; // Stores index of the nearest non-covalent neigbours of each monomer. As to not have to reallocate memory for this, the stored vectors
	// All have size 2, exept for the first and last (as this is the maximum number of noncovalent neigbours), and they are initialized to
	// -1, indicating no nearest neigbour.
	std::vector<std::vector<double>> interactionMatrix; // Stores the interaction matrix for the system
	


	ProteinSystem(); // Default constructor
	ProteinSystem(int n, std::vector<std::vector<double>> interactionMatrix, bool Random = false, bool is3Dinput = false); // Constructor of chain of N monomers. Bool Random can be set tu true if user wants
																								   // chain to not be completely straight when initialized.	
	~ProteinSystem(); // Destructor needed due to GSL pointers
	void setInteractionMatrix(std::vector<std::vector<double>> newInteractionMatrix); // Takes in a new interaction matrix
	void setAmino(int index, int aminoType); // Sets amino acid at given index equal to the given type
	void setTemp(double newT); //Sets temperature to newT
	std::vector<int> getAminoAcids(); //Gets all the amino acids of the system 
	void setAllAminoAcids(std::vector<int> newAcids); // Sets all the amino acids of the system


	
	void monteCarloDraw(); // Preforms a Monte Carlo draw
	MonteCarloData nMonteCarloSweeps(int n); // Preforms n Monte Carlo sweeps
	PhaseChangeData anealSystem(int itersBetweenTempchange, int tempsteps, bool showProgress = false, std::string filename = ""); // Aneals the system by slowly decreasing the temperature from 10 towards 0.  
	                                                                                                   // showProgress show how far along in the process we are. By ajusting the to ints this function solves problem 1.7 and 1.8. 
																									   // if filename is nonempty anealsystem wil write its MonteCarloData to a file 
	double calculateEnergy(); // Calculates the energy of the system 
	double calculateRoG();   // Calculates RoG assuming all monomers roughly equal mass
	double calculateEndToEndDistance(); // Calculates distance from start to end of the protein
	
	std::vector<std::vector<int>> generatePossiblePositions(int i);
	void writePositionsToFile(const std::string& filename); // Writes contents of positional data to file, such that it can be plotted in python
	void writeAminoAcidTypeToFile(const std::string& filename); // Writes positions to file 
	bool hasDiagonal();

};

std::vector<std::vector<int>> initializePosiblePositions(bool is3D);