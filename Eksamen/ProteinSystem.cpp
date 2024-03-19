
#include "ProteinSystem.h"


//Interaction matrix functions, and utility:
//------------------------------------------


// Intializes 4 or 6 x 3 array of zeros.
std::vector<std::vector<int>> initializePosiblePositions(bool is3D) {
	const std::vector<int> zeros(3, 0); // Utility vector
	if (is3D) {
		return std::vector<std::vector<int>>(6, zeros);
	}
	else {
		return std::vector<std::vector<int>>(4, zeros);
	}
}

// Generates a new 20x20 symmetric interaction matrix with random values between -2 and -4.
std::vector<std::vector<double>> generateNewInteractionMatrix()
{
	const int matrix_size = 20;
	double min_value = -4.0;
	double max_value = -2.0;

	// Initialize the random number generator
	gsl_rng_env_setup();
	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);
	gsl_rng_set(r, time(NULL));

	std::vector<std::vector<double>> interaction_matrix(matrix_size, std::vector<double>(matrix_size));

	// generates a symetric matrix M with elements picked from uniform dist between -2 and -4
	for (int i = 0; i < matrix_size; ++i) {
		for (int j = i; j < matrix_size; ++j) {
			double random_value = min_value + gsl_rng_uniform(r) * (max_value - min_value);
			interaction_matrix[i][j] = random_value;
			interaction_matrix[j][i] = random_value;
		}
	}

	// Freeing the RNG memmory
	gsl_rng_free(r);

	return interaction_matrix;
}


// Reads interaction matrix from file
std::vector<std::vector<double>> readInteractionMatrixFromFile(const std::string& filename)
{
	// Intitializes new matrix
	std::vector<std::vector<double>> interaction_matrix;
	std::ifstream file(filename);

	// Reads elements from file, assuming same format as write to file function is used
	if (file.is_open()) 
	{
		std::string line;
		
		while (std::getline(file, line)) 
		{
			std::vector<double> row;
			std::istringstream line_stream(line);
			double value;

			while (line_stream >> value) 
			{
				row.push_back(value);
			}

			interaction_matrix.push_back(row);
		}
		file.close();
	}
	else 
	{
		std::cerr << "Unable to open file for reading: " << filename << std::endl;
	}

	return interaction_matrix;
}


// Writes interaction matrix to file
void writeInteractionMatrixToFile(std::vector<std::vector<double>> interactionMatrix, std::string filename)
{
	std::ofstream file(filename);

	if (file.is_open()) {
		for (const auto& row : interactionMatrix) {
			for (size_t i = 0; i < row.size(); ++i) {
				file << row[i];
				if (i < row.size() - 1) {
					file << "\t";
				}
			}
			file << "\n";
		}
		file.close();
	}
	else {
		std::cerr << "Unable to open file for writing: " << filename << std::endl;
	}
}

// Verifies interaction matrix
bool isValidInteractionMatrix(std::vector<std::vector<double>> interactionMatrix)
{
	// Test if axis 1 is correct
	if (interactionMatrix.size() != 20)
		return false;

	// Test if axis 2 is correct
	for (int i = 0; i < interactionMatrix.size(); i++)
	{
		if (interactionMatrix[i].size() != 20)
			return false;
	}

	// Test if it is symmetric
	for (int i = 0; i < interactionMatrix.size(); i++)
	{
		for (int j = 0; i < interactionMatrix.size(); i++)
		{
			if (interactionMatrix[i][j] != interactionMatrix[j][i])
				return false;
		}
	}

	
	return true;
}


// Preformes averaging of vector
double averageVector(std::vector<double> vec, int startingIndex, int finalIndex)
{
	if (finalIndex == -1)
		finalIndex = vec.size() - 1;

	double result = 0;
	int elems = 0;

	// Sums over all indecies up to and including the final
	for (int i = startingIndex; i <= finalIndex; i++)
	{
		result += vec[i];
		elems++;
	}
	return result/static_cast<double>(elems);
}


//Logging data struct for MC
//-----------------------------------------

MonteCarloData::MonteCarloData(int n)
{
	dataSize = n;
	endToEndDistance.resize(n);
	RoG.resize(n);
	energy.resize(n);
	monteCarloSweeps.resize(n);
}

MonteCarloData::MonteCarloData()
{
	dataSize = 0;
}

// Writes contents of MonteCarloData to a file. 
void MonteCarloData::writeToFile(const std::string& fileName)
{
	std::ofstream outFile(fileName);

	if (outFile.is_open())
	{
		outFile << "Sweep EndToEnd RoG Energy\n";
		for (int i = 0; i < dataSize; ++i)
		{
			outFile << monteCarloSweeps[i] << " "
				<< endToEndDistance[i] << " "
				<< RoG[i] << " "
				<< energy[i] << "\n";

		}
		outFile.close();
	}
	else
	{
		std::cerr << "Unable to open file " << fileName << std::endl;
	}
}


// Overloaded + operator  makes new MonteCarloData struct which appends the data of the first onto the next.
MonteCarloData operator+(const MonteCarloData& lhs, const MonteCarloData& rhs) {
	MonteCarloData result = lhs;

	// Append the vectors of the second object (rhs) onto the first object (result)
	result.endToEndDistance.insert(result.endToEndDistance.end(), rhs.endToEndDistance.begin(), rhs.endToEndDistance.end());
	result.RoG.insert(result.RoG.end(), rhs.RoG.begin(), rhs.RoG.end());
	result.energy.insert(result.energy.end(), rhs.energy.begin(), rhs.energy.end());

	// Add an offset to the monteCarloSweeps of the second object (rhs) before appending them
	int lastSweep = lhs.monteCarloSweeps.back();
	std::vector<int> updatedSweeps(rhs.monteCarloSweeps.size());
	std::transform(rhs.monteCarloSweeps.begin(), rhs.monteCarloSweeps.end(), updatedSweeps.begin(), [lastSweep](int sweep) { return sweep + lastSweep; });

	result.monteCarloSweeps.insert(result.monteCarloSweeps.end(), updatedSweeps.begin(), updatedSweeps.end());

	// Update dataSize
	result.dataSize = result.endToEndDistance.size();

	return result;
}

// Overloaded += operator
MonteCarloData& MonteCarloData::operator+=(const MonteCarloData& rhs)
{
	// Append the vectors of the second object (rhs) onto the first object (*this)
	this->endToEndDistance.insert(this->endToEndDistance.end(), rhs.endToEndDistance.begin(), rhs.endToEndDistance.end());
	this->RoG.insert(this->RoG.end(), rhs.RoG.begin(), rhs.RoG.end());
	this->energy.insert(this->energy.end(), rhs.energy.begin(), rhs.energy.end());

	// Check if monteCarloSweeps is empty and set lastSweep accordingly
	int lastSweep = this->monteCarloSweeps.empty() ? 0 : this->monteCarloSweeps.back();

	// Add an offset to the monteCarloSweeps of the second object (rhs) before appending them
	std::vector<int> updatedSweeps(rhs.monteCarloSweeps.size());
	std::transform(rhs.monteCarloSweeps.begin(), rhs.monteCarloSweeps.end(), updatedSweeps.begin(), [lastSweep](int sweep) { return sweep + lastSweep; });

	this->monteCarloSweeps.insert(this->monteCarloSweeps.end(), updatedSweeps.begin(), updatedSweeps.end());

	// Update dataSize
	this->dataSize = this->endToEndDistance.size();

	return *this;
}


//Logging data struct for Anealing
//-----------------------------------------
PhaseChangeData::PhaseChangeData(int n)
{
	dataSize = n;
	avgEndToEndDistance.resize(n);
	avgRoG.resize(n);
	avgEnergy.resize(n);
	temperature.resize(n);
}

PhaseChangeData::PhaseChangeData()
{
	dataSize = 0;
}


// Writes contents of MonteCarloData to a file. 
void PhaseChangeData::writeToFile(const std::string& fileName)
{
	std::ofstream outFile(fileName);

	if (outFile.is_open())
	{
		outFile << "Temp avgEndToEnd avgRoG avgEnergy\n";
		for (int i = 0; i < dataSize; ++i)
		{
			outFile << temperature[i] << " "
				<< avgEndToEndDistance[i] << " "
				<< avgRoG[i] << " "
				<< avgEnergy[i] << "\n";

		}
		outFile.close();
	}
	else
	{
		std::cerr << "Unable to open file " << fileName << std::endl;
	}
}



//ProteinSystem:
//------------------------------------------

// Constructors and initializers
ProteinSystem::ProteinSystem() : is3D(false)
{
	gsl_rng_env_setup();
	const gsl_rng_type* T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	// Get a unique seed value for each run
	auto now = std::chrono::high_resolution_clock::now().time_since_epoch();
	unsigned long seed = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
	gsl_rng_set(r, seed);
	energy = 0;

}


ProteinSystem::ProteinSystem(int n, std::vector<std::vector<double>> iM, bool Random, bool is3Dinput) : is3D(is3Dinput), N(n), interactionMatrix(iM)
{
	xPositions = std::vector<int>(N, 0);
	yPositions = std::vector<int>(N, 0);
	zPositions = std::vector<int>(N, 0);  // Will remain a zero vector if is3D is false, as it is by default. 
	aminoAcidType = std::vector<int>(N, 0);

	// Sets up RNG
	gsl_rng_env_setup();
	const gsl_rng_type* T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	// Get a unique seed value for each run
	auto now = std::chrono::high_resolution_clock::now().time_since_epoch();
	unsigned long seed = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
	gsl_rng_set(r, seed);


	// Initializes positions and amino acid type
	for (int i = 0; i < N; i++)
	{
		xPositions[i] = i;
		aminoAcidType[i] = gsl_rng_uniform_int(r, 20);
	}

	// Makes a random folded protein if random is true
	if (Random)
	{
		// Sets up maximum nuber of atempts such that process does not run infinitly
		bool generationComplete = false;
		int maxAtempts = 50;

		while (!generationComplete && maxAtempts > 0)
		{
			generationComplete = true;

			for (int i = 1; i < N; i++)
			{
				// generate posible positions
				std::vector<std::vector<int>> possiblePositions = generatePossiblePositions(i);

				if (possiblePositions.empty())
				{
					generationComplete = false;
					break;
				}


				// Draws random number
				int index = gsl_rng_uniform_int(r, possiblePositions.size());

				xPositions[i] = possiblePositions[index][0];
				yPositions[i] = possiblePositions[index][1];
				zPositions[i] = possiblePositions[index][2];
			}

			--maxAtempts;
		}

		// Gives error message if max atempts was reached
		if (!generationComplete)
		{
			std::cerr << "Could not generate protein (50 attempts)" << std::endl;
		}
	}


	
	generateNearestNeigbours();
	
	double eg = calculateEnergy();
	
}


ProteinSystem::~ProteinSystem()
{
	gsl_rng_free(r);
}


std::vector<int> ProteinSystem::determineNearestNeigboursIndexI(int i, int posX, int posY, int posZ)
{
	

	std::vector<int> result(3, -1);
	if (is3D)
		result.resize(5, -1);

	// Itterates over all elements to find nearest neigbours of i
	for (int j = 0; j < N; j++)
	{
		// continues the for loop if i = j or i and j have covalent bonds (i and j are off by 1)
		if (abs(i - j) <= 1)
			continue;

		// Test if postiotion l-1 norm is 1 for pair j and i. For system in 2d mode the last term will allways be 0
		if ((abs(posX - xPositions[j]) + abs(posY - yPositions[j])) + abs(posZ - zPositions[j]) == 1)
		{
			// Itterates over the index list to find a free spot. When a free spot is found, it breaks. 
			for (int k = 0; k < result.size(); k++)
			{
				if (result[k] == -1)
				{
					result[k] = j;
					break;
				}
			}
		}
	}

	return result;
}

void ProteinSystem::generateNearestNeigbours()
{
	// Intitializes nearest neigbour to be Nx3/5 matrix where every element is -1. Negative 1 elements represent the particle not having a nearest neigbour.   
	
	// Adjusts nearest neigbours if system is in 3d mode. 
	int a = 3;
	if (is3D)
		a = 5;

	std::vector<int> emptyVec(a, -1);



	std::vector<std::vector<int>> nnEmpty(N, emptyVec);
	nearestNeigbours = nnEmpty;


	// Fill the nearest neigbour list:
	for (int i = 0; i < N; i++)
	{
		
		nearestNeigbours[i] = determineNearestNeigboursIndexI(i, xPositions[i], yPositions[i], zPositions[i]);
		
	}
}


// Getters and setters
void ProteinSystem::setInteractionMatrix(std::vector<std::vector<double>> newInteractionMatrix)
{
	interactionMatrix = newInteractionMatrix;
}


void ProteinSystem::setAmino(int index, int aminoType) // Sets amino acid at given index equal to the given type
{
	if (index < N && index >= 0)
		aminoAcidType[index] = aminoType;
}

std::vector<int> ProteinSystem::getAminoAcids() //Gets all the amino acids of the system 
{
	return aminoAcidType;
}


void ProteinSystem::setAllAminoAcids(std::vector<int> newAcids) // Sets all the amino acids of the system
{
	aminoAcidType = newAcids;
}



void ProteinSystem::setTemp(double newT)
{
	T = newT;
}




// Operators on the system


std::vector<std::vector<int>> ProteinSystem::determineLegalMoves(int index)
{
	int coIndex = N + 1;
	std::vector<std::vector<int>> legalMoves;


	if (index > 1 && index < N - 1)
	{
		// this together is a R^3 vector pointing to the oposite corner of the "vector" with index index. 
		int posibleXpos = xPositions[index - 1] + xPositions[index + 1] - xPositions[index];
		int posibleYpos = yPositions[index - 1] + yPositions[index + 1] - yPositions[index];
		int posibleZpos = zPositions[index - 1] + zPositions[index + 1] - zPositions[index]; // this will allways be 0 if we are in 2D mode. Will slow down computation by a hair. 

		// Test to se if this posible state is occupied
		for (int i = 0; i < N; i++)
		{
			// returns the empty legalMoves if it is occupied
			if (xPositions[i] == posibleXpos && yPositions[i] == posibleYpos && zPositions[i] == posibleZpos)
			{
				return legalMoves;
			}
		}


		// Adds the positions to the legal list if they have passed the tests
		legalMoves.push_back({ posibleXpos, posibleYpos, posibleZpos });

	}

	// Assigns value to coIndex if index is at one of the end points
	if (index == 0)
		coIndex = 1;

	if (index == N - 1)
		coIndex = N - 2;


	// Unit directions
	const std::vector<std::vector<int>> possibleMoves = {
	{1, 0, 0},
	{-1, 0, 0},
	{0, 1, 0},
	{0, -1, 0},
	{0, 0, 1},
	{0, 0, -1}
	};

	if (coIndex < N)
	{
		std::vector<std::vector<int>> posiblePositions;
		for (const auto& move : possibleMoves) {

			int newX = xPositions[coIndex] + move[0];
			int newY = yPositions[coIndex] + move[1];
			int newZ = zPositions[coIndex] + move[2];

			// If we are in 2D mode, skip moves in the Z direction
			if (!is3D && (move[2] != 0)) {
				continue;
			}

			// Check if the position is occupied
			bool isOccupied = false;
			for (int j = 0; j < N; j++) {
				if (xPositions[j] == newX && yPositions[j] == newY && zPositions[j] == newZ) {
					isOccupied = true;
					break;
				}
			}

			// Add the position to the list of legal moves if it is not occupied
			if (!isOccupied) {
				posiblePositions.push_back({ newX, newY, newZ });
			}
		}
		legalMoves.insert(legalMoves.end(), posiblePositions.begin(), posiblePositions.end());
	}
	return legalMoves;
}



void ProteinSystem::monteCarloDraw()
{
	
	
	int index = gsl_rng_uniform_int(r, N);
	std::vector<std::vector<int>> legalMoves = determineLegalMoves(index);
	
	
	if (legalMoves.empty())
	{
		
		return;
	}
	


	// Choose a random legal new position
	int randomIndex = gsl_rng_uniform_int(r, legalMoves.size());
	
	


	int newX = legalMoves[randomIndex][0];
	int newY = legalMoves[randomIndex][1];
	int newZ = legalMoves[randomIndex][2];


	// Finds difference in energy
	std::vector<int> newNeigbours = determineNearestNeigboursIndexI(index, newX, newY, newZ);
	std::vector<int> oldNeigbours = nearestNeigbours[index];
	double dE = 0;  // Difference in energy between new move and old protein


	
	for (int i = 0; i < oldNeigbours.size(); i++)
	{
		// Subtracts interaction energy of old neigbours
		if (oldNeigbours[i] != -1)
			dE -= interactionMatrix[aminoAcidType[index]][aminoAcidType[i]];
		
		// Adds interaction energy of new neigbours
		if (newNeigbours[i] != -1)
			dE += interactionMatrix[aminoAcidType[index]][aminoAcidType[i]];
	}



	// Metropolis algorithm
	if (dE <= 0 || gsl_rng_uniform(r) < exp(-dE / (k_b * T)))
	{
	
		
		// Accept move: Update the position
		xPositions[index] = newX;
		yPositions[index] = newY;
		zPositions[index] = newZ; // In 2D mode this should allways be set to zero by determineNearestNeigboursIndexI

		// Updates energy
		energy += dE;

		// Updates nearest neigbours
		nearestNeigbours[index] = newNeigbours;

		// Update the old and new neigbours nearest neigbours
		int k = 0; // Helper int to simplify expressions


		for (int i = 0; i < nearestNeigbours[0].size(); i++)
		{
	
			// Updates nearest neigbour list of new neigbour with index k
			k = nearestNeigbours[index][i];
			if (k != -1)
				nearestNeigbours[k] = determineNearestNeigboursIndexI(k, xPositions[k], yPositions[k], zPositions[k]);

			// Updates nearest neigbour list of old neigbour with index k
			k = oldNeigbours[i];
			if (k != -1)
				nearestNeigbours[k] = determineNearestNeigboursIndexI(k, xPositions[k], yPositions[k], zPositions[k]);
		
		}


	}
}

MonteCarloData ProteinSystem::nMonteCarloSweeps(int n)
{		
	MonteCarloData data(n);

	for (int sweep = 0; sweep < n; ++sweep)
	{
		// Perform N Monte Carlo draws (one for each node in the protein chain)
		for (int i = 0; i < N; ++i)
		{
			monteCarloDraw();
		}

		// Calculate system properties and store them in MonteCarloData struct
		data.endToEndDistance[sweep] = calculateEndToEndDistance();
		data.RoG[sweep] = calculateRoG();
		data.energy[sweep] = calculateEnergy();
		data.monteCarloSweeps[sweep] = sweep;
	}

	return data;
	
}

PhaseChangeData ProteinSystem::anealSystem(int itersBetweenTempchange, int tempsteps, bool ShowProgress, std::string filename)
{
	PhaseChangeData PCdata(tempsteps);
	setTemp(10);
	// Sets dT to ble slightly smaller than T/timesteps as to not get negative temperatures
	double dT = T / static_cast<double>(tempsteps + 1);
	std::vector<double> temperatures(tempsteps, T);


	MonteCarloData allMCData(0);

	// generates list of temperatures
	for (int i = 0; i < tempsteps; i++)
		temperatures[i] -= dT * static_cast<double>(i);

	for (int i = 0; i < tempsteps; i++)
	{
		setTemp(temperatures[i]);
		
		if (ShowProgress)
			std::cout << "Temperature step: " << i + 1 << " / " << tempsteps << " (T = " << temperatures[i] << ")\n";
		
		MonteCarloData MCdata = nMonteCarloSweeps(itersBetweenTempchange);

		// Collect MCData
		if (filename != "")
			allMCData += MCdata;

		
		// Averages the quantities over the last half of the elements
		PCdata.temperature[i] = temperatures[i];
		PCdata.avgEnergy[i] = averageVector(MCdata.energy,itersBetweenTempchange/2);
		PCdata.avgRoG[i] = averageVector(MCdata.RoG, itersBetweenTempchange / 2);
		PCdata.avgEndToEndDistance[i] = averageVector(MCdata.endToEndDistance, itersBetweenTempchange / 2);
	}
	
	// Writes MCData to file
	if (filename != "")
		allMCData.writeToFile(filename);

	return PCdata;

}




double ProteinSystem::calculateEnergy()
{
	
	double result = 0;
	// Run over all nearest neigbours, and use it to calculate the energy of the system
	for (int i = 0; i < N; i++)
	{
		
		for (int j = 0; j < nearestNeigbours[i].size(); j++)
		{
			
			// Gets index of nearest neigbour
			int k = nearestNeigbours[i][j];
			
			// Test if index is dummy index
			if (k == -1)
				continue;

			
			result += interactionMatrix[aminoAcidType[i]][aminoAcidType[k]];
		}
	}

	// Returns result divided by 2, as each nearest neigbour contribution is counted twice. Allso updates the saved system energy
	result /= 2;
	energy = result;
	return result;

}

double ProteinSystem::calculateRoG()
{

	// Calculate the center of mass (COM)
	double comX = 0;
	double comY = 0;
	double comZ = 0;
	for (int i = 0; i < N; ++i) {
		comX += static_cast<double>(xPositions[i]);
		comY += static_cast<double>(yPositions[i]);
		comZ += static_cast<double>(zPositions[i]); // Allways in 2D mode
	}
	comX /= N;
	comY /= N;
	comZ /= N;

	// Calculate the mean squared distance of each lattice point from the COM
	double meanSquaredDistance = 0;
	for (int i = 0; i < N; ++i) {
		double deltaX = xPositions[i] - comX;
		double deltaY = yPositions[i] - comY;
		double deltaZ = zPositions[i] - comZ;
		meanSquaredDistance += deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
	}
	meanSquaredDistance /= N;

	// Calculate the radius of gyration (RoG)
	double radiusOfGyration = std::sqrt(meanSquaredDistance);
	return radiusOfGyration;
}

double ProteinSystem::calculateEndToEndDistance()
{
	// Note about 3D: since in 3d mode zPositions[i] = 0 for all i the code works for 2d as well
	
	// Converts start and end positions to doubles
	double posXstart = static_cast<double>(xPositions[0]);
	double posYstart = static_cast<double>(yPositions[0]);
	double posZstart = static_cast<double>(zPositions[0]);

	
	double posXend = static_cast<double>(xPositions.back());
	double posYend = static_cast<double>(yPositions.back());
	double posZend = static_cast<double>(zPositions.back());

	// Returns absolute value of the vector betwen the points
	return sqrt(pow(posXstart - posXend, 2) + pow(posYstart - posYend, 2) + pow(posZstart - posZend, 2));

}

std::vector<std::vector<int>> ProteinSystem::generatePossiblePositions(int i)
{
	
	// Initializing set of vectors
	
	
	std::vector<std::vector<int>> possiblePositions = initializePosiblePositions(is3D);
	
	

	// Finds the possible positions
	possiblePositions[0] = std::vector<int>{ xPositions[i - 1] + 1, yPositions[i - 1], zPositions[i - 1]};
	possiblePositions[1] = std::vector<int>{ xPositions[i - 1] - 1, yPositions[i - 1], zPositions[i - 1]};
	possiblePositions[2] = std::vector<int>{ xPositions[i - 1], yPositions[i - 1] + 1, zPositions[i - 1]};
	possiblePositions[3] = std::vector<int>{ xPositions[i - 1], yPositions[i - 1] - 1, zPositions[i - 1]};

	if (is3D)
	{	
		possiblePositions[4] = std::vector<int>{ xPositions[i - 1], yPositions[i - 1], zPositions[i - 1] + 1};
		possiblePositions[5] = std::vector<int>{ xPositions[i - 1], yPositions[i - 1], zPositions[i - 1] - 1};

	}

	

	// Erases any disalowed positions
	// Loops over all positions that have allready been set, i.e from 0 to i
	for (int j = 0; j < i; j++)
	{
		// Loops over possible positions
		for (int l = 0; l < possiblePositions.size() ; l++)
		{
			// Checks if any positions are occupied
			if (xPositions[j] == possiblePositions[l][0] && yPositions[j] == possiblePositions[l][1] && zPositions[j] == possiblePositions[l][2])
			{
				
				possiblePositions.erase(possiblePositions.begin() + l);
				
				// reduces the index by 1, as the lth element was removed. 
				l--;
				
			}
		}		
	
	}

	

	return possiblePositions;
}

void ProteinSystem::writePositionsToFile(const std::string& filename)
{
	std::ofstream file(filename);

	if (file.is_open()) {
		for (int i = 0; i < N; ++i) {
			file << xPositions[i] << "\t" << yPositions[i];
			
			if (is3D)
				file << "\t" << zPositions[i];
			
			if (i < N - 1) 
				file << "\n";
			
		}
		file.close();
	}
	else {
		std::cerr << "Unable to open file for writing: " << filename << std::endl;
	}
}

void ProteinSystem::writeAminoAcidTypeToFile(const std::string& filename) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Failed to open the file: " << filename << std::endl;
		return;
	}

	for (int i = 0; i < aminoAcidType.size(); ++i) {
		outFile << aminoAcidType[i];
		if (i < aminoAcidType.size() - 1) {
			outFile << ",";
		}
	}

	outFile.close();
}


bool ProteinSystem::hasDiagonal() {
	for (int i = 0; i < N - 1; ++i) {
		int dx = std::abs(xPositions[i] - xPositions[i + 1]);
		int dy = std::abs(yPositions[i] - yPositions[i + 1]);
		int dz = std::abs(zPositions[i] - zPositions[i + 1]);

		// Check for diagonal movement in 2D and 3D cases
		if (dx + dy + dz > 1) {
			return true;
		}
	}
	return false;
}