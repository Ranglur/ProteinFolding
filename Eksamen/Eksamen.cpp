// Eksamen.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ProteinSystem.h"
#include "threaders.h"
#include <chrono> 
#include <thread>




int main()
{
    std::vector<std::vector<double>> interactionMatrix = readInteractionMatrixFromFile("interaction_matrix.txt");
    if (!isValidInteractionMatrix(interactionMatrix))
        std::cerr << "Invalid interaction matrix";

    //Problem 1.3
    //--------------------------------------------------------------
    //int N = 15;
    //ProteinSystem systemA(N, interactionMatrix, true);
    //systemA.writePositionsToFile("ProteinA.txt");
    //std::cout << systemA.calculateEnergy() << std::endl;

    //ProteinSystem systemB(N, interactionMatrix, true);
    //systemB.writePositionsToFile("ProteinB.txt");
    //std::cout << systemB.calculateEnergy() << std::endl;

    //ProteinSystem systemC(N, interactionMatrix, true);
    //systemC.writePositionsToFile("ProteinC.txt");
    //std::cout << systemC.calculateEnergy() << std::endl;



    // Problem 1.5
    // ---------------------------------------------------------------
    //int N = 15;
    //
    //ProteinSystem system(N, interactionMatrix);
    //system.setTemp(10); 

    //MonteCarloData data = system.nMonteCarloSweeps(1);
    //system.writePositionsToFile("Positions1sweepT10.txt");

    //MonteCarloData data2 = data + system.nMonteCarloSweeps(9);
    //system.writePositionsToFile("Positions10sweepT10.txt");
    //
    //MonteCarloData data3 = data2 + system.nMonteCarloSweeps(100);
    //system.writePositionsToFile("Positions100sweepT10.txt");
    //
    //data3.writeToFile("MonteCarloData.txt");




    // Problem 1.6
    //-----------------------------------------------------------------
    //int N = 15;

    //ProteinSystem system(N, interactionMatrix, false, false);
    //system.setTemp(1);
    //MonteCarloData data = system.nMonteCarloSweeps(1000000);
    //data.writeToFile("MonteCarloDataN15T1.txt");
    //system.writePositionsToFile("Positions10000000sweepT1.txt");
    


    // Problem 1.7
    //----------------------------------------------------------------- 
    // Initialized stretched out


    //std::thread thread1(problem171, interactionMatrix);
    //std::thread thread2(problem172, interactionMatrix);
    //std::thread thread3(problem173, interactionMatrix);

    
    
    // Version where proteins are initialized as not stretched out
    //----------------------------------------------------------------- 

    //std::thread thread4(problem174, interactionMatrix);
    //std::thread thread5(problem175, interactionMatrix);
    //std::thread thread6(problem176, interactionMatrix);


    // problem 1.8
    //----------------------------------------------------------------- 

    //int n = 30;
    //ProteinSystem systemNotAnealed1(30, interactionMatrix);
    //ProteinSystem systemNotAnealed2(30, interactionMatrix);
    //ProteinSystem systemAnealed(30, interactionMatrix);

    //systemNotAnealed1.setTemp(1);
    //systemNotAnealed2.setTemp(1);
    //systemAnealed.setTemp(10);
    //
    //// Sets structure equal
    //std::vector<int> aminoAcids = systemNotAnealed1.getAminoAcids();
    //systemNotAnealed2.setAllAminoAcids(aminoAcids);
    //systemAnealed.setAllAminoAcids(aminoAcids);

    //MonteCarloData notAnealed1 = systemNotAnealed1.nMonteCarloSweeps(10000000);
    //MonteCarloData notAnealed2 = systemNotAnealed2.nMonteCarloSweeps(10000000);
    //
    //// Writes Monte carlo data to file
    //notAnealed1.writeToFile("McDataNotAnealed1.txt");
    //notAnealed1.writeToFile("McDataNotAnealed2.txt");
    //systemAnealed.anealSystem(2, 5000000, false, "McDataAnealed.txt");

    //// Writes positional data to file
    //systemNotAnealed1.writePositionsToFile("PositionsNotAnealed1.txt");
    //systemNotAnealed2.writePositionsToFile("PositionsNotAnealed2.txt");
    //systemAnealed.writePositionsToFile("PositionAnealed.txt");



    // problem 1.9
    //----------------------------------------------------------------- 
    int N = 50;

    std::vector<std::vector<double>> newInteractionMatrix = interactionMatrix;
    // Let theese amino acids positively repell

    newInteractionMatrix[2][4] = - newInteractionMatrix[4][2]; 
    newInteractionMatrix[4][2] = newInteractionMatrix[2][4];

    newInteractionMatrix[16][9] = -newInteractionMatrix[9][16];
    newInteractionMatrix[9][16] = newInteractionMatrix[16][9];

    newInteractionMatrix[5][5] = -newInteractionMatrix[5][5];
    newInteractionMatrix[7][7] = -newInteractionMatrix[7][7];
    newInteractionMatrix[9][9] = -newInteractionMatrix[9][9];

    newInteractionMatrix[11][1] = -newInteractionMatrix[1][11];
    newInteractionMatrix[1][11] = newInteractionMatrix[11][1];

    ProteinSystem systemReppulcion1(N,newInteractionMatrix);
    ProteinSystem systemReppulcion2(N,newInteractionMatrix);

    systemReppulcion1.anealSystem(2, 50000000, false);
    systemReppulcion2.anealSystem(2, 50000000, false);

    systemReppulcion1.writePositionsToFile("AnealedRepulsivePositions1.txt");
    systemReppulcion2.writePositionsToFile("AnealedRepulsivePositions2.txt");


    // problem 2.1
    //----------------------------------------------------------------- 
    //int N = 15;

    //ProteinSystem systemA3D(N, interactionMatrix, true, true);
    //systemA3D.writePositionsToFile("ProteinA3D.txt");
    //std::cout << systemA3D.calculateEnergy() << std::endl;

    //ProteinSystem systemB3D(N, interactionMatrix, true, true);
    //systemB3D.writePositionsToFile("ProteinB3D.txt");
    //std::cout << systemB3D.calculateEnergy() << std::endl;

    //ProteinSystem systemC3D(N, interactionMatrix, true, true);
    //systemC3D.writePositionsToFile("ProteinC3D.txt");
    //std::cout << systemC3D.calculateEnergy() << std::endl;

    


    // problem 2.2
    //----------------------------------------------------------------- 
    //int N = 15;

    //ProteinSystem system3D(N, interactionMatrix,false,true);
    //system3D.setTemp(10); 

    //MonteCarloData data13D = system3D.nMonteCarloSweeps(1);
    //system3D.writePositionsToFile("Positions1sweepT103D.txt");

    //MonteCarloData data23D = data13D + system3D.nMonteCarloSweeps(9);
    //system3D.writePositionsToFile("Positions10sweepT103D.txt");
    //
    //MonteCarloData data33D = data23D + system3D.nMonteCarloSweeps(100);
    //system3D.writePositionsToFile("Positions100sweepT103D.txt");
    //
    //data33D.writeToFile("MonteCarloData3D.txt");


    // problem 2.3
    //----------------------------------------------------------------- 

    //std::thread thread7(problem231, interactionMatrix);
    //std::thread thread8(problem232, interactionMatrix);
    //std::thread thread9(problem233, interactionMatrix);


    
    // Version where proteins are initialized as not stretched out
    //----------------------------------------------------------------- 

    //std::thread thread10(problem234, interactionMatrix);
    //std::thread thread11(problem235, interactionMatrix);
    //std::thread thread12(problem236, interactionMatrix);

    //thread1.join();
    //thread2.join();
    //thread3.join();
    //thread4.join();
    //thread5.join();
    //thread6.join();
    //thread7.join();
    //thread8.join();
    //thread9.join();
    //thread10.join();
    //thread11.join();
    //thread12.join();

    return 0;
}



