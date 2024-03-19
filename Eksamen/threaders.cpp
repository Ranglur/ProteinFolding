#include "threaders.h"


// Creates functions for threading. 
// problem 1.7:
//----------------------------------
void problem171(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(15, interactionMatrix);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN15.txt");
    system.writePositionsToFile("PositionsAnealedN15.txt");

}

void problem172(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(50, interactionMatrix);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN50.txt");
    system.writePositionsToFile("PositionsAnealedN50.txt");

}

void problem173(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(100, interactionMatrix);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN100.txt");
    system.writePositionsToFile("PositionsAnealedN100.txt");

}

void problem174(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(15, interactionMatrix, true, false);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN15RNG.txt");
    system.writePositionsToFile("PositionsAnealedN15RNG.txt");

}

void problem175(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(50, interactionMatrix, true, false);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN50RNG.txt");
    system.writePositionsToFile("PositionsAnealedN50RNG.txt");


}

void problem176(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(100, interactionMatrix, true, false);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN100RNG.txt");
    system.writePositionsToFile("PositionsAnealedN100RNG.txt");

}


// problem 2.3:
//----------------------------------

void problem231(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(15, interactionMatrix, false, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN153D.txt");
    system.writePositionsToFile("PositionsAnealedN153D.txt");

}

void problem232(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(50, interactionMatrix, false, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN503D.txt");
    system.writePositionsToFile("PositionsAnealedN503D.txt");

}

void problem233(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(100, interactionMatrix, false, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN1003D.txt");
    system.writePositionsToFile("PositionsAnealedN1003D.txt");

}

void problem234(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(15, interactionMatrix, true, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN153DRNG.txt");
    system.writePositionsToFile("PositionsAnealedN153DRNG.txt");

}

void problem235(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(50, interactionMatrix, true, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20);
    PCdata.writeToFile("PhaseChangeDataN503DRNG.txt");
    system.writePositionsToFile("PositionsAnealedN503DRNG.txt");

}

void problem236(std::vector<std::vector<double>> interactionMatrix)
{
    ProteinSystem system(100, interactionMatrix, true, true);
    system.setTemp(10);
    PhaseChangeData PCdata = system.anealSystem(3000000, 20, true);
    PCdata.writeToFile("PhaseChangeDataN1003DRNG.txt");
    system.writePositionsToFile("PositionsAnealedN1003DRNG.txt");

}