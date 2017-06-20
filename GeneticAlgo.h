//
// Created by bohdan on 6/19/17.
//

#ifndef MYGA_GENETICALGO_H
#define MYGA_GENETICALGO_H

#include "Population.h"
#include "Log.h"
#include <functional>
using namespace std::placeholders;

class GeneticAlgorithm
{
public:

    GeneticAlgorithm(void);
    ~GeneticAlgorithm(void);

    void Initialize( const int& enc,
                     const int& crate,
                     const int& mrate,
                     const int& psize,
                     const int& iter,
                     const int& csize,
                     const int& tsize,
                     std::function<double(double, double)> fit_func,
                     const std::string& path );
    void Run();

private:

    void CreatePopulation();
    double Evaluate();
    void Crossover();
    void Mutate();
    void Select();

    void SetParameters( const int& enc,
                        const int& crate,
                        const int& mrate,
                        const int& psize,
                        const int& iter,
                        const int& csize,
                        const int& tsize,
                        std::function<double(double,double)> fit_func);

    void LogResult( const double& result,
                    const int& iter,
                    const int& count );

private:

    int encoding;
    int mutationRate;
    int crossoverRate;
    int populationSize;
    int numberIterations;
    int chromosomeSize;
    int tournamentSize;

    int bestFitnessIndex;
    double bestFitness;
    float best_x;
    float best_y;

    Population pop;
    Log log;

};
#endif //MYGA_GENETICALGO_H
