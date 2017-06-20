//
// Created by bohdan on 6/19/17.
//

#include "GeneticAlgo.h"
#include <sstream>
#include <stdlib.h>
#include <math.h>


GeneticAlgorithm::GeneticAlgorithm(void)
{
    encoding = 0;
    mutationRate = 5;
    crossoverRate = 90;
    populationSize = 100;
    numberIterations = 10000;
    chromosomeSize = 64;
    tournamentSize = populationSize / 5;
    bestFitness = infinity;
}

GeneticAlgorithm::~GeneticAlgorithm(void)
{}

void GeneticAlgorithm::Initialize( const int& enc,
                                   const int& crate,
                                   const int& mrate,
                                   const int& psize,
                                   const int& iter,
                                   const int& csize,
                                   const int& tsize,
                                   std::function<double(double, double)> fit_func,
                                   const std::string& path )
{
    SetParameters( enc, crate, mrate, psize, iter, csize, tsize, fit_func );
    CreatePopulation();
    log.Open( path.c_str() );
}

void GeneticAlgorithm::Run()
{
    for ( int i = 0; i < numberIterations; i++ )
    {
        LogResult( Evaluate(), i, 10 );
        Select();
        Crossover();
        Mutate();
    }
}

double GeneticAlgorithm::Evaluate()
{
    float bx = -1;
    float by = -1;
    double best = pop.EvaluatePopulation( bx, by );
    if ( best < bestFitness )
    {
        bestFitness = best;
        best_x = bx;
        best_y = by;
    }
    return bestFitness;
}

void GeneticAlgorithm::Crossover()
{
    for ( int i = 0; i < populationSize; i++ )
    {
        int r = rand() % 100;

        if ( r < crossoverRate )
        {
            int index1 = rand() % populationSize;
            int index2 = rand() % populationSize;
            while ( index1 == index2 )
            {
                index2 = rand() % populationSize;
            }

            // Get crossover points
            // Point1: 0 - 31
            int point1 = rand() % chromosomeSize / 2;

            // Point1: 32 - 64
            int point2 = chromosomeSize / 2 +
                         rand() % ( chromosomeSize / 2 );

            while ( point1 == point2 ) {
                point2 = chromosomeSize / 2 +
                         rand() % ( chromosomeSize / 2 );
            }

            if ( point1 > point2 )
            {
                int tmp = point1;
                point1 = point2;
                point2 = tmp;
            }

            // Do 1-point crossover
            pop.Crossover( index1, index2, point1, point2 );
        }
    }
}

void GeneticAlgorithm::Mutate()
{
    for ( int i = 0; i < populationSize; i++ )
    {
        int r = rand() % 100;
        if ( r < mutationRate )
        {
            pop.Mutation( i );
        }
    }
}

// Select population chromosomes according to fitness
// Using a pairwise tournament selection mechanism
void GeneticAlgorithm::Select()
{
    // For each pair of chromosomes selected
    int i = 0;
    while ( i < tournamentSize )
    {
        // Get the chromosome pair for tournament selection
        int index1 = rand() % populationSize;
        int index2 = rand() % populationSize;

        while ( index1 == index2 )
        {
            index2 = rand() % populationSize;
        }

        double fitness1 = fabs( pop.GetChromosomeFitness( index1 ) );
        double fitness2 = fabs( pop.GetChromosomeFitness( index2 ) );

        // We seek to find [x,y] that minimizes this function
        // The bigget the value returned, the lower its fitness
        if ( fitness1 > fitness2 )
        {
            // Copy chromosome 1 elements into chromosome 2
            pop.CopyChromosome( index2, index1 );
        }
        else
        {
            // Copy chromosome 2 elements into chromosome 1
            pop.CopyChromosome( index1, index2 );
        }
        i++;
    }
}

void GeneticAlgorithm::SetParameters( const int& enc,
                                      const int& crate,
                                      const int& mrate,
                                      const int& psize,
                                      const int& iter,
                                      const int& csize,
                                      const int& tsize,
                                      std::function<double(double, double)> fit_func)
{
    mutationRate = mrate;
    crossoverRate = crate;
    populationSize = psize;
    numberIterations = iter;
    chromosomeSize = csize;
    tournamentSize = tsize;
    pop.SetChromosomeEncoding( enc );
    pop.SetFitFunction(fit_func);

}

void GeneticAlgorithm::CreatePopulation()
{
    pop.CreateRandomPopulation( populationSize );
}

// Write results into a file
void GeneticAlgorithm::LogResult( const double& result,
                                  const int& iter,
                                  const int& count )
{
    if ( iter % count == 0 )
    {
        std::stringstream ss;
        //ss << iter << "t" << result << "t" << best_x << "t" << best_y;
        ss <<result<< ", " <<  best_x<<", "<< best_y << "\n";
        log.Write( (char*) ss.str().c_str() );
    }
}