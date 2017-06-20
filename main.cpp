#include "Log.h"
#include <sstream>
#include <cmath>
#include "GeneticAlgo.h"

const int encoding_type = 1;
const int crossover_rate = 70;
const int mutation_rate = 3;
const int population_size = 1000;
const int number_iterations = 1300;
const int chromosome_size = 64;
const int tournament_size = population_size / 10;

double Rosenbrock_function(const double& x, const double& y) {
    // Rosenbrock: (1-x)^2 + 100(y-x*x)^2

    double fitness = ( std::pow( (double)( 1.0 - x ), 2 ) ) +
                     100 * ( pow( (double) ( y - ( x * x ) ), 2 ) ) ;

    return fitness;
}

double Sin_x_y(const double& x, const double& y) {
    double fitness = sin(x * x + y * y);
    return fitness;
}

double x_y_squared(const double& x, const double& y) {
    double fitness = x * x + y * y + 100;
    return fitness;
}

int main()
{
    // Run the GA!
    GeneticAlgorithm ga;

    ga.Initialize( encoding_type,
                   crossover_rate,
                   mutation_rate,
                   population_size,
                   number_iterations,
                   chromosome_size,
                   tournament_size,
                   Rosenbrock_function,
                   "/home/bohdan/Desktop/lool.txt");
    ga.Run();

    return 0;
}