#include <vector>
#include <string>
#include <functional>
#include <mutex>
#include "Chromosome.h"

const double infinity = 9999999999999;

class Population
{
public:
    enum Encoding {
        GRAY,
        IEEE_754,
    };

    Population(void);
    ~Population(void);

    void SetChromosomeEncoding( const int& type );
    void SetChromosomeSize( const int& size );
    void SetFitFunction(std::function<double(double, double)> fit_f);
    void CreateRandomPopulation( const int& size );
    void Crossover( const int& index1,
                    const int& index2,
                    const int& point );

    void Crossover( const int& index1,
                    const int& index2,
                    const int& point1,
                    const int& point2 );

    void Mutation( const int& index );
    double EvaluatePopulation( float& bx,
                               float& by );

    double CalcChromosomeFitness( const int& index,
                                  float& xv,
                                  float& yv);

    double GetChromosomeFitness( const int& index ) const;
    void CopyChromosome( const int& source,
                         const int& dest );

private:
    Chromosome* CreateRandomChromosome();
    std::string GetXstring( Chromosome* chr );
    std::string GetYstring( Chromosome* chr );
    float GetFloat32_IEEE754( std::string Binary );
    float GetFloat32_Gray( std::string Binary );
    int Binary32ToHex( std::string Binary );
    double CalculateFitnessFunction( const float& x,
                                     const float& y );
    double CalcChromosomeFitness_IEEE754( const int& index,
                                          float& xv, float& yv);

    double CalcChromosomeFitnessGray( const int& index,
                                      float& xv, float& yv );
    void additional_function_to_fitness_calc(int i, double& totalFitness, double& bestFitness,
                                         float& bx, float& by, int& bestFitnessIndex,
                                         std::mutex& mut);
    double EvaluatePopulation1( float& bx, float& by );

    std::string Gray2Bin( std::string gray );
    long Bin2Dec( std::string bin );
private:
    std::vector< Chromosome* > pop;
    std::function<double(double, double)> fit_func;
    int chrSize;
    int encoding;
};