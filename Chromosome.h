//
// Created by bohdan on 6/19/17.
//

#ifndef MYGA_CHROMOSOME_H
#define MYGA_CHROMOSOME_H
#pragma once

const int chrSize = 64;

class Chromosome
{
public:

    Chromosome(void);
    ~Chromosome(void);

    void SetFitness( const double& value );
    void SetChromosome( const int& index, const unsigned char& value );
    unsigned char GetChromosome( const int& index );
    double GetFitness() const;
    int size() const;
    void Print( const int& index ) const;

private:

    unsigned char chr[ chrSize ];
    double fitness;
};
#endif //MYGA_CHROMOSOME_H
