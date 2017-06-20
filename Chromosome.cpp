//
// Created by bohdan on 6/19/17.
//

#include "Chromosome.h"
#include <iostream>

// Constructor
Chromosome::Chromosome(void) {}

// Destructor
Chromosome::~Chromosome(void) {}

// Set gen of the chromosome
void Chromosome::SetChromosome( const int& index, const unsigned char& value )
{
    if ( index < 0 || index >= chrSize ) return;
    chr[ index ] = value;
}

// Get the gen of the chromosome
unsigned char Chromosome::GetChromosome( const int& index )
{
    unsigned char element = chr[ index ];
    return element;
}

void Chromosome::SetFitness( const double& value )
{
    fitness = value;
}

double Chromosome::GetFitness() const
{
    return fitness;
}

int Chromosome::size() const
{
    return chrSize;
}

// Print the chromosome with it`s fitness
void Chromosome::Print( const int& index ) const
{
    std::string str;
    for ( int i = 0; i < chrSize; i++ ) {
        unsigned char value = chr[ i ];
        str.append( value == 0 ? "0" : "1" );
    }
    std::cout << index << "\t" << str.c_str() << "\t" << fitness << std::endl;
}