//
// Created by bohdan on 6/19/17.
//

#include "Log.h"

Log::Log()
{}

Log::Log( char* filename )
{
    m_stream.open( filename, std::fstream::app );
}

void Log::Write( char* txt )
{
    m_stream << txt << std::endl;
}

Log::~Log()
{
    m_stream.close();
}

// Open the log file for writing
void Log::Open( const char* filename )
{
    m_stream.open( filename, std::fstream::app );
}