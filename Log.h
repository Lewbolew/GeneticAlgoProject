//
// Created by bohdan on 6/19/17.
//

#ifndef MYGA_LOG_H
#define MYGA_LOG_H

#include <fstream>

class Log
{
public:

    Log();
    Log( char* );
    ~  Log();
    void Write( char* txt );
    void Open( const char* );

private:
    std::ofstream m_stream;
};

#endif //MYGA_LOG_H
