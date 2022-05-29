//
// Created by Jan on 25.05.2022.
//

#ifndef PERFTEST_CPP_CSVWRITER_H
#define PERFTEST_CPP_CSVWRITER_H

#include <string>
#include <utility>
#include <fstream>
#include <iostream>

class CSVWriter {
private:
    std::string _file_path;
    std::string _delimiter;
    int _line_count;
    std::fstream _file;
public:
    explicit CSVWriter(std::string filepath, std::string delimiter = ",") : _file_path(std::move(filepath)), _delimiter(std::move(delimiter)), _line_count(0) {
        // Open the file in truncate mode if first line else in Append Mode
        _file.open(_file_path, std::ios::out | (_line_count ? std::ios::app : std::ios::trunc));
    };

    ~CSVWriter(){
        //std::cout << "closing file stream\n";
        // Close the file on destruction
        _file.close();
    }

    template<typename T>
    void addRow(T first, T last);

};

template<typename T>
void CSVWriter::addRow(T first, T last) {
    // Iterate over the range and add each lement to file seperated by delimeter.
    for (; first != last; )
    {
        _file << *first;
        if (++first != last)
            _file << _delimiter;
    }
    _file << "\n";
    _line_count++;

}


#endif //PERFTEST_CPP_CSVWRITER_H
