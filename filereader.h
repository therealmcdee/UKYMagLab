#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "geometry.h"

class FileReader
{
private:
    std::string path;
    std::ifstream infile;
    std::string outputMessage;
    std::vector<Vec3d> coordinateList;

public:
    bool isSuccessful;

    FileReader(const std::string& n_path) : path(n_path)
    {
        infile.open(path);
        if (infile.fail())
        {
            outputMessage = "Failed to load : " + path;
            isSuccessful = false;
        }
        else
        {
            while (!infile.eof())
            {
                Vec3d temp;

                infile >> temp.x;
                infile >> temp.y;
                infile >> temp.z;
                coordinateList.push_back(temp);

            }
            outputMessage = "File loaded : " + path;
            isSuccessful = true;
        }


        infile.close();

    }

    std::vector<Vec3d> getCoordinates()
    {
        return coordinateList;
    }

    std::string getOutput()
    {
        return outputMessage;
    }
};
