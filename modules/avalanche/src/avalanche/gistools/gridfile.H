/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Matthias Rauter
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    gridfile

Description
    A class for handling ESRI-gridfiles.

SourceFiles
    gridfile.C

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#ifndef Foam_avalanche_gridfile_H
#define Foam_avalanche_gridfile_H

#include <iostream>

class gridfile
{
public:

    gridfile();

    gridfile(double xllcenter, double yllcenter,
             double dx, double dy,
             int ncols,  int nrows);

    ~gridfile();

    // Logging (eg, stdout)
    std::ostream& log() const;

    void clear();

    std::string info() const;

    int read(const std::string& fname);

    int write(const std::string& fname);

    //Bilinear interpolation
    double interpolate(const double x, const double y) const;

    //Nearest neighbor interpolation
    double interpolateNN(const double x, const double y) const;

    const std::string& filename() const noexcept { return filename_; }

    double xllcenter() const noexcept { return xllcenter_; }

    double yllcenter() const noexcept { return yllcenter_; }

    double dx() const noexcept { return dx_; }

    double dy() const noexcept { return dy_; }

    unsigned int ncols() const noexcept {return ncols_;}

    unsigned int nrows() const noexcept { return nrows_; }

    unsigned int xCellCount() const noexcept { return ncols_; }

    unsigned int yCellCount() const noexcept { return nrows_ ; }

    double v(int xindex, int yindex) const { return v_[xindex][yindex]; }

    double& vRef(int xindex, int yindex) { return v_[xindex][yindex]; }

    inline double** vRef() { return v_; }

    double NODATA_value() const noexcept { return NODATA_value_; }


private:

    // Debug level
    int debug_;

    //Last known filename
    std::string filename_;

    //lower left pixel centre coordinates
    double xllcenter_;
    double yllcenter_;

    //cell size of the raster x-direction
    double dx_;

    //cell size of the raster y-direction
    double dy_;

    //number of columns
    unsigned int ncols_;

    //number of rows
    unsigned int nrows_;

    //proxy value indicating that no data is available
    double NODATA_value_;

    //field values in a nrows x ncols array
    double **v_;
};

#endif

// ************************************************************************* //
