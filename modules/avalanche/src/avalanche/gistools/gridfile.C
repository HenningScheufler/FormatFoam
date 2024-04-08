/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 Matthias Rauter
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

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "gridfile.H"
#include <fstream>
#include <iterator>
#include <sstream>
#include <limits>
#include <cctype>
#include <cmath>
#include <algorithm>

// Debugging output
#undef  DebugInfo
#define DebugInfo  if (debug_) this->log()


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gridfile::gridfile()
:
    debug_(0),
    v_(nullptr)
{
    clear();
}


gridfile::gridfile
(
    double xllcenter,
    double yllcenter,
    double dx,
    double dy,
    int ncols,
    int nrows
)
:
    debug_(0),
    xllcenter_(xllcenter),
    yllcenter_(yllcenter),
    dx_(dx),
    dy_(dy),
    ncols_(ncols),
    nrows_(nrows),
    NODATA_value_(-9999.),
    v_(nullptr)
{
    if (this->nrows_ > 0 && this->ncols_ > 0)
    {
        double * vv = new double [this->nrows_*this->ncols_];
        this->v_ = new double*[this->nrows_];
        for (unsigned int i=0; i<this->nrows_; i++)
        {
            v_[i] = vv+i*ncols_;
        }
        for (unsigned int i=0; i<this->nrows_*this->ncols_; i++)
        {
            vv[i] = this->NODATA_value_;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gridfile::~gridfile()
{
    if (v_ != nullptr)
    {
        delete v_[0];
        delete v_;
    }
    v_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& gridfile::log() const
{
    return std::cout;
}


void gridfile::clear()
{
    this->ncols_ = 0;
    this->nrows_ = 0;
    this->dx_ = NAN;
    this->dy_ = NAN;
    this->xllcenter_ = NAN;
    this->yllcenter_ = NAN;
    this->NODATA_value_ = NAN;

    if (v_ != nullptr)
    {
        delete v_[0];
        delete v_;
    }
    v_ = nullptr;
}


std::string gridfile::info() const
{
    std::ostringstream ss;

    ss << "Gridfile " << this->filename_ << ":" << std::endl
       << "nrows = " << this->nrows_ << std::endl
       << "ncols = " << this->ncols_ << std::endl
       << "dx = " << this->dx_ << std::endl
       << "dy = " << this->dy_ << std::endl
       << "xllcorder = " << this->xllcenter_ << std::endl
       << "yllcorder = " << this->yllcenter_ << std::endl
       << "NODATA_value = " << this->NODATA_value_ << std::endl << std::endl;

    return ss.str();
}


int gridfile::read(const std::string& fname)
{
    clear();

    this->filename_ = fname;

    std::ifstream pFile(filename_);
    if (!pFile.is_open())
    {
        DebugInfo << "Cannot open file \"" << filename_ << "\"" << std::endl;
        return 0;
    }

    std::string line;

    bool readheader = true;

    bool offsetFromCornerToCenter = false;
    while((!pFile.eof()) && readheader )
    {
        std::getline(pFile, line);
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // OR:
        // stringOps::inplaceLower(line);

        DebugInfo << "line=" << line << std::endl;

        std::string::size_type i = std::string::npos;

        if ((i = line.find("ncols")) != std::string::npos)
        {
            this->ncols_ = std::stod(line.substr(i + 5));
        }
        else if ((i = line.find("nrows")) != std::string::npos)
        {
            this->nrows_ = std::stod(line.substr(i + 5));
        }
        else if ((i = line.find("xllcorner")) != std::string::npos)
        {
            offsetFromCornerToCenter = true;
            this->xllcenter_ = std::stod(line.substr(i + 9));
        }
        else if ((i = line.find("xllcenter")) != std::string::npos)
        {
            this->xllcenter_ = std::stod(line.substr(i + 9));
        }
        else if ((i = line.find("yllcorner")) != std::string::npos)
        {
            offsetFromCornerToCenter = true;
            this->yllcenter_ = std::stod(line.substr(i + 9));
        }
        else if ((i = line.find("yllcenter")) != std::string::npos)
        {
            this->yllcenter_ = std::stod(line.substr(i + 9));
        }
        else if ((i = line.find("cellsize")) != std::string::npos)
        {
            this->dx_ = std::stod(line.substr(i + 8));
            this->dy_ = std::stod(line.substr(i + 8));
        }
        else if ((i = line.find("dx")) != std::string::npos)
        {
            this->dx_ = std::stod(line.substr(i + 2));
        }
        else if ((i = line.find("dy")) != std::string::npos)
        {
            this->dy_ = std::stod(line.substr(i + 2));
        }
        else if ((i = line.find("nodata_value")) != std::string::npos)
        {
            this->NODATA_value_ = std::stod(line.substr(i + 12));
            readheader = false;
        }
        else
        {
            DebugInfo << "Skipping line  " << line << std::endl;
        }
    }

    if (this->ncols_ < 1 || this->nrows_ < 1)
    {
        DebugInfo << "Number of columns (" << this->ncols_  << ") or number of rows (" << this->nrows_ << ") invalid." << std::endl;
        return 0;
    }

    if (offsetFromCornerToCenter)
    {
        this->xllcenter_ += this->dx_/2.;
        this->yllcenter_ += this->dy_/2.;
    }

    double * vv = new double [this->nrows_*this->ncols_];
    this->v_ = new double*[this->nrows_];
    for (unsigned int i=0; i<this->nrows_; i++)
    {
        v_[i] = vv+i*ncols_;
    }

    unsigned int i = 0, j = 0;

    while(!pFile.eof() && i < this->nrows_)
    {
        DebugInfo << "Reading line " << i << "/" << this->nrows_ << std::endl;

        std::getline(pFile, line);
        j = 0;
#if 0
        std::regex e("(\\S+)");
        // The following part is very flexible but slow
        while (std::regex_search (line, match, e) and j < this->ncols) {
            try
            {
                this->v[i][j] = std::stod(match[1]);
                j++;
            }
            catch(const std::invalid_argument &ia)
            {
                DebugInfo << "Invalid argument in record: " << ia.what() << std::endl;
            }
            line = match.suffix().str();
        }
#endif

        //This is faster
        std::istringstream ss(line);

        while (j < this->ncols_) {
            double val;
            ss >> val;
            this->v_[i][j] = val;
            j++;
        }
        i++;
    }

    return 1;
}


int gridfile::write(const std::string& fname)
{
    this->filename_ = fname;

    std::ofstream pFile(filename_);
    if (!pFile.is_open())
    {
        return 0;
    }

    pFile << "ncols " << this->ncols_ << std::endl
          << "nrows " << this->nrows_ << std::endl
          << "xllcenter " << this->xllcenter_ << std::endl
          << "yllcenter " << this->yllcenter_ << std::endl;

    if (this->dx_ == this->dy_)
    {
        pFile << "cellsize " << dx_ << std::endl;
    }
    else
    {
        pFile << "dx " << dx_ << std::endl
              << "dy " << dy_ << std::endl;
    }

    pFile << "NODATA_value " << this->NODATA_value_ << std::endl;

    for (unsigned int j=0; j<this->nrows_; j++)
    {
        for (unsigned int i=0; i<this->ncols_; i++)
        {
            //pFile << "(" << i << "," << j << ")";
            pFile << v_[j][i] << " ";
        }
        pFile << std::endl;
    }

    pFile.close();

    return 1;
}

double gridfile::interpolate(const double x, const double y) const
{
    double dx = x-this->xllcenter_;
    double dy = y-this->yllcenter_;

    int itemp = floor(dx/this->dx_);
    int jtemp = ceil(dy/this->dy_);

    unsigned int i = std::min<double>(std::max<double>(itemp, 0), ncols_-1);
    unsigned int j = std::min<double>(std::max<double>(jtemp, 0), nrows_-1);

    unsigned int ipp = std::min<double>(std::max<double>(itemp+1, 0), ncols_-1);
    unsigned int jpp = std::min<double>(std::max<double>(jtemp-1, 0), nrows_-1);

    dx = dx/this->dx_-i;
    dy = dy/this->dy_-j;

    j = this->nrows_-1-j;
    jpp = this->nrows_-1-jpp;
    dy = -dy;

    return v_[j][i]*(1-dy)*(1-dx)+
           v_[jpp][i]*(dy)*(1-dx)+
           v_[j][ipp]*(1-dy)*(dx)+
           v_[jpp][ipp]*(dy)*(dx);
}


double gridfile::interpolateNN(const double x, const double y) const
{
    double dx = x-this->xllcenter_;
    double dy = y-this->yllcenter_;

    int itemp = round(dx/this->dx_);
    int jtemp = round(dy/this->dy_);

    unsigned int i = std::min<double>(std::max<double>(itemp, 0), ncols_-1);
    unsigned int j = std::min<double>(std::max<double>(jtemp, 0), nrows_-1);

    //unsigned int ipp = std::min<double>(std::max<double>(itemp+1, 0), ncols-1);
    //unsigned int jpp = std::min<double>(std::max<double>(jtemp-1, 0), nrows-1);

    //dx = dx/cellsize-i;
    //dy = dy/cellsize-j;

    j = this->nrows_-1-j;
    //jpp = this->nrows-1-jpp;
    //dy = -dy;


    return v_[j][i];
}


// ************************************************************************* //
