/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 Matthias Rauter
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "shapefile.H"
#include <fstream>
#include <iterator>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

// Debugging output
#undef  DebugInfo
#define DebugInfo  if (debug_) this->log()


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Local implementation of ntohl to avoid dependency on
// <arpa/inet.h> or <winsock.h>
//
// Compare with OpenFOAM endian::swap32()

#ifdef ntohl
#undef ntohl
#endif

inline uint32_t ntohl(uint32_t u)
{
#if (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    return u;
#elif (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
# ifdef __GNUC__
    return __builtin_bswap32(u);
# else
    return
    (
        (((u) & 0xff000000) >> 24)      // 3 -> 0
      | (((u) & 0x00ff0000) >>  8)      // 2 -> 1
      | (((u) & 0x0000ff00) <<  8)      // 2 <- 1
      | (((u) & 0x000000ff) << 24)      // 3 <- 0
    );
# endif
#else
# error "__BYTE_ORDER__ is not BIG or LITTLE endian"
#endif
}


template<class T>
std::string vec2string(const std::vector<T>& v)
{
    std::ostringstream ss;

    unsigned i = 0;

    ss << v.size() << '(';
    for (const T& val : v)
    {
        if (i++) ss << ',';
        ss << val;
    }
    ss << ')';

    return ss.str();
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

shapefile::shapefile(type shptype)
{
    debug_ = 0;
    clear();
    this->shptype_ = shptype;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int shapefile::pointcount(int recordIndex, int partIndex) const
{
       if (partIndex < this->partcount(recordIndex)-1)
           return this->partStart(recordIndex, partIndex+1)-this->partStart(recordIndex, partIndex);
       return this->pointcount(recordIndex)-this->partStart(recordIndex, partIndex);
}


int shapefile::addRecord()
{
    this->recordcount_++;

    this->shapeparts_.push_back(0);
    this->shapepoints_.push_back(0);

    this->px_.emplace_back();
    this->py_.emplace_back();
    this->pz_.emplace_back();
    this->pm_.emplace_back();
    this->pn_.emplace_back();
    this->pt_.emplace_back();

    this->bxmin_.push_back(0);
    this->bxmax_.push_back(1);
    this->bymin_.push_back(0);
    this->bymax_.push_back(1);
    this->bzmin_.push_back(0);
    this->bzmax_.push_back(1);
    this->bmmin_.push_back(0);
    this->bmmax_.push_back(1);

    this->v_.emplace_back();
    this->c_.emplace_back();

    return this->recordcount_-1;
}


int shapefile::addPart(int record)
{
    if (record < 0)
        record = this->px_.size()-1;
    this->pn_[record].push_back(this->px_[record].size());
    this->shapeparts_[record]++;


    return this->shapeparts_[record]-1;
}


int shapefile::addPoint(int record, double x, double y, double z, double m)
{
    if (record < 0)
        record = this->px_.size()-1;
    px_[record].push_back(x);
    py_[record].push_back(y);
    pz_[record].push_back(z);
    pm_[record].push_back(m);
    this->shapepoints_[record]++;

    return this->shapepoints_[record]-1;
}


int shapefile::addField(std::string fieldname, char fieldtype, int length, int decimals)
{
    this->fieldcount_++;
    this->fn_.push_back(fieldname);
    this->fl_.push_back(length);
    this->fd_.push_back(decimals);
    this->ft_.push_back(fieldtype);
    if (fieldtype == 'N' || fieldtype == 'F')
    {
        this->isnumeric_.push_back(true);
        for (std::vector<double> &vRef : v_)
        {
            vRef.push_back(100);
        }
        this->fieldIndex_.push_back(this->v_[0].size()-1);
    }
    else
    {
        this->isnumeric_.push_back(false);
        for (std::vector<std::string> &cRef : c_)
        {
            cRef.emplace_back();
        }
        this->fieldIndex_.push_back(this->c_[0].size()-1);
    }
    return this->fieldcount_-1;
}


void shapefile::setField(int recordIndex, int fieldIndex, double value)
{
    if(!this->isnumeric_[fieldIndex])
        throw std::invalid_argument("field is not numeric");
    this->v_[recordIndex][this->fieldIndex_[fieldIndex]] = value;
}


void shapefile::setField(int recordIndex, int fieldIndex, const std::string &value)
{
    if(this->isnumeric_[fieldIndex])
        throw std::invalid_argument("field is numeric");
    this->c_[recordIndex][this->fieldIndex_[fieldIndex]] = value;
}


double shapefile::numericfield(int recordIndex, int fieldIndex) const
{
    if(!this->isnumeric_[fieldIndex])
        throw std::invalid_argument("field is not numeric");
    return this->v_[recordIndex][this->fieldIndex_[fieldIndex]];
}


const std::string& shapefile::stringfield(int recordIndex, int fieldIndex) const
{
    if(this->isnumeric_[fieldIndex])
        throw std::invalid_argument("field is numeric");
    return this->c_[recordIndex][this->fieldIndex_[fieldIndex]];
}


std::ostream& shapefile::log() const
{
    return std::cout;
}


void shapefile::clear()
{
    DebugInfo << "shapefile cleared." << std::endl;

    filename_.clear();

    shptype_ = NULLSHP;
    recordcount_ = 0;

    shapeparts_.clear();
    pn_.clear();
    pt_.clear();

    shapepoints_.clear();
    px_.clear();
    py_.clear();
    pz_.clear();
    pm_.clear();
    bxmin_.clear();
    bxmax_.clear();
    bymin_.clear();
    bymax_.clear();
    bzmin_.clear();
    bzmax_.clear();
    bmmin_.clear();
    bmmax_.clear();

    fieldcount_ = 0;

    fn_.clear();
    ft_.clear();
    fl_.clear();
    fd_.clear();
    v_.clear();
    c_.clear();
    isnumeric_.clear();
    fieldIndex_.clear();

    dbase_datestring_[0] = 0;
    dbase_datestring_[1] = 0;
    dbase_datestring_[2] = 0;

    dbase_versionbyte_ = 2;
}


void shapefile::calcBoundingBox()
{
    xmin_ = std::numeric_limits<double>::infinity();
    ymin_ = std::numeric_limits<double>::infinity();
    zmin_ = std::numeric_limits<double>::infinity();
    mmin_ = std::numeric_limits<double>::infinity();
    xmax_ = -std::numeric_limits<double>::infinity();
    ymax_ = -std::numeric_limits<double>::infinity();
    zmax_ = -std::numeric_limits<double>::infinity();
    mmax_ = -std::numeric_limits<double>::infinity();

    bxmin_.clear();
    bxmax_.clear();
    bymin_.clear();
    bymax_.clear();
    bzmin_.clear();
    bzmax_.clear();
    bmmin_.clear();
    bmmax_.clear();

    for (unsigned int i=0; i<recordcount_; i++)
    {
        for (unsigned int j=0; j<px_[i].size(); j++)
        {
            bxmin_[i] = std::min(bxmin_[i], px_[i][j]);
            bxmax_[i] = std::max(bxmax_[i], px_[i][j]);
            bymin_[i] = std::min(bymin_[i], py_[i][j]);
            bymax_[i] = std::max(bymax_[i], py_[i][j]);

            xmin_ = std::min(xmin_, px_[i][j]);
            xmax_ = std::max(xmax_, px_[i][j]);
            ymin_ = std::min(ymin_, py_[i][j]);
            ymax_ = std::max(ymax_, py_[i][j]);
        }
        for (unsigned int j=0; j<pz_[i].size(); j++)
        {
            bzmin_[i] = std::min(bzmin_[i], pz_[i][j]);
            bzmax_[i] = std::max(bzmax_[i], pz_[i][j]);

            zmin_ = std::min(zmin_, pz_[i][j]);
            zmax_ = std::max(zmax_, pz_[i][j]);
        }
        for (unsigned int j=0; j<pm_[i].size(); j++)
        {
            bmmin_[i] = std::min(bmmin_[i], pm_[i][j]);
            bmmax_[i] = std::max(bmmax_[i], pm_[i][j]);

            mmin_ = std::min(mmin_, pm_[i][j]);
            mmax_ = std::max(mmax_, pm_[i][j]);
        }
    }

    if (recordcount_ == 0)
    {
        xmin_ = ymin_ = zmin_ = mmin_ = xmax_ = ymax_ = zmax_ = mmax_ = 0;
    }
    DebugInfo << "BB ("<<  xmin_ << ", " << ymin_ << ") -- (" << xmax_ << ", " << ymax_ << ")" << std::endl;
}


std::string shapefile::info() const
{
    std::ostringstream ss;

    ss << "Shapefile " << filename_ << ":" << std::endl
       << "Type = " << shptype_ << std::endl
       << "Shape count = " << recordcount_ << std::endl
       << "Parts count = " << vec2string<int>(shapeparts_) << std::endl
       << std::endl;

    for (unsigned int i=0; i<recordcount_; i++)
    {
        ss << "Record " << i << ": " << std::endl
           << "Parts start = " << vec2string<int>(pn_[i]) << std::endl
           << "Part types = " << vec2string<int>(pt_[i]) << std::endl
           << "x-coords = " << vec2string<double>(px_[i]) << std::endl
           << "y-coords = " << vec2string<double>(py_[i]) << std::endl
           << std::endl;
    }

    ss << "Field count = " << fieldcount_ << std:: endl
       << "FieldNames = " << vec2string<std::string>(fn_) << std:: endl
       << "FieldTypes = " << vec2string<char>(ft_) << std:: endl
       << "FieldLenghts = " << vec2string<int>(fl_) << std:: endl
       << "FieldDecimals = " << vec2string<int>(fd_) << std:: endl;


    for (unsigned int i=0; i<fieldcount_; i++)
    {
        ss << "Field " << i << ": " << std::endl
           << "Field Name = " << fn_[i] << std::endl
           << "Field Type = " << ft_[i] << std::endl
           << "Field Length/Decimals = " << fl_[i] << "," << fd_[i] << std::endl
           << std::endl;
    }
    for (unsigned int i=0; i<recordcount_; i++)
    {
        ss << "Field values records #" << i << " = " << vec2string<double>(v_[i]) << std::endl;
        ss << "Field string records #" << i << " = " << vec2string<std::string>(c_[i]) << std::endl;
    }

    return ss.str();
}


int shapefile::read(std::string filename)
{
    clear();

    this->filename_ = filename;

    std::vector<int> recordPositions;

    std::ifstream input;

    input.open(filename+".shx", std::ios::binary);
    if (!input.is_open())
    {
        return 0;
    }
    std::vector<char>buffer(std::istreambuf_iterator<char>(input), {});
    char* c = &buffer[0];

    int bufferlength = buffer.size();
    DebugInfo << "\n\nShx-file is " << bufferlength << " byte long." << std::endl;

    int pos = 100;

    do
    {
        uint32_t recordoffset = ntohl(*reinterpret_cast<uint32_t*>(c+pos));
        uint32_t recordlength = ntohl(*reinterpret_cast<uint32_t*>(c+pos+4));
        pos += 8;

        DebugInfo
            << "Recordoffset: " << recordoffset << std::endl
            << "Recordlength: " << recordlength << std::endl
            << "read " << pos << "/" << bufferlength << std::endl;

        recordPositions.push_back(recordoffset);
    }
    while(pos < bufferlength);

    input.close();

    input.open(filename+".shp", std::ios::binary);
    if (!input.is_open())
    {
        return 0;
    }
    buffer = std::vector<char>(std::istreambuf_iterator<char>(input), {});
    c = &buffer[0];

    bufferlength = buffer.size();

    DebugInfo << "\n\nShp-file is " << bufferlength << " byte long." << std::endl;

    uint32_t filecode = ntohl(*reinterpret_cast<uint32_t*>(c));
    uint32_t filelength = ntohl(*reinterpret_cast<uint32_t*>(c+24));
    uint32_t version = *reinterpret_cast<uint32_t*>(c+28);
    this->shptype_ = *reinterpret_cast<uint32_t*>(c+32);

    this->xmin_ = *reinterpret_cast<double*>(c+36);
    this->ymin_ = *reinterpret_cast<double*>(c+44);
    this->xmax_ = *reinterpret_cast<double*>(c+52);
    this->ymax_ = *reinterpret_cast<double*>(c+60);
    this->zmin_ = *reinterpret_cast<double*>(c+68);
    this->zmax_ = *reinterpret_cast<double*>(c+76);
    this->mmin_ = *reinterpret_cast<double*>(c+84);
    this->mmax_ = *reinterpret_cast<double*>(c+92);

    DebugInfo
        << "Filecode: " << filecode << std::endl
        << "Filelength: " << filelength << std::endl
        << "Version: " << version << std::endl
        << "Shypetype: " << shptype_ << std::endl
        << "BB ("
            << this->xmin_ << ", " << this->ymin_ << ", "
            << this->zmin_ << ", " << this->mmin_ << ") -- ("
            << this->xmax_ << ", " << this->ymax_ << ", "
            << this->zmin_ << ", " << this->mmin_ << ")"
        << std::endl;

    unsigned int i = 0;
    pos = 100;

    for(i = 0; i<recordPositions.size(); i++)
    {
        pos = recordPositions[i]*2;

        char* pos2 = c+pos;

        uint32_t recordnumber = ntohl(*reinterpret_cast<uint32_t*>(pos2));
        pos2 += 4;
        uint32_t recordlength = ntohl(*reinterpret_cast<uint32_t*>(pos2));
        pos2 += 4;

        uint32_t shptype2 = *reinterpret_cast<uint32_t*>(pos2);
        pos2 += 4;

        DebugInfo << "Shapetype (global/local): " << this->shptype_ << "/" << shptype2 << std::endl;
        if (this->shptype_ != shptype2)
        {
            log() << "Uexpected shapetype. Ignoring local shapetype";
        }

        DebugInfo
            << std::endl
            << "Recordnumber: " << recordnumber << std::endl
            << "Recordoffset: " << pos/2 << std::endl
            << "Recordlength: " << recordlength << std::endl
            << "Current pos is " << pos << std::endl;

        if (recordlength < 3){

            DebugInfo << "Recording is empty - not reading further." << std::endl;

            this->px_.emplace_back();
            this->py_.emplace_back();
            this->pn_.emplace_back();
            this->shapeparts_.push_back(0);
            this->shapepoints_.push_back(0);

            continue;
        }

        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> m;
        std::vector<int> parts;
        std::vector<int> parttypes;

        if (this->shptype_ == NULLSHP)
        {
            this->shapeparts_.push_back(1);
            this->shapepoints_.push_back(0);
        }
        else if (this->shptype_ == POINT ||
                 this->shptype_ == POINTM ||
                 this->shptype_ == POINTZ)
        {
            this->shapeparts_.push_back(1);
            this->shapepoints_.push_back(1);

            double xx = *reinterpret_cast<double*>(pos2);
            pos2 += 8;
            double yy = *reinterpret_cast<double*>(pos2);
            pos2 += 8;

            DebugInfo << "p(" << xx << ", " << yy << ")" << std::endl;

            x.push_back(xx);
            y.push_back(yy);

            if (this->shptype_ == POINTZ)
            {
                double zz = *reinterpret_cast<double*>(pos2);
                pos2 += 8;

                DebugInfo << ", z " << zz << std::endl;
                z.push_back(yy);
            }

            if (this->shptype_ == POINTM ||
                this->shptype_ == POINTZ)
            {
                double mm = *reinterpret_cast<double*>(pos2);
                pos2 += 8;

                DebugInfo << "m " << mm << std::endl;
                m.push_back(mm);
            }
        }
        else if (this->shptype_ == MULTIPOINT ||
                 this->shptype_ == MULTIPOINTM ||
                 this->shptype_ == MULTIPOINTZ)
        {
            this->bxmin_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bymin_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bxmax_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bymax_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;

            uint32_t numpoints = *reinterpret_cast<uint32_t*>(pos2);
            pos2 += 4;

            DebugInfo << "# Point: " << numpoints << std::endl;

            this->shapeparts_.push_back(1);
            this->shapepoints_.push_back(numpoints);

            for (unsigned int k=0; k<numpoints; k++)
            {
                double xx = *reinterpret_cast<double*>(pos2);
                pos2 += 8;
                double yy = *reinterpret_cast<double*>(pos2);
                pos2 += 8;

                DebugInfo << "p(" << xx << ", " << yy << ")" << std::endl;

                x.push_back(xx);
                y.push_back(yy);
            }

            if (this->shptype_ == MULTIPOINTZ)
            {
                this->bzmin_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;
                this->bzmax_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;

                for (unsigned int k=0; k<numpoints; k++)
                {
                    double zz = *reinterpret_cast<double*>(pos2);
                    pos2 += 8;

                    DebugInfo << "z " << zz << std::endl;

                    z.push_back(zz);
                }
            }

            if (this->shptype_ == MULTIPOINTM ||
                this->shptype_ == MULTIPOINTZ)
            {
                this->bmmin_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;
                this->bmmax_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;

                for (unsigned int k=0; k<numpoints; k++)
                {
                    double mm = *reinterpret_cast<double*>(pos2);
                    pos2 += 8;

                    DebugInfo << "m " << mm << std::endl;

                    m.push_back(mm);
                }
            }
        }
        else if (this->shptype_ == POLYLINE ||
                 this->shptype_ == POLYLINEM ||
                 this->shptype_ == POLYLINEZ ||
                 this->shptype_ == POLYGON ||
                 this->shptype_ == POLYGONM ||
                 this->shptype_ == POLYGONZ ||
                 this->shptype_ == MULTIPATCH)
        {
            // these types are all very similar and
            // therefore parsed with the same code

            this->bxmin_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bymin_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bxmax_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;
            this->bymax_.push_back(*reinterpret_cast<double*>(pos2));
            pos2 += 8;

            uint32_t numparts = *reinterpret_cast<uint32_t*>(pos2);
            pos2 += 4;
            uint32_t numpoints = *reinterpret_cast<uint32_t*>(pos2);
            pos2 += 4;

            DebugInfo
                << "# Parts: " << numparts << std::endl
                << "# Point: " << numpoints << std::endl;

            this->shapeparts_.push_back(numparts);
            this->shapepoints_.push_back(numpoints);

            for (unsigned int k=0; k<numparts; k++)
            {
                int j = *reinterpret_cast<uint32_t*>(pos2);
                pos2 += 4;

                DebugInfo << "Part start: " << j << std::endl;

                parts.push_back(j);
            }

            if (this->shptype_ == MULTIPATCH)
            {
                for (unsigned int k=0; k<numparts; k++)
                {
                    int j = *reinterpret_cast<uint32_t*>(pos2);
                    pos2 += 4;

                    DebugInfo << "Part type: " << j << std::endl;

                    parttypes.push_back(j);
                }
            }

            for (unsigned int k=0; k<numpoints; k++)
            {
                double xx = *reinterpret_cast<double*>(pos2);
                pos2 += 8;
                double yy = *reinterpret_cast<double*>(pos2);
                pos2 += 8;

                DebugInfo << "p(" << xx << ", " << yy << ")" << std::endl;

                x.push_back(xx);
                y.push_back(yy);
            }

            if (this->shptype_ == POLYLINEZ ||
                this->shptype_ == POLYGONZ ||
                this->shptype_ == MULTIPATCH)
            {
                this->bzmin_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;
                this->bzmax_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;

                for (unsigned int k=0; k<numpoints; k++)
                {
                    double zz = *reinterpret_cast<double*>(pos2);
                    pos2 += 8;

                    DebugInfo << "z " << zz << std::endl;

                    z.push_back(zz);
                }
            }
            if (this->shptype_ == POLYLINEZ ||
                this->shptype_ == POLYLINEM ||
                this->shptype_ == POLYGONZ ||
                this->shptype_ == POLYGONM ||
                this->shptype_ == MULTIPATCH)
            {
                this->bmmin_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;
                this->bmmax_.push_back(*reinterpret_cast<double*>(pos2));
                pos2 += 8;

                for (unsigned int k=0; k<numpoints; k++)
                {
                    double mm = *reinterpret_cast<double*>(pos2);
                    pos2 += 8;

                    DebugInfo << "m " << mm << std::endl;

                    m.push_back(mm);
                }
            }
        }
        else
        {
            //Not implemented shape
            DebugInfo << "Shape type \"" << shptype_ << "\" not implemented - not reading further." << std::endl;
            return 1;
        }

        this->px_.push_back(x);
        this->py_.push_back(y);
        this->pz_.push_back(z);
        this->pm_.push_back(m);
        this->pn_.push_back(parts);
        this->pt_.push_back(parttypes);

        pos += 8+recordlength*2;
        if (pos > bufferlength -2)
        {
            break;
        }

        DebugInfo << "read " << pos << "/" << bufferlength << std::endl;
    }

    this->recordcount_ = recordPositions.size();

    DebugInfo << "\nFinished reading shp-file. Found " << this->recordcount_ << ". Closing file.\n" << std::endl;

    input.close();


    input.open(filename+".dbf", std::ios::binary);
    if (!input.is_open())
    {
        return 0;
    }
    buffer = std::vector<char>(std::istreambuf_iterator<char>(input), {});
    c = &buffer[0];

    bufferlength = buffer.size();
    DebugInfo << "\n\nDbf-file is " << bufferlength << " byte long." << std::endl;

    this->dbase_versionbyte_ = c[0];
    this->dbase_datestring_[0] = c[1];
    this->dbase_datestring_[1] = c[2];
    this->dbase_datestring_[2] = c[3];

    uint32_t recordnumber = *reinterpret_cast<uint32_t*>(c+4);
    uint16_t headerlength = *reinterpret_cast<uint16_t*>(c+8);
    uint16_t recordlength = *reinterpret_cast<uint16_t*>(c+10);

    DebugInfo
        << "Recordnumber: " << recordnumber << std::endl
        << "Headerlength: " << headerlength << std::endl
        << "Recordlenght: " << recordlength << std::endl;


    i = 0;
    pos = 32;

    int ic = 0;
    int in = 0;
    do
    {
        char fieldname[12] = {'\0'};
        char fieldtype[2] = {'\0'};

        strncpy(fieldname, c+pos, 11);
        strncpy(fieldtype, c+pos+11, 1);
        DebugInfo
            << "\n\n\nFieldname: " << fieldname << std::endl
            << "Fieldtype: " << fieldtype << std::endl;


        uint32_t fieldpos = *reinterpret_cast<uint32_t*>(c+pos+12);
        uint16_t fieldlength = *reinterpret_cast<uint8_t*>(c+pos+16);
        uint16_t fielddec = *reinterpret_cast<uint8_t*>(c+pos+17);

        DebugInfo
            << "Fieldpos: " << fieldpos << std::endl
            << "Fieldlength: " << fieldlength << std::endl
            << "Fielddec: " << fielddec << std::endl;


        uint16_t termbit = *reinterpret_cast<uint8_t*>(c+pos+32);
        DebugInfo << "termbit: " << termbit << std::endl;

        pos += 32;

        this->fn_.push_back(std::string(fieldname));
        this->ft_.push_back(fieldtype[0]);
        this->fl_.push_back(fieldlength);
        this->fd_.push_back(fielddec);
        if (fieldtype[0] == 'N' || fieldtype[0] == 'F')
        {
            this->isnumeric_.push_back(true);
            this->fieldIndex_.push_back(in++);
        }
        else
        {
            this->fieldIndex_.push_back(false);
            this->fieldIndex_.push_back(ic++);
        }
        i++;

        if (termbit == 13)
            break;
    }
    while(pos < bufferlength);

    this->fieldcount_ = i;

    DebugInfo
        << std::endl << std::endl
        << "Fieldcount: " << fieldcount_ << std::endl;


    pos = headerlength;

    i = 0;
    do
    {
        pos += 1;

        std::vector<double> values;
        std::vector<std::string> strings;

        DebugInfo
            << std::endl
            << "Reading record #" << i << std::endl << std::endl;

        for(unsigned int j=0; j<fieldcount_; j++)
        {

            int buflen = fl_[j];

            //FIX from Mark
            std::string bufstring;
            bufstring.resize(buflen+2);
            char* buf = &(bufstring[0]);
            //char buf[buflen+2] = {'\0'};

            strncpy(buf, c+pos, buflen);

            pos += fl_[j];

            if (this->ft_[j] == 'N' || this->ft_[j] == 'F')
            {
                double val = NAN;
                try
                {
                    val = std::stod(buf);
                }
                catch(const std::invalid_argument &ia)
                {
                    DebugInfo << "Invalid argument in record: " << ia.what() << std::endl;
                }
                strings.emplace_back();
                values.push_back(val);

                DebugInfo << fn_[j] << " = " << "\"" << buf << "\"" << "/" << values[values.size()-1] << std::endl;
            }
            else if (this->ft_[j] == 'C')
            {
                strings.emplace_back(buf);
                values.push_back(NAN);

                DebugInfo << fn_[j] << " = " << "\"" << buf << "\"" << "/" << strings[strings.size()-1] << std::endl;
            }
            else if (this->ft_[j] == 'D')
            {
                //Date in format YYYYMMDD
                strings.emplace_back(buf);
                values.push_back(NAN);

                DebugInfo << fn_[j] << " = " << "\"" << buf << "\"" << "/" << strings[strings.size()-1] << std::endl;
            }
            else if (this->ft_[j] == 'L')
            {
                //Boolean value YyTt -> 1 NnFf -> 0
                strings.emplace_back(buf);
                if (buf[0] == 'Y' || buf[0] == 'y' || buf[0] == 'T' || buf[0] == 't')
                    values.push_back(1);
                else
                    values.push_back(0);
            }
            else
            {
                DebugInfo << "Value type \"" << this->ft_[j] << "\" not implemented - not reading further." << std::endl;
                return 0;
            }
        }

        this->v_.push_back(values);
        this->c_.push_back(strings);
        i++;

        DebugInfo << "read " << pos << "/" << bufferlength << std::endl;

    }
    while(i <= recordcount_ && pos < bufferlength);

    input.close();

    DebugInfo << "Finished reading file." << std::endl;
    return 1;
}


int shapefile::write(std::string filename)
{
    //calcBoundingBox();

    this->filename_ = filename;

    std::ofstream output;
    output.open(filename+".shp", std::ios::binary);

    int totlength = 100;

    std::vector<uint32_t> recordLengths;
    std::vector<uint32_t> recordPos;

    for (unsigned int i=0; i<this->recordcount_; i++)
    {
        if (this->shptype_ == NULLSHP)
        {
            recordLengths.push_back(0);
        }
        else if (this->shptype_ == POINT)
        {
            recordLengths.push_back(2*8);//x,y
        }
        else if (this->shptype_ == POINTM)
        {
            recordLengths.push_back(3*8); //x,y,m
        }
        else if (this->shptype_ == POINTZ)
        {
            recordLengths.push_back(4*8); //x,y,z,m
        }
        else if (this->shptype_ == MULTIPOINT)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +4 //npoints
            +2*8*this->shapepoints_[i]); //x,y
        }
        else if (this->shptype_ == MULTIPOINTM)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +4 //npoints
            +2*8*this->shapepoints_[i] //x,y
            +2*8 //mmin, mmax
            +8*this->shapepoints_[i]); //m
        }
        else if (this->shptype_ == MULTIPOINTZ)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +4 //npoints
            +2*8*this->shapepoints_[i] //x,y
            +2*8 //zmin, zmax
            +8*this->shapepoints_[i] //z
            +2*8 //mmin, mmax
            +8*this->shapepoints_[i]); //m
        }
        else if (this->shptype_ == POLYLINE || this->shptype_ == POLYGON)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +2*4 //nparts,npoints
            +4*this->shapeparts_[i] //parts
            +2*8*this->shapepoints_[i]); //x,y
        }
        else if (this->shptype_ == POLYLINEM || this->shptype_ == POLYGONM)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +2*4 //nparts,npoints
            +4*this->shapeparts_[i] //parts
            +2*8*this->shapepoints_[i] //x,y
            +2*8 //mmin, mmax
            +8*this->shapepoints_[i]); //m
        }
        else if (this->shptype_ == POLYLINEZ || this->shptype_ == POLYGONZ)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +2*4 //nparts,npoints
            +4*this->shapeparts_[i] //parts
            +2*8*this->shapepoints_[i] //x,y
            +2*8 //zmin, zmax
            +8*this->shapepoints_[i] //z
            +2*8 //mmin, mmax
            +8*this->shapepoints_[i]); //m
        }
        else if (this->shptype_ == MULTIPATCH)
        {
            recordLengths.push_back(4*8 //xmin, ymin, xmax, ymax
            +2*4 //nparts,npoints
            +4*this->shapeparts_[i] //parts
            +4*this->shapeparts_[i] //parttypes
            +2*8*this->shapepoints_[i] //x,y
            +2*8 //zmin, zmax
            +8*this->shapepoints_[i] //z
            +2*8 //mmin, mmax
            +8*this->shapepoints_[i]); //m
        }
        else
        {
            DebugInfo << "Shape type \"" << this->shptype_ << "\" not implemented - not writing further." << std::endl;
            return 0;
        }
    }


    for (unsigned int i=0; i<this->recordcount_; i++)
    {
        totlength += 4+4+4+recordLengths[i];
        recordLengths[i] = recordLengths[i]/2+4/2;
        if (i > 0)
            recordPos.push_back(recordPos[i-1]+recordLengths[i-1]+8/2);
        else
            recordPos.push_back(100/2);
    }


    uint32_t ft(ntohl(9994));
    uint32_t empty(0);
    uint32_t fl(ntohl(totlength/2));
    uint32_t version(1000);

    // 100 first bytes of header //

    output.write(reinterpret_cast<char*>(&ft), 4); //file code
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&fl), 4); //file length
    output.write(reinterpret_cast<char*>(&version), 4); //version
    output.write(reinterpret_cast<char*>(&this->shptype_), 4); //shapetype
    output.write(reinterpret_cast<char*>(&this->xmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->ymin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->xmax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->ymax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->zmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->zmax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->mmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->mmax_), 8); //boundingbox


    for (unsigned int i=0; i<this->recordcount_; i++)
    {
        uint32_t recordnumber(ntohl(i+1));
        uint32_t contentlength(ntohl(recordLengths[i]));

        output.write(reinterpret_cast<char*>(&recordnumber), 4);  //recordnumber
        output.write(reinterpret_cast<char*>(&contentlength), 4);  //length of the record in byte
        output.write(reinterpret_cast<char*>(&this->shptype_), 4);  //shapetype (again)

        if (this->shptype_ == NULLSHP)
        {
            //Nothing to write for NULLSHP
        }
        else if (this->shptype_ == POINT ||
                 this->shptype_ == POINTM ||
                 this->shptype_ == POINTZ)
        {
            output.write(reinterpret_cast<char*>(&this->px_[i][0]), 8);  //x-coordinate point
            output.write(reinterpret_cast<char*>(&this->py_[i][0]), 8);  //y-coordinate point

            if (this->shptype_ == POINTZ)
            {
                output.write(reinterpret_cast<char*>(&this->pz_[i][0]), 8);  //z-value point
            }

            if (this->shptype_ == POINTM ||
                this->shptype_ == POINTZ)
            {
                output.write(reinterpret_cast<char*>(&this->pm_[i][0]), 8);  //m-value point
            }
        }
        else if (this->shptype_ == MULTIPOINT ||
                 this->shptype_ == MULTIPOINTM ||
                 this->shptype_ == MULTIPOINTZ)
        {
            output.write(reinterpret_cast<char*>(&this->bxmin_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bymin_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bxmax_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bymax_[i]), 8);  //boundingbox


            output.write(reinterpret_cast<char*>(&shapepoints_[i]), 4);  //number of parts
            for (int j=0; j<this->shapepoints_[i]; j++)
            {
                output.write(reinterpret_cast<char*>(&this->px_[i][j]), 8);  //x-coordinate point
                output.write(reinterpret_cast<char*>(&this->py_[i][j]), 8);  //y-coordinate point
            }

            if (this->shptype_ == MULTIPOINTZ)
            {
                output.write(reinterpret_cast<char*>(&this->bzmin_[i]), 8);  //zmin
                output.write(reinterpret_cast<char*>(&this->bzmax_[i]), 8);  //zmax
                for (int j=0; j<this->shapepoints_[i]; j++)
                {
                    output.write(reinterpret_cast<char*>(&this->pz_[i][j]), 8);  //z-coordinate point
                }
            }

            if (this->shptype_ == MULTIPOINTM ||
                this->shptype_ == MULTIPOINTZ)
            {
                output.write(reinterpret_cast<char*>(&this->bmmin_[i]), 8);  //mmin
                output.write(reinterpret_cast<char*>(&this->bmmax_[i]), 8);  //mmax
                for (int j=0; j<this->shapepoints_[i]; j++)
                {
                    output.write(reinterpret_cast<char*>(&this->pm_[i][j]), 8);  //m value
                }
            }

        }
        else if (this->shptype_ == POLYLINE ||
                 this->shptype_ == POLYLINEM ||
                 this->shptype_ == POLYLINEZ ||
                 this->shptype_ == POLYGON ||
                 this->shptype_ == POLYGONM ||
                 this->shptype_ == POLYGONZ ||
                 this->shptype_ == MULTIPATCH)
        {
            output.write(reinterpret_cast<char*>(&this->bxmin_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bymin_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bxmax_[i]), 8);  //boundingbox
            output.write(reinterpret_cast<char*>(&this->bymax_[i]), 8);  //boundingbox

            output.write(reinterpret_cast<char*>(&shapeparts_[i]), 4);  //number of parts
            output.write(reinterpret_cast<char*>(&shapepoints_[i]), 4);  //number of parts

            for (int j=0; j<this->shapeparts_[i]; j++)
            {
                output.write(reinterpret_cast<char*>(&this->pn_[i][j]), 4);  //start of part
            }

            if (this->shptype_ == MULTIPATCH)
            {
                for (int j=0; j<this->shapeparts_[i]; j++)
                {
                    output.write(reinterpret_cast<char*>(&this->pt_[i][j]), 4);  //part type
                }
            }

            for (int j=0; j<this->shapepoints_[i]; j++)
            {
                output.write(reinterpret_cast<char*>(&this->px_[i][j]), 8);  //x-coordinate point
                output.write(reinterpret_cast<char*>(&this->py_[i][j]), 8);  //y-coordinate point
            }

            if (this->shptype_ == POLYLINEZ ||
                this->shptype_ == POLYGONZ ||
                this->shptype_ == MULTIPATCH)
            {
                output.write(reinterpret_cast<char*>(&this->bzmin_[i]), 8);  //zmin
                output.write(reinterpret_cast<char*>(&this->bzmax_[i]), 8);  //zmax
                for (int j=0; j<this->shapepoints_[i]; j++)
                {
                    output.write(reinterpret_cast<char*>(&this->pz_[i][j]), 8);  //z-coordinate point
                }
            }

            if (this->shptype_ == POLYLINEM ||
                this->shptype_ == POLYLINEZ ||
                this->shptype_ == POLYGONM ||
                this->shptype_ == POLYGONZ ||
                this->shptype_ == MULTIPATCH)
            {
                output.write(reinterpret_cast<char*>(&this->bmmin_[i]), 8);  //mmin
                output.write(reinterpret_cast<char*>(&this->bmmax_[i]), 8);  //mmax
                for (int j=0; j<this->shapepoints_[i]; j++)
                {
                    output.write(reinterpret_cast<char*>(&this->pm_[i][j]), 8);  //m value
                }
            }

        }
        else
        {
            DebugInfo << "Shape type \"" << this->shptype_ << "\" not implemented - not writing further." << std::endl;
            return 0;
        }

    }
    output.close();

    output.open(filename+".shx", std::ios::binary);

    totlength = 100+this->recordcount_*8;
    fl = ntohl(totlength/2);
    output.write(reinterpret_cast<char*>(&ft), 4); //file code
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&empty), 4); //unused
    output.write(reinterpret_cast<char*>(&fl), 4); //file length
    output.write(reinterpret_cast<char*>(&version), 4); //version
    output.write(reinterpret_cast<char*>(&this->shptype_), 4); //shapetype
    output.write(reinterpret_cast<char*>(&this->xmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->ymin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->xmax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->ymax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->zmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->zmax_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->mmin_), 8); //boundingbox
    output.write(reinterpret_cast<char*>(&this->mmax_), 8); //boundingbox

    for (unsigned int i=0; i<recordLengths.size(); i++)
    {
        DebugInfo << "Writing indices: " << recordPos[i] << "/" << recordLengths[i] << std::endl;

       recordPos[i] = ntohl(recordPos[i]);
       recordLengths[i] = ntohl(recordLengths[i]);
       output.write(reinterpret_cast<char*>(&recordPos[i]), 4);
       output.write(reinterpret_cast<char*>(&recordLengths[i]), 4);
   }

    output.close();


    output.open(filename+".dbf", std::ios::binary);

    uint32_t numrecords(this->recordcount_);
    uint16_t headersize(33 + this->fieldcount_*32);

    uint16_t recordsize(1);

    for (unsigned int i=0; i<this->fieldcount_; i++)
    {
        recordsize += this->fl_[i];
    }

    output.write(&dbase_versionbyte_, 1); //Version byte
    output.write(dbase_datestring_, 3); //Date string
    output.write(reinterpret_cast<char*>(&numrecords), 4); //Number of records in table
    output.write(reinterpret_cast<char*>(&headersize), 2); //Number of bytes in header
    output.write(reinterpret_cast<char*>(&recordsize), 2); //Number of bytes in a record

    char filler[20] = {'\0'};
    char spacefiller[20] = {' '};

    output.write(filler, 20); //unused

    for (unsigned int i=0; i<this->fieldcount_; i++)
    {
        char fieldName[12] = {' '};
        strncpy(fieldName, this->fn_[i].c_str(), 11);

        uint8_t fieldlength(this->fl_[i]);
        uint8_t fielddecimal(this->fd_[i]);

        output.write(fieldName, 11); //field name
        output.write(&this->ft_[i], 1); //field type
        output.write(filler, 4); //address in memory
        output.write(reinterpret_cast<char*>(&fieldlength), 1); //field length
        output.write(reinterpret_cast<char*>(&fielddecimal), 1); //field decimal count
        output.write(filler, 14); //unused
    }

    uint8_t termbit(13);
    output.write(reinterpret_cast<char*>(&termbit), 1);//termination bit == 13

    for (unsigned int i=0; i<this->recordcount_; i++)
    {
        output.write(spacefiller, 1); //write all values

        for (unsigned int j=0; j<this->fieldcount_; j++)
        {
            if (this->ft_[j] == 'N' || this->ft_[j] == 'F')
            {
                int buflen = this->fl_[j];
                //FIX from Mark:
                std::string bufstring;
                bufstring.resize(buflen+1);
                char* buf = &(bufstring[0]);
                //char buf[buflen+1] = {'\0'};

                char formatstr[7] = {' '};
                formatstr[0] = '%';

                if (std::isnan(v_[i][j]))
                {
                    for (int k=0; k<buflen; k++)
                        buf[k] = '*';
                }
                else if(this->fd_[j] > 0)
                {
                    // TDB: use string concatenation?
                    // eg,
                    // std::string fmtstr
                    // (
                    //     "%" + std::to_string(this->fl_[j])
                    //   + "." + std::to_string(this->fd_[j]) + "f"
                    // );

                    if (snprintf(formatstr+1, 6, "%d.%df", this->fl_[j], this->fd_[j]) < 0)
                    {
                        //should never happen but prevents gcc warning
                        DebugInfo << "Format string exceeds 6 characters - not writing further." << std::endl;
                        return 0;
                    }
                    snprintf(buf, buflen+1, formatstr, v_[i][j]);
                }
                else
                {
                    snprintf(formatstr+1, 6, "%dd", this->fl_[j]);
                    snprintf(buf, buflen+1, formatstr, int(v_[i][j]));
                }

                DebugInfo
                    << formatstr << std::endl
                    << v_[i][j] << "->" << "\"" << buf << "\"" << std::endl;

                /*for(int k=0; k<buflen; k++)
                {
                    if (buf[k] == ' ')
                        buf[k] = '0';
                    if (buf[k] == '\0')
                        buf[k] = '0';
                }*/
                DebugInfo << v_[i][j] << "->" << "\"" << buf << "\"" << std::endl;

                output.write(buf, buflen); //write all values
            }
            else if (this->ft_[j] == 'C' || this->ft_[j] == 'D')
            {
                int buflen = this->fl_[j];
                //FIX from Mark
                std::string bufstring;
                bufstring.resize(buflen+1);
                char* buf = &(bufstring[0]);
                //char buf[buflen+1] = {'\0'};

                DebugInfo
                    << c_[i].size() << std::endl
                    << c_[i][j] << std::endl;

                strncpy(buf, c_[i][j].c_str(), this->fl_[j]);
                output.write(buf, buflen); //write all values
            }
            else
            {
                DebugInfo << "Value type \"" << this->ft_[j] << "\" not implemented - not reading further." << std::endl;
                return 0;
            }
        }
    }

    output.close();

    return 1;
}


// ************************************************************************* //
