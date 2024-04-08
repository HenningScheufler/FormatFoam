/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | avalanche module
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Matthias Rauter
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

\*---------------------------------------------------------------------------*/

#include "isoLine.H"
#include "Time.H"
#include "areaFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "shapefile.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(isoLine, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        isoLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::isoLine::isoLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    #if (OPENFOAM < 2310)
    aMesh_(obr_.lookupObject<faMesh>("faMesh")),
    #else
    aMesh_(faMesh::mesh(refCast<const polyMesh>(obr_))),
    #endif
    fileName_(),
    fieldName_(),
    isoValues_(),
    offset_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::isoLine::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    dict.readEntry("fileName", fileName_);
    dict.readEntry("field", fieldName_);
    dict.readEntry("values", isoValues_);

    offset_ = Zero;
    dict.readIfPresent("offset", offset_);

    return true;
}


bool Foam::functionObjects::isoLine::execute()
{
    return true;
}


bool Foam::functionObjects::isoLine::write()
{
    if (!obr_.time().writeTime())
    {
        return false;
    }

    Info << type() << " " << name() << " write:" << nl;

    const areaScalarField& phi =
        aMesh_.thisDb().lookupObject<areaScalarField>(fieldName_);

    Info<< "field " << fieldName_
        << ": limits: " << min(phi).value() << ", " << max(phi).value() << nl;


    const labelList& own = aMesh_.edgeOwner();
    const labelList& nei = aMesh_.edgeNeighbour();
    const edgeList& edges = aMesh_.edges();

    vectorList edgePoint(aMesh_.nEdges(), Zero);

    labelList ncc(aMesh_.nFaces(), Zero);
    labelList npc(aMesh_.nPoints(), Zero);

    forAll(nei, edgei)
    {
        ncc[nei[edgei]]++;
    }
    forAll(own, edgei)
    {
        ncc[own[edgei]]++;

        npc[edges[edgei][0]]++;
        npc[edges[edgei][1]]++;
    }

    labelListList faceEdges(aMesh_.nFaces());
    labelListList pointEdges(aMesh_.nPoints());

    forAll(faceEdges, facei)
    {
        faceEdges[facei].resize(ncc[facei]);
    }
    ncc = 0;

    forAll(pointEdges, pointi)
    {
        pointEdges[pointi].resize(npc[pointi]);
    }
    npc = 0;

    forAll(nei, edgei)
    {
        label neiFacei = nei[edgei];
        faceEdges[neiFacei][ncc[neiFacei]++] = edgei;
    }

    forAll(own, edgei)
    {
        label ownFacei = own[edgei];
        faceEdges[ownFacei][ncc[ownFacei]++] = edgei;

        label point1i = edges[edgei][0];
        label point2i = edges[edgei][1];

        pointEdges[point1i][npc[point1i]++] = edgei;
        pointEdges[point2i][npc[point2i]++] = edgei;
    }

    //- The shapefile which is written
    shapefile shp(shapefile::POLYLINEZ);

    DynamicList<scalar> isoValuesForRecords;

    forAll(isoValues_, isoi)
    {
        scalar isoValue = isoValues_[isoi];
        scalarList isoV(aMesh_.nEdges(), -1);

        Info<< "creating isoLine for " << isoValue << endl;

        forAll(own, edgei)
        {
            if (edgei < nei.size())
            {
                scalar a = phi.internalField()[own[edgei]];
                scalar b = phi.internalField()[nei[edgei]];
                if (a-b != 0)
                {
                    isoV[edgei] = (a-isoValue)/(a-b);
                }
                else if (a == isoValue)
                {
                    isoV[edgei] = 0.5;
                }

                edgePoint[edgei] = aMesh_.areaCentres()[nei[edgei]]*isoV[edgei]+aMesh_.areaCentres()[own[edgei]]*(1.-isoV[edgei]);
            }
            else
            {
                if (phi.internalField()[own[edgei]] > isoValue)
                {
                    isoV[edgei] = 0.5;
                }
                edgePoint[edgei] =  0.5*(aMesh_.points()[edges[edgei][0]]+aMesh_.points()[edges[edgei][1]]);
            }
        }

        forAll(nei, edgei)
        {
            if (isoV[edgei] > 0 && isoV[edgei] <= 1)
            {
                vectorList p;
                p.append(edgePoint[edgei]);

                isoV[edgei] = -1;

                vector n(1,0,0);

                label edgej = edgei;
                do
                {
                    labelList candidates;
                    Switch closable = false;


                    for (int toli=0; toli<2; toli++)
                    {
                        scalar tol = toli*1e-4;

                        for(int pi=0; pi<2; pi++)
                        {
                            forAll(pointEdges[edges[edgej][pi]], edgeii)
                            {
                                label cand = pointEdges[edges[edgej][pi]][edgeii];
                                if (cand != edgej)
                                    if (isoV[cand] > -tol && isoV[cand] <= 1+tol)
                                        candidates.append(cand);
                                if (cand == edgei)
                                    closable = true;
                            }
                        }
                        if (closable || candidates.size() > 0)
                            break;
                    }

                    scalar minAngle = 2*M_PI;
                    label bestCand = -1;
                    vector bestN;

                    forAll(candidates, ci)
                    {
                        vector newN = edgePoint[candidates[ci]]-edgePoint[edgej];
                        scalar newD = mag(newN);

                        scalar angle = 0.0001;
                        if (newD != 0)
                        {
                            newN = newN/newD;

                            scalar angle1 = std::atan2(-n.y(), n.x());
                            scalar angle2 = std::atan2(-newN.y(), newN.x());

                            angle = angle2-angle1-M_PI;
                            while(angle > 2*M_PI)
                                angle -= 2*M_PI;
                            while (angle < 0)
                                angle += 2*M_PI;
                        }
                        else
                        {
                            newN = n;
                        }

                        if (angle < minAngle)
                        {
                            minAngle = angle;
                            bestN = newN;
                            bestCand = ci;
                        }
                    }
                    if (bestCand > -1)
                    {
                        edgej = candidates[bestCand];
                        p.append(edgePoint[edgej]);
                        isoV[edgej] = -1;
                        n = bestN;
                    }
                    else
                    {
                        if (!closable)
                        {
                            Info << "Incomplete isolone. Last edge " << edgej << endl;
                        }
                        break;
                    }

                }
                while(edgej != edgei);

                int rI = shp.addRecord();
                shp.addPart(rI);
                forAll(p, pi)
                {
                    shp.addPoint(rI, p[pi].x()-offset_.x(),
                                     p[pi].y()-offset_.y(),
                                     p[pi].z()-offset_.z());
                }
                shp.addPoint(rI, p[0].x()-offset_.x(),
                                 p[0].y()-offset_.y(),
                                 p[0].z()-offset_.z());

                isoValuesForRecords.append(isoValue);
            }
        }
    }

    if (!isoValuesForRecords.empty())
    {
        int fieldIndex = shp.addField(fieldName_, 'F', 12, 6);
        forAll(isoValuesForRecords, ri)
        {
            shp.setField(ri, fieldIndex, isoValuesForRecords[ri]);
        }
    }


    fileName fn
    (
        fileName::concat(obr_.time().timePath(), fileName_)
    );

    Info<< "Writing to " << fn << nl << nl;
    shp.calcBoundingBox();
    shp.write(fn);

    return true;
}


// ************************************************************************* //
