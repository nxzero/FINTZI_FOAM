/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "Turb.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "cartesianCS.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Turb, 0);
    addToRunTimeSelectionTable(functionObject, Turb, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::Turb::magUbarAve(const vectorField& U) 
{
    vector magUbarAve(0,0,0);
    const scalarField& cv = mesh_.V();
    scalar sumVol = 0.0;

    // average of U
    forAll(cv, i)
    {
        sumVol += cv[i];
        scalar volCell = cv[i];
        magUbarAve += U[i]*volCell;
    }

    reduce(sumVol, sumOp<scalar>());
    reduce(magUbarAve, sumOp<vector>());

    magUbarAve /= sumVol;
    Uavg_ = magUbarAve;
}

void Foam::functionObjects::Turb::Tensorcalc(const vectorField& U) 
{
    tensor  Ttmp(Zero);
    vector Udiff;
    const scalarField& cv = mesh_.V();
    scalar sumVol = 0.0;

    // average of U
    forAll(cv, i)
    {
        sumVol += cv[i];
        scalar volCell = cv[i];
        Udiff =  U[i] - Uavg_;
        Ttmp += Udiff * Udiff * volCell;
    }

    reduce(sumVol, sumOp<scalar>());
    reduce(Ttmp, sumOp<tensor>());

    Ttmp /= sumVol;
    T_ = Ttmp;
}


void Foam::functionObjects::Turb::writeFileHeader(Ostream& os)
{
    if (Pstream::master())
    {
    writeHeader(os, "Pseudo-Tubulent tensor = <Ui' Uj'>");

    // writeCommented(os,"time");

    os  << tab << 
     "(xx yy zz)" << tab <<
       "(xy xz zy)"<<nl;
    }
}

void Foam::functionObjects::Turb::WriteTensor()
{
    Log << "Writing Pseudo-Turbulent Tensor" << endl;
    // autoPtr<OFstream> osPtr = createFile(name(), time_.value());

    OFstream& os = osPtr_.ref();

    writeCurrentTime(os);
    
    os << T_.xx() << tab
    << T_.yy() << tab
    << T_.zz() << tab
    << T_.xy() << tab 
    << T_.xz() << tab
    << T_.zy() << tab;
    
    os << nl;


    Log << "Pseudo-turbulence"<< tab <<"T = "
    << T_.xx() << tab
    << T_.yy() << tab
    << T_.zz() << tab
    << T_.xy() << tab 
    << T_.xz() << tab
    << T_.zy() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Turb::Turb
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    UName_("U")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Turb::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    dict.readIfPresent("U", UName_);    
    OFstream& os = osPtr_.ref();
    writeFileHeader(os);
    // this->file_("postProcessing/Pseudo_turblulent_tensor.dat");
    return true;
}


bool Foam::functionObjects::Turb::execute()
{   
    return true;
}


bool Foam::functionObjects::Turb::write()
{
    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    magUbarAve(U);    
    Tensorcalc(U);
    if (Pstream::master())
    {
    WriteTensor();
    }
    return true;
}


// ************************************************************************* //
