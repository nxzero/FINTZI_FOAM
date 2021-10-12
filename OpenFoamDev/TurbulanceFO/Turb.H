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

Class
    Foam::functionObjects::Turb

Group
    grpFieldFunctionObjects

Description
    Calculates the energy spectrum for a structured IJK mesh

Usage
    Example of function object specification:
    \verbatim
    Turb1
    {
        type        Turb;
        libs        (fieldFunctionObjects);
        ...
        log         yes;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description               | Required    | Default value
        type         | type name: Turb | yes         |
        log          | write info to standard output | no      | yes
    \endtable

    Output data is written to the file \<timeDir\>/Turb.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    Turb.C
    TurbTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Turb_H
#define functionObjects_Turb_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"
#include <iostream>
#include <fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class Turb Declaration
\*---------------------------------------------------------------------------*/

class Turb
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Name of velocity field, default = U
        word UName_;
        vector Uavg_;
        // tensor
        tensor T_;
        std::ofstream file_;
        autoPtr<OFstream> osPtr_ = createFile(name(), time_.value());

    // Protected Member Functions
        void magUbarAve(const vectorField& U);
        void Tensorcalc(const vectorField& U);

        void WriteTensor();
        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- No copy construct
        Turb(const Turb&) = delete;

        //- No copy assignment
        void operator=(const Turb&) = delete;


public:

    //- Runtime type information
    TypeName("Turb");


    // Constructors

        //- Construct from Time and dictionary
        Turb
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Turb() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the Turb
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
