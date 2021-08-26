/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::uniformVelAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclicAMI
    condition as a base.  The jump is specified as a time-varying uniform
    value across the patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclicAMI| yes |
        jumpTable    | jump value              | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformJumpAMI;
        patchType       cyclicAMI;
        jumpTable       constant 10;
    }
    \endverbatim

    The above example shows the use of a fixed jump of '10'.

Note
    The uniformValue entry is a Function1 type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

    The underlying \c patchType should be set to \c cyclic

See also
    Foam::jumpCyclicAMIFvPatchField
    Foam::Function1Types

SourceFiles
    uniformVelAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformVelAMIFvPatchField_H
#define uniformVelAMIFvPatchField_H

#include "fixedJumpAMIFvPatchField.H"
#include "Function1.H"
#include <iostream>
#include <fstream>
// #include "Jumpforces.H"
// #include "patchMeanVelocityForce.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class uniformVelAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformVelAMIFvPatchField
:
    public fixedJumpAMIFvPatchField<Type>
{
// friend class Jumpforces
protected:

    // Protected data

        //- "jump" table
        autoPtr<Function1<Type>> jumpTable_;
        
        //- Desired mean velocity
        scalar Ubarx_;
        scalar Ubary_;
        scalar Ubarz_;
        scalar Utol_;
        scalar Relax_;
        int i_ = 0;
        int NOC_ = 0;
        bool a_ = true;
        Vector<scalar> UbarVec_;
        scalar Ubar_;
        vector flowDir_;
        // autoPtr<OFstream> forceFilePtr_;
        scalar jumpTable2_;
        virtual scalar magUbarAve(const volVectorField& U) const;

        // const Type& Val_;

        // Foam::fv::patchMeanVelocityForce Mean_;

public:

    //- Runtime type information
    TypeName("uniformVelAMI");

    // Constructors

        //- Construct from patch and internal field
        uniformVelAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformVelAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformVelAMIFvPatchField onto a
        //  new patch
        uniformVelAMIFvPatchField
        (
            const uniformVelAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformVelAMIFvPatchField
        (
            const uniformVelAMIFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new uniformVelAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformVelAMIFvPatchField
        (
            const uniformVelAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new uniformVelAMIFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformVelAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
