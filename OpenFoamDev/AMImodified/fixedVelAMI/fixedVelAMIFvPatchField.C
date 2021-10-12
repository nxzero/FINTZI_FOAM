/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "fixedVelAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedVelAMIFvPatchField<Type>::fixedVelAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    VelCyclicAMIFvPatchField<Type>(p, iF),
    jump_(this->size(), Zero)
{}


template<class Type>
Foam::fixedVelAMIFvPatchField<Type>::fixedVelAMIFvPatchField
(
    const fixedVelAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    VelCyclicAMIFvPatchField<Type>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper)
{}


template<class Type>
Foam::fixedVelAMIFvPatchField<Type>::fixedVelAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    VelCyclicAMIFvPatchField<Type>(p, iF),
    jump_(p.size(), Zero)
{
    if (this->cyclicAMIPatch().owner())
    {
        jump_ = Field<Type>("jump", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );

    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::fixedVelAMIFvPatchField<Type>::fixedVelAMIFvPatchField
(
    const fixedVelAMIFvPatchField<Type>& ptf
)
:
    VelCyclicAMIFvPatchField<Type>(ptf),
    jump_(ptf.jump_)
{}


template<class Type>
Foam::fixedVelAMIFvPatchField<Type>::fixedVelAMIFvPatchField
(
    const fixedVelAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    VelCyclicAMIFvPatchField<Type>(ptf, iF),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fixedVelAMIFvPatchField<Type>::jump() const
{
    if (this->cyclicAMIPatch().owner())
    {
        // Info<< "Jump = " << jump_ <<"Owner"<< nl << endl;
        return jump_;
    }
    else
    {
        const fixedVelAMIFvPatchField& nbrPatch =
            refCast<const fixedVelAMIFvPatchField<Type>>
            (
                this->neighbourPatchField()
            );

        if (this->cyclicAMIPatch().applyLowWeightCorrection())
        {
            return this->cyclicAMIPatch().interpolate
            (
                nbrPatch.jump(),
                Field<Type>(this->size(), Zero)

            );
        }
        else
        {
            
            //Info<< "Jump = " << jump_ <<"Correction"<< nl << endl;
            return this->cyclicAMIPatch().interpolate(nbrPatch.jump());

        }
    }
}


template<class Type>
void Foam::fixedVelAMIFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    VelCyclicAMIFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
}


template<class Type>
void Foam::fixedVelAMIFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    VelCyclicAMIFvPatchField<Type>::rmap(ptf, addr);

    const fixedVelAMIFvPatchField<Type>& tiptf =
        refCast<const fixedVelAMIFvPatchField<Type>>(ptf);
    jump_.rmap(tiptf.jump_, addr);

}


template<class Type>
void Foam::fixedVelAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("patchType", this->interfaceFieldType());

    if (this->cyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
        Info<< "Jump = " << jump_ <<"Correction"<< nl << endl;

    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
