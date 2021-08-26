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

#include "uniformVelAMIFvPatchField.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include <iostream>
#include <fstream>
using namespace std;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformVelAMIFvPatchField<Type>::uniformVelAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<Type>(p, iF),
    jumpTable_(),
    Ubarx_(0),
    Ubary_(0),
    Ubarz_(0),
    Utol_(0.0001),
    Relax_(1),
    jumpTable2_(1)
{}


template<class Type>
Foam::uniformVelAMIFvPatchField<Type>::uniformVelAMIFvPatchField
(
    const uniformVelAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<Type>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_.clone()),
    Ubarx_(ptf.Ubarx_),
    Ubary_(ptf.Ubary_),
    Ubarz_(ptf.Ubarz_),
    Utol_(ptf.Utol_),
    Relax_(ptf.Relax_),
    jumpTable2_(ptf.jumpTable2_)
{}


template<class Type>
Foam::uniformVelAMIFvPatchField<Type>::uniformVelAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<Type>(p, iF),
    jumpTable_(),
    Ubarx_(readScalar(dict.lookup("Ubarx"))),
    Ubary_(readScalar(dict.lookup("Ubary"))),
    Ubarz_(readScalar(dict.lookup("Ubarz"))),
    Utol_(dict.lookupOrDefault<scalar>("Utol",0.00001)),
    Relax_(dict.lookupOrDefault<scalar>("Relax",1)),
    UbarVec_(Ubarx_,Ubary_,Ubarz_),
    Ubar_(mag(UbarVec_)),
    flowDir_(UbarVec_/mag(UbarVec_)),
    jumpTable2_(dict.lookupOrDefault<scalar>("jumpTable2",1.0))
{

    if (this->cyclicAMIPatch().owner())
    {
        jumpTable_ = Function1<Type>::New("jumpTable", dict);

    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::uniformVelAMIFvPatchField<Type>::uniformVelAMIFvPatchField
(
    const uniformVelAMIFvPatchField<Type>& ptf
)
:
    fixedJumpAMIFvPatchField<Type>(ptf),
    jumpTable_(ptf.jumpTable_.clone()),
    Ubarx_(ptf.Ubarx_),
    Ubary_(ptf.Ubary_),
    Ubarz_(ptf.Ubarz_),
    Utol_(ptf.Utol_),
    Relax_(ptf.Relax_),
    jumpTable2_(ptf.jumpTable2_)
{}


template<class Type>
Foam::uniformVelAMIFvPatchField<Type>::uniformVelAMIFvPatchField
(
    const uniformVelAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<Type>(ptf, iF),
    jumpTable_(ptf.jumpTable_.clone()),
    Ubarx_(ptf.Ubarx_),
    Ubary_(ptf.Ubary_),
    Ubarz_(ptf.Ubarz_),
    Utol_(ptf.Utol_),
    Relax_(ptf.Relax_),
    jumpTable2_(ptf.jumpTable2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
Foam::scalar Foam::uniformVelAMIFvPatchField<Type>::magUbarAve
(
    const volVectorField& U
) const
{
    scalar magUbarAve = 0.0;
    // const fvMesh & BCpatch = this->patch().boundaryMesh().mesh();
    const scalarField& cv = this->patch().boundaryMesh().mesh().V();
    scalar sumVol = 0.0;

    // average of U
    forAll(cv, i)
    {
        sumVol += cv[i];
        scalar volCell = cv[i];
        magUbarAve += (flowDir_ & U[i])*volCell;
    }

    reduce(sumVol, sumOp<scalar>());
    reduce(magUbarAve, sumOp<scalar>());

    magUbarAve /= sumVol;

    return magUbarAve;
}

template<class Type>
void Foam::uniformVelAMIFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    if (this->cyclicAMIPatch().owner())
    {
        // const fvPatch & BCmesh = this->patch();
        // const scalarField & Surf = BCmesh.magSf();

        this->jump_ = jumpTable_->value(this->db().time().value());

        if(this->db().objectRegistry::foundObject<volVectorField>("U")){
            // if(a_){
            //     volVectorField & Uinit = this->db().objectRegistry::lookupObjectRef<volVectorField>("U");
            //     forAll(Uinit, cellI) {
            //        Uinit[cellI].x() = Ubarx_; 
            //        Uinit[cellI].y() = Ubary_; 
            //        Uinit[cellI].z() = Ubarz_; 
            //        a_ = false;
            //     }      
            // }

            const volVectorField & U = this->db().objectRegistry::lookupObject<volVectorField>("U");
            
            IOdictionary solDict
            (
            	IOobject
            	(
            	"fvSolution",
            	U.time().system(),
            	U.mesh(),
            	IOobject::MUST_READ,
            	IOobject::AUTO_WRITE
            	)
            );

            dictionary& solverControl = solDict.subDict("SIMPLE");
            int NOC_ = readScalar(solverControl.lookup("nNonOrthogonalCorrectors"))+1;
            i_ = (i_+1) % (NOC_);

            ofstream dp("system/Dp");
            dp << "Dp  "<< jumpTable2_<<";";
            dp.close();
            if(i_ == 0){

                scalar magU = this->magUbarAve(U);
                Info<<"velocity average in the flow dir :" <<magU<< nl;
                Info<<"Utol :" <<Utol_<< nl;
                float a = 1/Utol_;

                if(round((1   -  (magU  -  Ubar_)/Ubar_)*a)/a != 1 ){
                    jumpTable2_ *= (1   - Relax_ *  (magU  -  Ubar_)/Ubar_);
                    Info<<"The jump has been increased by  :" << 1   -  (magU  -  Ubar_)/Ubar_ << nl;
                    Info<<"Jump  :" << jumpTable2_ << nl;
                }else{
                    jumpTable2_ *=  round((1   - Relax_ *  (magU  -  Ubar_)/Ubar_)*a)/a;
                    Info<<"The Jump is now constant at  :" << jumpTable2_ << " with an error of :"<< (magU  -  Ubar_)<<" on the velocity" << nl;
                }

                // file << magU << " " << jumpTable2_ << std::endl;
                this->jump_ *= jumpTable2_;
                // if(jumpTable2_ < 0){
                //     Info<<"The new Jump is  :" << 0 << nl;
                //     this->jump_ = this->jump_*0;
                // }else if(jumpTable2_ > 1000){
                //     jumpTable2_ = 1000;
                //     this->jump_ = this->jump_*jumpTable2_;
                //     Info<<"The new Jump is  :" << jumpTable2_ << nl;
                // }
            
            }else{
                this->jump_ *= jumpTable2_;
            }

        }
    }
    fixedJumpAMIFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformVelAMIFvPatchField<Type>::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<Type>::write(os);
    if (this->cyclicAMIPatch().owner())
    {
        jumpTable_->writeData(os);
    }
    os.writeKeyword("Ubarx")<< Ubarx_ << token::END_STATEMENT<<nl ;
    os.writeKeyword("Ubary")<< Ubary_ << token::END_STATEMENT<<nl ;
    os.writeKeyword("Ubarz")<< Ubarz_ << token::END_STATEMENT<<nl ;
    os.writeKeyword("Utol")<< Utol_ << token::END_STATEMENT<<nl ;
    os.writeKeyword("Relax")<< Relax_ << token::END_STATEMENT<<nl ;
    os.writeKeyword("jumpTable2")<< jumpTable2_ << token::END_STATEMENT<<nl ;
}


// ************************************************************************* //
