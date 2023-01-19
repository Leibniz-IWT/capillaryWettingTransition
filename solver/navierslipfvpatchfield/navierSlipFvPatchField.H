/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) held by original author
     \\/     M anipulation  |
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

Author:
    Dirk Gruending

\*---------------------------------------------------------------------------*/

#ifndef navierSlipFvPatchField_H
#define navierSlipFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class navierSlipFvPatch Declaration
\*---------------------------------------------------------------------------*/

class navierSlipFvPatchField
:
    public fvPatchField<vector>
{
        scalarField slipLength_;
        vectorField wallVelocity_;
public:

    TypeName("navierSlip");

        navierSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        navierSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        navierSlipFvPatchField
        (
            const navierSlipFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        navierSlipFvPatchField
        (
            const navierSlipFvPatchField&
        );

        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new navierSlipFvPatchField(*this)
            );
        }

        navierSlipFvPatchField
        (
            const navierSlipFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new navierSlipFvPatchField(*this, iF)
            );
        }
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

            virtual bool fixesValue() const
            {
                return true;
            }

            const scalarField& slipLength() const
            {
                return slipLength_;
            }

            void setSlipLength(tmp<scalarField> slipLength);

            const vectorField& wallVelocity() const
            {
                return wallVelocity_;
            }

            virtual tmp<Field<vector> > snGrad() const;

            virtual void evaluate
            (
                const UPstream::commsTypes commsType = UPstream::commsTypes::blocking
            );

            virtual tmp<Field<vector> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            virtual tmp<Field<vector> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            virtual tmp<Field<vector> > gradientInternalCoeffs() const;

            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        virtual void write(Ostream&) const;
        virtual void operator=(const UList<vector>&) {}
        virtual void operator=(const fvPatchField<vector>&) {}
        virtual void operator+=(const fvPatchField<vector>&) {}
        virtual void operator-=(const fvPatchField<vector>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}
        virtual void operator+=(const Field<vector>&) {}
        virtual void operator-=(const Field<vector>&) {}
        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}
        virtual void operator=(const vector&) {}
        virtual void operator+=(const vector&) {}
        virtual void operator-=(const vector&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
