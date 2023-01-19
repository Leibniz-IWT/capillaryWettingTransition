
    #include <cmath>


    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& deltaN(1e-16/pow(average(alpha1.mesh().V()), 1.0/3.0));

    //const dimensionedScalar& deltaN(1e-8);
    //const fvMesh& mesh = alpha1.mesh();
    //const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1));

    // Interpolated face-gradient of alpha
    //surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
    // Face unit interface normal
    //surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));
    volVectorField nHatfv(gradAlpha/(mag(gradAlpha) + deltaN));

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    dimensionedScalar nu1
    (
        transportProperties.subDict("water").lookup("nu")
    );

    dimensionedScalar nu2
    (
        transportProperties.subDict("air").lookup("nu")
    );

    // Need to store rho for ddt(rho, U)
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1*nu1 + alpha2*rho2*nu2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


volTensorField gradU = fvc::grad(U);   // is not written, and needs not to be declared in createFields.H if declaration (volTensorField) is done here

int a0;
int a1;
int a2;

//int b0;
//int b1;
//int b2;

int i = 0;
int j = 0;

scalar kp0;
scalar kp1;
scalar kp2;

scalar kn0;
scalar kn1;
scalar kn2;


const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

//--------------------------------------------------------------------------------------------------------------------
   D = (gradU+gradU.T());

	forAll(mesh.C(), celli)
	{

		if (nHatfv[celli].component(0) < 0)
		{
			a0 = -1;
		}
		else
		{
			a0 = 1;
		}

		if (nHatfv[celli].component(1) < 0)
		{
			a1 = -1;
		}
		else
		{
			a1 = 1;
		}

		if (nHatfv[celli].component(2) < 0)
		{
			a2 = -1;
		}
		else
		{
			a2 = 1;
		}



		// getting cell size
    		const cell & cc = mesh.cells()[celli];
    		labelList pLabels(cc.labels(ff));
    		pointField pLocal(pLabels.size(), vector::zero);

    		forAll (pLabels, pointi)
          		pLocal[pointi] = pp[pLabels[pointi]];

    		scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));

    		scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));

    		scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));

		label labelOfClosestCellP = celli;
		label labelOfClosestCellN = celli;
		label labelNOld = celli;
		label labelPOld = celli;

		//getting cell coordinates + shift in direction of normal vector
		if (alpha1[celli] < 0.9 && alpha1[celli] > 0.1) // 0.9 0.1
			{
			while (alpha1[labelOfClosestCellP] <= 0.99 && labelOfClosestCellP != celli) // 0.999
			{

				kp0 = mesh.C()[celli].component(0) + (i+1)*xDim*nHatfv[celli].component(0);
				kp1 = mesh.C()[celli].component(1) + (i+1)*yDim*nHatfv[celli].component(1);
				kp2 = mesh.C()[celli].component(2) + (i+1)*zDim*nHatfv[celli].component(2);

				labelPOld = labelOfClosestCellP;
				labelOfClosestCellP = mesh.findNearestCell(point(kp0, kp1, kp2));

				i++;

				if (alpha1[labelOfClosestCellP] < alpha1[labelPOld] && labelOfClosestCellP != celli) break;
				if (abs(alpha1[labelOfClosestCellP] - alpha1[labelPOld]) < 0.000001 && labelOfClosestCellP != celli) break;
			}

			i = 0;

			while (alpha1[labelOfClosestCellN] >= 0.01 && labelOfClosestCellN != celli) // 0.001
			{

				kn0 = mesh.C()[celli].component(0) - (j+1)*xDim*nHatfv[celli].component(0);

				kn1 = mesh.C()[celli].component(1) - (j+1)*yDim*nHatfv[celli].component(1);

				kn2 = mesh.C()[celli].component(2) - (j+1)*zDim*nHatfv[celli].component(2);
				labelNOld = labelOfClosestCellN;
				labelOfClosestCellN = mesh.findNearestCell(point(kn0, kn1, kn2));

				j++;

				if (alpha1[labelOfClosestCellN] > alpha1[labelNOld] && labelOfClosestCellN != celli) break;
				if (abs(alpha1[labelOfClosestCellN] - alpha1[labelNOld]) < 0.000001 && labelOfClosestCellN != celli) break;

			}

			j = 0;

				scalar dist = (pow((mesh.C()[labelOfClosestCellP].component(0)-mesh.C()[labelOfClosestCellN].component(0)),2) + pow((mesh.C()[labelOfClosestCellP].component(1)-mesh.C()[labelOfClosestCellN].component(1)),2) + pow((mesh.C()[labelOfClosestCellP].component(2)-mesh.C()[labelOfClosestCellN].component(2)),2));

				scalar dist_sqrt = Foam::pow(dist,0.5);

				scalar distCell = (pow((mesh.C()[labelOfClosestCellN].component(0)-mesh.C()[celli].component(0)),2) + pow((mesh.C()[labelOfClosestCellN].component(1)-mesh.C()[celli].component(1)),2) + pow((mesh.C()[labelOfClosestCellN].component(2)-mesh.C()[celli].component(2)),2));

				scalar distCell_sqrt = Foam::pow(distCell,0.5);

				//shearRate[celli]= (nHatfv[celli] & D[labelOfClosestCellN]) + (((nHatfv[celli] & D[labelOfClosestCellP]) - (nHatfv[celli] & D[labelOfClosestCellN]))/dist_sqrt)*distCell_sqrt;


			if (dist_sqrt == 0)
			{
				// shear
				shearRate[celli] = nHatfv[celli] & D[celli];

				// strain
				strainRate[celli].component(0) = (a0*D[celli].component(0)*(1-abs(nHatfv[celli].component(0))) + a1*D[celli].component(1)*(1-abs(nHatfv[celli].component(1))) + a2*D[celli].component(2)*(1-abs(nHatfv[celli].component(2))));
				strainRate[celli].component(1) = (a0*D[celli].component(3)*(1-abs(nHatfv[celli].component(0))) + a1*D[celli].component(4)*(1-abs(nHatfv[celli].component(1))) + a2*D[celli].component(5)*(1-abs(nHatfv[celli].component(2))));
				strainRate[celli].component(2) = (a0*D[celli].component(6)*(1-abs(nHatfv[celli].component(0))) + a1*D[celli].component(7)*(1-abs(nHatfv[celli].component(1))) + a2*D[celli].component(8)*(1-abs(nHatfv[celli].component(2))));
			}
			else
			{

				// shear rate at interface

				// x-component
				shearRate[celli].component(0) = (D[labelOfClosestCellN].component(0)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(1)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(2)*nHatfv[celli].component(2)) + (((D[labelOfClosestCellP].component(0)*nHatfv[celli].component(0) + D[labelOfClosestCellP].component(1)*nHatfv[celli].component(1) + D[labelOfClosestCellP].component(2)*nHatfv[celli].component(2)) - (D[labelOfClosestCellN].component(0)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(1)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(2)*nHatfv[celli].component(2)))/dist_sqrt)*distCell_sqrt;

				// y-component
				shearRate[celli].component(1) = (D[labelOfClosestCellN].component(3)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(4)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(5)*nHatfv[celli].component(2)) + (((D[labelOfClosestCellP].component(3)*nHatfv[celli].component(0) + D[labelOfClosestCellP].component(4)*nHatfv[celli].component(1) + D[labelOfClosestCellP].component(5)*nHatfv[celli].component(2)) - (D[labelOfClosestCellN].component(3)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(4)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(5)*nHatfv[celli].component(2)))/dist_sqrt)*distCell_sqrt;

				// z-component
				shearRate[celli].component(2) = (D[labelOfClosestCellN].component(6)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(7)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(8)*nHatfv[celli].component(2)) + (((D[labelOfClosestCellP].component(6)*nHatfv[celli].component(0) + D[labelOfClosestCellP].component(7)*nHatfv[celli].component(1) + D[labelOfClosestCellP].component(8)*nHatfv[celli].component(2)) - (D[labelOfClosestCellN].component(6)*nHatfv[celli].component(0) + D[labelOfClosestCellN].component(7)*nHatfv[celli].component(1) + D[labelOfClosestCellN].component(8)*nHatfv[celli].component(2)))/dist_sqrt)*distCell_sqrt;



				// strain rate at interface

				// x-component
				strainRate[celli].component(0) = (a0*D[labelOfClosestCellN].component(0)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(1)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(2)*(1-abs(nHatfv[celli].component(2)))) + (((a0*D[labelOfClosestCellP].component(0)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellP].component(1)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellP].component(2)*(1-abs(nHatfv[celli].component(2)))) - (a0*D[labelOfClosestCellN].component(0)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(1)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(2)*(1-abs(nHatfv[celli].component(2)))))/(alpha1[labelOfClosestCellP]-alpha1[labelOfClosestCellN]))*alpha1[celli];

				// y-component
				strainRate[celli].component(1) = (a0*D[labelOfClosestCellN].component(3)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(4)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(5)*(1-abs(nHatfv[celli].component(2)))) + (((a0*D[labelOfClosestCellP].component(3)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellP].component(4)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellP].component(5)*(1-abs(nHatfv[celli].component(2)))) - (a0*D[labelOfClosestCellN].component(3)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(4)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(5)*(1-abs(nHatfv[celli].component(2)))))/(alpha1[labelOfClosestCellP]-alpha1[labelOfClosestCellN]))*alpha1[celli];

				// z-component
				strainRate[celli].component(2) = (a0*D[labelOfClosestCellN].component(6)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(7)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(8)*(1-abs(nHatfv[celli].component(2)))) + (((a0*D[labelOfClosestCellP].component(6)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellP].component(7)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellP].component(8)*(1-abs(nHatfv[celli].component(2)))) - (a0*D[labelOfClosestCellN].component(6)*(1-abs(nHatfv[celli].component(0))) + a1*D[labelOfClosestCellN].component(7)*(1-abs(nHatfv[celli].component(1))) + a2*D[labelOfClosestCellN].component(8)*(1-abs(nHatfv[celli].component(2)))))/(alpha1[labelOfClosestCellP]-alpha1[labelOfClosestCellN]))*alpha1[celli];
			}


				// TauShear at interface
				tauShear[celli] = mu[celli]*shearRate[celli];

				// TauStrain at interface
				tauStrain[celli] = mu[celli]*strainRate[celli];

			}
		else
			{
		shearRate[celli].component(0) = 0;
		shearRate[celli].component(1) = 0;
		shearRate[celli].component(2) = 0;

		strainRate[celli].component(0) = 0;
		strainRate[celli].component(1) = 0;
		strainRate[celli].component(2) = 0;

		tauShear[celli].component(0) = 0;
		tauShear[celli].component(1) = 0;
		tauShear[celli].component(2) = 0;

		tauStrain[celli].component(0) = 0;
		tauStrain[celli].component(1) = 0;
		tauStrain[celli].component(2) = 0;
			}


	}

/*
//get interface shear stress at wall

#include "fvPatch.H"


	const fvPatchList& patches = mesh.boundary();

	surfaceTensorField gradUs(fvc::interpolate(gradU));
	Ds = (gradUs + gradUs.T());

	forAll(patches,patchi)
	{
				if (alpha1[patchi] < 0.5)
		{
 			const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
        		const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

			// wall shear rate
			vectorField ssp = (-Sfp/magSfp) & Ds;


			//const surfaceVectorField& wallShearRate1 = ssp;
			// x-component
			wallShearRate1[patchi].component(0) = ssp[patchi].component(0);
			// y-component
			wallShearRate1[patchi].component(1) = ssp[patchi].component(1);
			// z-component
			wallShearRate1[patchi].component(2) = ssp[patchi].component(2);




			// wall shear stress
			tauWallShear[patchi] = mu[patchi]*wallShearRate[patchi];


			// wall shear strain
			vectorField nSurf = -Sfp/magSfp;

		if (nSurf[patchi].component(0) < 0)
		{
			b0 = -1;
		}
		else
		{
			b0 = 1;
		}

		if (nSurf[patchi].component(1) < 0)
		{
			b1 = -1;
		}
		else
		{
			b1 = 1;
		}

		if (nSurf[patchi].component(2) < 0)
		{
			b2 = -1;
		}
		else
		{
			b2 = 1;
		}


			wallStrainRate[patchi].component(0) = b0*((vector(1,1,1).component(0) * Ds[patchi].component(0)) - abs(wallShearRate[patchi].component(0)));
			wallStrainRate[patchi].component(1) = b1*((vector(1,1,1).component(1) * Ds[patchi].component(1)) - abs(wallShearRate[patchi].component(1)));
			wallStrainRate[patchi].component(2) = b2*((vector(1,1,1).component(2) * Ds[patchi].component(2)) - abs(wallShearRate[patchi].component(2)));



			tauWallStrain = mu[patchi]*wallStrainRate[patchi];

		}
		else
		{
				wallShearRate[patchi].component(0) = 0;
				wallShearRate[patchi].component(1) = 0;
				wallShearRate[patchi].component(2) = 0;

				tauWallShear[patchi].component(0) = 0;
				tauWallShear[patchi].component(1) = 0;
				tauWallShear[patchi].component(2) = 0;

				wallStrainRate[patchi].component(0) = 0;
				wallStrainRate[patchi].component(1) = 0;
				wallStrainRate[patchi].component(2) = 0;

				tauWallStrain[patchi].component(0) = 0;
				tauWallStrain[patchi].component(1) = 0;
				tauWallStrain[patchi].component(2) = 0;
		}



	}
*/
