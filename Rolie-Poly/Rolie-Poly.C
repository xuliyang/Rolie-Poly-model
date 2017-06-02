/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian 
    fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

#	include "ramWriteInit.H"

    //mesh.clearPrimitives();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

             
        GradU = fvc::grad(U);
        volTensorField gradU = fvc::grad(U);

      for (int corr = 0; corr < nCorr; corr++)
      {
//define LinearPTT equation
 /*      tmp<fvTensorMatrix> TEqn
        (
	 fvm::ddt(sigmap)
        + fvm::div(phi,sigmap)
	==
	 nup/lambda1*(gradU + gradU.T())
        + (sigmap & gradU) +  (gradU.T() & sigmap)
	- (fvm::Sp(epsilon/nup*tr(sigmap), sigmap) + fvm::Sp(1/lambda1,sigmap))
        - (xi/1)*( (gradU & sigmap) + (sigmap & gradU.T()) + (gradU.T() & sigmap) + (sigmap & gradU)) 
        );
*/ 
 //solve LinearPTT equation      
       //  TEqn().relax();
//	 solve(TEqn);
//	gradU.write();
//	(gradU.T()).write();

//define FENE-CD-JS equation
/*       tmp<fvTensorMatrix> TEqn
        (
        fvm::ddt(TensorA)
        + fvm::div(phi,TensorA)
       ==
       (TensorA & gradU) + (gradU.T() & TensorA)
	-L_2/((L_2 - tr(TensorA))*lambda1*sqr(1-k+k*sqrt(traceA/2))) * (fvm::Sp(1,TensorA) - uT)
	-xi/2*( (gradU & TensorA) + (TensorA & gradU.T()) + (gradU.T() & TensorA) + (TensorA & gradU))
        );
//solve FENE-CD-JS equation to obtain TensorA
	solve(TEqn);
        traceA = tr(TensorA);
*/
//define Rolie-Poly equation
	traceA = tr(TensorA);

	tmp<fvTensorMatrix> TEqn
        (
        fvm::ddt(TensorA)
        + fvm::div(phi,TensorA)
        ==
        (TensorA & gradU) + (gradU.T() & TensorA)
	-(TensorA-uT)/lambda_D
	-2*(1-sqrt(2/traceA))/lambda_R*(TensorA+beta_CCR*sqrt(traceA/2)*(TensorA-uT))
        );
	
	solve(TEqn);
	
	sigmap = nup/lambda_D * TensorA; //FENE-CD-JS
	sigmas = nus* (gradU + gradU.T());
        // calculate sigma (extra stress tensor)
        sigma = sigmap + sigmas;
	traceA = tr(TensorA);

  //      solve(TEqn == nup/lambda1*(gradU + gradU.T()) );

//      for (int corr = 0; corr < nCorr; corr++)
 //     {

//define momentum equation
            tmp<fvVectorMatrix> UEqn
            (
                fvm::ddt(U)
	       +fvm::div(phi, U)
               +fvc::laplacian(kappa/rho, U)
	       -fvm::laplacian((nus+kappa)/rho, U)
               -fvc::div(sigmap/rho)
            );
//solve for U
            UEqn().relax();
            Info << "END U.Eqn relax=" << nl << endl;
	    solve(UEqn() == -fvc::grad(p));

	    p.boundaryField().updateCoeffs();
            volScalarField AU = UEqn().A();
            U = UEqn().H()/AU;
            UEqn.clear();
            phi = fvc::interpolate(U) & mesh.Sf();
            adjustPhi(phi, U, p);

	    Info << "END solveU=" << nl << endl;
            p.storePrevIter();

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/AU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U -= fvc::grad(p)/AU;
            U.correctBoundaryConditions();
       } 

	#include "ramWrite.H"
		
        runTime.write();

        Info<< "ExecutionTime+++++++++ = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    #include "ramWriteClean.H"
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
