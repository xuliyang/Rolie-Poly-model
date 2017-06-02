    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	
    dimensionedScalar nus
    (
        transportProperties.lookup("nus")
    );


    dimensionedScalar epsilon
    (
        transportProperties.lookup("epsilon")
    );

    dimensionedScalar nup
    (
        transportProperties.lookup("nup")
    );

    dimensionedScalar lambda_D
    (
        transportProperties.lookup("lambda_D")
    );


    dimensionedScalar lambda_R
    (
        transportProperties.lookup("lambda_R")
    );

    dimensionedScalar L_2 
    (
        transportProperties.lookup("L_2")
    );

    dimensionedTensor uT
    (
     	transportProperties.lookup("uT")

    );

    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );


    dimensionedScalar kappa
    (
        transportProperties.lookup("kappa")

    );

   dimensionedScalar rho
    (
        transportProperties.lookup("rho")

    );

   dimensionedScalar beta_CCR
    (
        transportProperties.lookup("beta_CCR")

    );


   dimensionedScalar xi
    (
        transportProperties.lookup("xi")
 
    );


    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


/*    Info<< "Reading field p2\n" << endl;
    volScalarField p2
    (
        IOobject
        (
            "p2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

*/
 /*   Info<< "Reading field U2\n" << endl;
    volVectorField U2
    (
        IOobject
        (
            "U2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

/*    Info<< "Reading field deltap\n" << endl;
    volScalarField deltap
    (
        IOobject
        (
            "deltap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

 /*   Info<< "Reading field deltaU\n" << endl;
    volVectorField deltaU
    (
        IOobject
        (
            "deltaU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

    Info<< "Reading field sigmap\n" << endl;
    volTensorField sigmap 
    (
        IOobject
        (
            "sigmap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field sigmas\n" << endl;
    volTensorField sigmas
    (
        IOobject
        (
            "sigmas",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("zero",dimensionSet(1,-1,-2,0,0,0,0),tensor::zero)
    );





/*    Info<< "Reading field sigmas\n" << endl;
    volTensorField sigmas
    (
        IOobject
        (
            "sigmas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

/*    Info<< "Reading field sigmas2\n" << endl;
    volTensorField sigmas2
    (
        IOobject
        (
            "sigmas2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/
 /*   Info<< "Reading field deltasigmas\n" << endl;
    volTensorField deltasigmas
    (
        IOobject
        (
            "deltasigmas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/


    Info<< "Reading field sigma\n" << endl;
    volTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field TensorA\n" << endl;
    volTensorField TensorA
    (
        IOobject
        (
            "TensorA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );




    Info<< "Reading field GradU\n" << endl;
    volTensorField GradU
    (
        IOobject
        (
            "GradU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*    Info<< "Reading field GradU2\n" << endl;
    volTensorField GradU2
    (
        IOobject
        (
            "GradU2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/
    Info<< "Reading field Rel\n" << endl;
    volScalarField Rel
    (
        IOobject
        (
            "Rel",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*    Info<< "Reading field traceA\n" << endl;
    volScalarField aTRACEA
    (
        IOobject
        (
            "aTRACEA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field traceA\n" << endl;
    volScalarField bTRACEA
    (
        IOobject
        (
            "bTRACEA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

    Info<< "Reading field traceA\n" << endl;
    volScalarField traceA
    (
        IOobject
        (
            "traceA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*    Info<< "Reading field traceA\n" << endl;
    volScalarField TRACEA
    (
        IOobject
        (
            "TRACEA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

		//Add by Liyang Xu
        Info<< "    Calculating vorticity" << endl;
        volVectorField vorticity
        (
            IOobject
            (
                "vorticity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::curl(U)
        );
/*
        Info<< "    Calculating vorticity2" << endl;
        volVectorField vorticity2
        (
            IOobject
            (
                "vorticity2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::curl(U2)
        );
*/

/*	Info<< "Reading/calculating face flux field phi2\n" << endl;

	surfaceScalarField phi2
	(
    	IOobject
   	 (
       		 "phi2",
       		runTime.timeName(),
        	mesh,
        	IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    	),
    	linearInterpolate(U) & mesh.Sf()
	);
    label pRefCell2 = 0;
    scalar pRefValue2 = 0.0;
*/
 //   setRefCell(p2, mesh.solutionDict().subDict("PISO"), pRefCell2, pRefValue2);

