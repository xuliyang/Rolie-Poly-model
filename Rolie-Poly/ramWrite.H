/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Description


\*---------------------------------------------------------------------------*/

//store to ram
//static scalar timestamp = lambda1.value()/1000;

//if((timestamp - runTime.value()) < runTime.deltaTValue())
//{
//    do{
//	timestamp += lambda1.value()/1000;
//    }while((timestamp + lambda1.value()/1000 - runTime.value()) < runTime.deltaTValue());
    
    vorticity = fvc :: curl(U);
 //   vorticity2 = fvc :: curl(U2);

    for(int i = 0; i< ramWriteList_size; i++){
	(*(ramWriteList[i].ramP))[step] = p[ramWriteList[i].cellID];
       // (*(ramWriteList[i].ramP2))[step] = p2[ramWriteList[i].cellID];
//        (*(ramWriteList[i].ramdeltaP))[step] = deltap[ramWriteList[i].cellID];
	(*(ramWriteList[i].ramU))[step] = U[ramWriteList[i].cellID];
//        (*(ramWriteList[i].ramU2))[step] = U2[ramWriteList[i].cellID];
 //       (*(ramWriteList[i].ramdeltaU))[step] = deltaU[ramWriteList[i].cellID];
	(*(ramWriteList[i].ramVorticity))[step] = vorticity[ramWriteList[i].cellID];
  //      (*(ramWriteList[i].ramVorticity2))[step] = vorticity2[ramWriteList[i].cellID];
	(*(ramWriteList[i].ramTensorA))[step] = TensorA[ramWriteList[i].cellID];
	(*(ramWriteList[i].ramSigmaP))[step] = sigmap[ramWriteList[i].cellID];
	(*(ramWriteList[i].ramSigma))[step] = sigma[ramWriteList[i].cellID];
        (*(ramWriteList[i].ramSigmas))[step] = sigmas[ramWriteList[i].cellID];
    //    (*(ramWriteList[i].ramSigmas2))[step] = sigmas2[ramWriteList[i].cellID];
//        (*(ramWriteList[i].ramdeltaSigmas))[step] = deltasigmas[ramWriteList[i].cellID];
        (*(ramWriteList[i].ramGradU))[step] = GradU[ramWriteList[i].cellID];
      //  (*(ramWriteList[i].ramGradU2))[step] = GradU2[ramWriteList[i].cellID];
	//Pout << "write P = "<<(*(ramWriteList[i].ramP))[step]<<" to cell "<< i <<" at step "<<step<<endl;
    }
    step ++;

    //ofstream to file
    if(step == outputStep){
	
	for(int i = 0; i< ramWriteList_size; i++){
		//Pout <<"filep: "<<(*(ramWriteList[i].fileP))<<" **************\n";
		*(ramWriteList[i].fileP) << SubList<scalar>(*(ramWriteList[i].ramP),step)<<endl;
		ramWriteList[i].fileP->flush();

 //               *(ramWriteList[i].fileP2) << SubList<scalar>(*(ramWriteList[i].ramP2),step)<<endl;
//                ramWriteList[i].fileP2->flush();

//                *(ramWriteList[i].filedeltaP) << SubList<scalar>(*(ramWriteList[i].ramdeltaP),step)<<endl;
//                ramWriteList[i].filedeltaP->flush();
		
		*(ramWriteList[i].fileU) << *(ramWriteList[i].ramU)<<endl;
		ramWriteList[i].fileU->flush();

//                *(ramWriteList[i].fileU2) << *(ramWriteList[i].ramU2)<<endl;
//                ramWriteList[i].fileU2->flush();

   //             *(ramWriteList[i].filedeltaU) << *(ramWriteList[i].ramdeltaU)<<endl;
     //           ramWriteList[i].filedeltaU->flush();

		*(ramWriteList[i].fileVorticity) << *(ramWriteList[i].ramVorticity)<<endl;
		ramWriteList[i].fileVorticity->flush();
		
//                *(ramWriteList[i].fileVorticity2) << *(ramWriteList[i].ramVorticity2)<<endl;
  //              ramWriteList[i].fileVorticity2->flush();

		*(ramWriteList[i].fileTensorA) << *(ramWriteList[i].ramTensorA)<<endl;
		ramWriteList[i].fileTensorA->flush();

		*(ramWriteList[i].fileSigmaP) << *(ramWriteList[i].ramSigmaP)<<endl;
		ramWriteList[i].fileSigmaP->flush();

		
		*(ramWriteList[i].fileSigma) << *(ramWriteList[i].ramSigma)<<endl;
		ramWriteList[i].fileSigma->flush();


                *(ramWriteList[i].fileSigmas) << *(ramWriteList[i].ramSigmas)<<endl;
                ramWriteList[i].fileSigmas->flush();


 //               *(ramWriteList[i].fileSigmas2) << *(ramWriteList[i].ramSigmas2)<<endl;
//                ramWriteList[i].fileSigmas2->flush();


 //               *(ramWriteList[i].filedeltaSigmas) << *(ramWriteList[i].ramdeltaSigmas)<<endl;
   //             ramWriteList[i].filedeltaSigmas->flush();


                *(ramWriteList[i].fileGradU) << *(ramWriteList[i].ramGradU)<<endl;
                ramWriteList[i].fileGradU->flush();

 //               *(ramWriteList[i].fileGradU2) << *(ramWriteList[i].ramGradU2)<<endl;
//                ramWriteList[i].fileGradU2->flush();
	}
	step = 0;
    }
//}
