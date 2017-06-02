//ofstream to file
if(step >0){
	
	for(int i = 0; i< ramWriteList_size; i++){
		//Pout <<"filep: "<<(*(ramWriteList[i].fileP))<<" **************\n";
		(ramWriteList[i].ramP)->setSize(step);
		*(ramWriteList[i].fileP) << *(ramWriteList[i].ramP)<<endl;
		ramWriteList[i].fileP->flush();

                (ramWriteList[i].ramP2)->setSize(step);
                *(ramWriteList[i].fileP2) << *(ramWriteList[i].ramP2)<<endl;
                ramWriteList[i].fileP2->flush();

                (ramWriteList[i].ramdeltaP)->setSize(step);
                *(ramWriteList[i].filedeltaP) << *(ramWriteList[i].ramdeltaP)<<endl;
                ramWriteList[i].filedeltaP->flush();

		(ramWriteList[i].ramU)->setSize(step);
		*(ramWriteList[i].fileU) << *(ramWriteList[i].ramU)<<endl;
		ramWriteList[i].fileU->flush();

                (ramWriteList[i].ramU2)->setSize(step);
                *(ramWriteList[i].fileU2) << *(ramWriteList[i].ramU2)<<endl;
                ramWriteList[i].fileU2->flush();

                (ramWriteList[i].ramdeltaU)->setSize(step);
                *(ramWriteList[i].filedeltaU) << *(ramWriteList[i].ramdeltaU)<<endl;
                ramWriteList[i].filedeltaU->flush();

		(ramWriteList[i].ramVorticity)->setSize(step);
		*(ramWriteList[i].fileVorticity) << *(ramWriteList[i].ramVorticity)<<endl;
		ramWriteList[i].fileVorticity->flush();

                (ramWriteList[i].ramVorticity2)->setSize(step);
                *(ramWriteList[i].fileVorticity2) << *(ramWriteList[i].ramVorticity2)<<endl;
                ramWriteList[i].fileVorticity2->flush();

		(ramWriteList[i].ramTensorA)->setSize(step);
		*(ramWriteList[i].fileTensorA) << *(ramWriteList[i].ramTensorA)<<endl;
		ramWriteList[i].fileTensorA->flush();

		(ramWriteList[i].ramSigmaP)->setSize(step);
		*(ramWriteList[i].fileSigmaP) << *(ramWriteList[i].ramSigmaP)<<endl;
		ramWriteList[i].fileSigmaP->flush();


		(ramWriteList[i].ramSigma)->setSize(step);
		*(ramWriteList[i].fileSigma) << *(ramWriteList[i].ramSigma)<<endl;
		ramWriteList[i].fileSigma->flush();

                (ramWriteList[i].ramSigmas)->setSize(step);
                *(ramWriteList[i].fileSigmas) << *(ramWriteList[i].ramSigmas)<<endl;
                ramWriteList[i].fileSigmas->flush();

                (ramWriteList[i].ramSigmas2)->setSize(step);
                *(ramWriteList[i].fileSigmas2) << *(ramWriteList[i].ramSigmas2)<<endl;
                ramWriteList[i].fileSigmas2->flush();

                (ramWriteList[i].ramdeltaSigmas)->setSize(step);
                *(ramWriteList[i].filedeltaSigmas) << *(ramWriteList[i].ramdeltaSigmas)<<endl;
                ramWriteList[i].filedeltaSigmas->flush();

                (ramWriteList[i].ramGradU)->setSize(step);
                *(ramWriteList[i].fileGradU) << *(ramWriteList[i].ramGradU)<<endl;
                ramWriteList[i].fileGradU->flush();

                (ramWriteList[i].ramGradU2)->setSize(step);
                *(ramWriteList[i].fileGradU2) << *(ramWriteList[i].ramGradU2)<<endl;
                ramWriteList[i].fileGradU2->flush();

	}
	step = 0;
}

//clean the newed data
for(int i = 0; i< ramWriteList_size; i++){
	delete ramWriteList[i].fileP;
	delete ramWriteList[i].ramP;

        delete ramWriteList[i].fileP2;
        delete ramWriteList[i].ramP2;

        delete ramWriteList[i].filedeltaP;
        delete ramWriteList[i].ramdeltaP;

	delete ramWriteList[i].fileU;
	delete ramWriteList[i].ramU;

        delete ramWriteList[i].fileU2;
        delete ramWriteList[i].ramU2;

        delete ramWriteList[i].filedeltaU;
        delete ramWriteList[i].ramdeltaU;

	delete ramWriteList[i].fileVorticity;
	delete ramWriteList[i].ramVorticity;


        delete ramWriteList[i].fileVorticity2;
        delete ramWriteList[i].ramVorticity2;

	delete ramWriteList[i].fileTensorA;
	delete ramWriteList[i].ramTensorA;

	delete ramWriteList[i].fileSigmaP;
	delete ramWriteList[i].ramSigmaP;


	delete ramWriteList[i].fileSigma;
	delete ramWriteList[i].ramSigma;

        delete ramWriteList[i].fileSigmas;
        delete ramWriteList[i].ramSigmas;

        delete ramWriteList[i].fileSigmas2;
        delete ramWriteList[i].ramSigmas2;

        delete ramWriteList[i].filedeltaSigmas;
        delete ramWriteList[i].ramdeltaSigmas;

        delete ramWriteList[i].fileGradU;
        delete ramWriteList[i].ramGradU;

        delete ramWriteList[i].fileGradU2;
        delete ramWriteList[i].ramGradU2;
}

