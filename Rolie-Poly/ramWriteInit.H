int nCell =35;
int outputStep = 1000;
int step = 0;;

struct ramWriteCell{
	label cellID;
	label name;
	point location;
	
	OFstream* fileP;
        OFstream* fileP2;
        OFstream* filedeltaP;
	OFstream* fileU;
        OFstream* fileU2;
        OFstream* filedeltaU;
	OFstream* fileVorticity;
        OFstream* fileVorticity2;
	OFstream* fileTensorA;
	OFstream* fileSigmaP;
	OFstream* fileSigma;
        OFstream* fileSigmas;
        OFstream* fileSigmas2;
        OFstream* filedeltaSigmas;
        OFstream* fileGradU;
        OFstream* fileGradU2;
	//OFstream* fileMaxTrA;

	List<scalar>*  ramP;
        List<scalar>*  ramP2;
        List<scalar>*  ramdeltaP;
	List<vector>*  ramU;
        List<vector>*  ramU2;
        List<vector>*  ramdeltaU;
	List<vector>*  ramVorticity;
        List<vector>*  ramVorticity2;
	List<tensor>*  ramTensorA;
	List<tensor>*  ramSigmaP;
	List<tensor>*  ramSigma;
        List<tensor>*  ramSigmas;
        List<tensor>*  ramSigmas2;
        List<tensor>*  ramdeltaSigmas;
        List<tensor>*  ramGradU;
        List<tensor>*  ramGradU2;
	//List<scalar>*  ramMaxTrA;
};

ramWriteCell ramWriteList[nCell];
int ramWriteList_size = 0;


//16:1
//get cell point list
/*
List<point> ramWritePoints(nCell);
ramWritePoints[0]=point(0.303825, 0.0158429, 2e-05);
ramWritePoints[1]=point(0.303825, 0.0118527, 2e-05);
ramWritePoints[2]=point(0.303825, 0.00801027, 2e-05);
ramWritePoints[3]=point(0.303825, 0.00438403, 2e-05);
ramWritePoints[4]=point(0.303825, 0.000286184, 2e-05);
ramWritePoints[5]=point(0.312589, 0.0158429, 2e-05);
ramWritePoints[6]=point(0.312589, 0.0118527, 2e-05);
ramWritePoints[7]=point(0.312589, 0.00801027, 2e-05);
ramWritePoints[8]=point(0.312589, 0.00438403, 2e-05);
ramWritePoints[9]=point(0.312589, 0.000286184, 2e-05);
ramWritePoints[10]=point(0.319005, 0.0158429, 2e-05);
ramWritePoints[11]=point(0.319005, 0.0118527, 2e-05);
ramWritePoints[12]=point(0.319005, 0.00801027, 2e-05);
ramWritePoints[13]=point(0.319005, 0.00438403, 2e-05);
ramWritePoints[14]=point(0.319005, 0.000286184, 2e-05);
ramWritePoints[15]=point(0.320004, 0.00751752, 2e-05);
ramWritePoints[16]=point(0.320004, 0.0076874, 2e-05);
ramWritePoints[17]=point(0.320004, 0.00801027, 2e-05);
ramWritePoints[18]=point(0.320004, 0.00832447, 2e-05);
ramWritePoints[19]=point(0.320004, 0.00848248, 2e-05);
ramWritePoints[20]=point(0.320997, 0.00751752, 2e-05);
ramWritePoints[21]=point(0.320997, 0.0076874, 2e-05);
ramWritePoints[22]=point(0.320997, 0.00801027, 2e-05);
ramWritePoints[23]=point(0.320997, 0.00832447, 2e-05);
ramWritePoints[24]=point(0.320997, 0.00848248, 2e-05);
ramWritePoints[25]=point(0.322019, 0.00751752, 2e-05);
ramWritePoints[26]=point(0.322019, 0.0076874, 2e-05);
ramWritePoints[27]=point(0.322019, 0.00801027, 2e-05);
ramWritePoints[28]=point(0.322019, 0.00832447, 2e-05);
ramWritePoints[29]=point(0.322019, 0.00848248, 2e-05);
ramWritePoints[30]=point(0.323141, 0.00751752, 2e-05);
ramWritePoints[31]=point(0.323141, 0.0076874, 2e-05);
ramWritePoints[32]=point(0.323141, 0.00801027, 2e-05);
ramWritePoints[33]=point(0.323141, 0.00832447, 2e-05);
ramWritePoints[34]=point(0.323141, 0.00848248, 2e-05);

*/

//8:1
/*
List<point> ramWritePoints(nCell);                                                                            
ramWritePoints[0]=point(0.148515, 0.00398059, 5e-05);                                                         
ramWritePoints[1]=point(0.148515, 0.00320394, 5e-05);                                                         
ramWritePoints[2]=point(0.148515, 0.00200891, 5e-05);                                                          
ramWritePoints[3]=point(0.148515, 0.00088408, 5e-05);                                                         
ramWritePoints[4]=point(0.148515, 1.94139e-05, 5e-05);                                                        
ramWritePoints[5]=point(0.149002, 0.00398059, 5e-05);                                                         
ramWritePoints[6]=point(0.149002, 0.00320394, 5e-05);                                                         
ramWritePoints[7]=point(0.149002, 0.00200891, 5e-05);                                                         
ramWritePoints[8]=point(0.149002, 0.00088408, 5e-05);                                                         
ramWritePoints[9]=point(0.149002, 1.94139e-05, 5e-05);                                                        
ramWritePoints[10]=point(0.14951, 0.00398059, 5e-05);                                                         
ramWritePoints[11]=point(0.14951, 0.00320394, 5e-05);                                                         
ramWritePoints[12]=point(0.14951, 0.00200891, 5e-05);                                                         
ramWritePoints[13]=point(0.14951, 0.00088408, 5e-05);                                                         
ramWritePoints[14]=point(0.14951, 1.94139e-05, 5e-05);                                                        
ramWritePoints[15]=point(0.150002, 0.00224064, 5e-05);                                                        
ramWritePoints[16]=point(0.150002, 0.00214323, 5e-05);                                                        
ramWritePoints[17]=point(0.150002, 0.00200891, 5e-05);                                                        
ramWritePoints[18]=point(0.150002, 0.00186689, 5e-05);                                                        
ramWritePoints[19]=point(0.150002, 0.00175936, 5e-05);                                                        
ramWritePoints[20]=point(0.150508, 0.00224064, 5e-05);                                                        
ramWritePoints[21]=point(0.150508, 0.00214323, 5e-05);                                                        
ramWritePoints[22]=point(0.150508, 0.00200891, 5e-05);                                                        
ramWritePoints[23]=point(0.150508, 0.00186689, 5e-05);                                                        
ramWritePoints[24]=point(0.150508, 0.00175936, 5e-05);                                                        
ramWritePoints[25]=point(0.151011, 0.00224064, 5e-05);                                                        
ramWritePoints[26]=point(0.151011, 0.00214323, 5e-05);                                                        
ramWritePoints[27]=point(0.151011, 0.00200891, 5e-05);                                                        
ramWritePoints[28]=point(0.151011, 0.00186689, 5e-05);                                                        
ramWritePoints[29]=point(0.151011, 0.00175936, 5e-05);                                                        
ramWritePoints[30]=point(0.151484, 0.00224064, 5e-05);                                                        
ramWritePoints[31]=point(0.151484, 0.00214323, 5e-05);                                                        
ramWritePoints[32]=point(0.151484, 0.00200891, 5e-05);                                                        
ramWritePoints[33]=point(0.151484, 0.00186689, 5e-05);                                                        
ramWritePoints[34]=point(0.151484, 0.00175936, 5e-05);     
*/
//4:1
//
//
//

List<point> ramWritePoints(nCell);
ramWritePoints[0]=point(0.148502, 0.00197933, 5e-05);
ramWritePoints[1]=point(0.148502, 0.00150244, 5e-05);
ramWritePoints[2]=point(0.148502, 0.0010082, 5e-05);
ramWritePoints[3]=point(0.148502, 0.000515144, 5e-05);
ramWritePoints[4]=point(0.148502, 2.01043e-05, 5e-05);
ramWritePoints[5]=point(0.148988, 0.00197933, 5e-05);
ramWritePoints[6]=point(0.148988, 0.00150244, 5e-05);
ramWritePoints[7]=point(0.148988, 0.0010082, 5e-05);
ramWritePoints[8]=point(0.148988, 0.000515144, 5e-05);
ramWritePoints[9]=point(0.148988, 2.01043e-05, 5e-05);
ramWritePoints[10]=point(0.149503, 0.00197933, 5e-05);
ramWritePoints[11]=point(0.149503, 0.00150244, 5e-05);
ramWritePoints[12]=point(0.149503, 0.0010082, 5e-05);
ramWritePoints[13]=point(0.149503, 0.000515144, 5e-05);
ramWritePoints[14]=point(0.149503, 2.01043e-05, 5e-05);
ramWritePoints[15]=point(0.150002, 0.00123931, 5e-05);
ramWritePoints[16]=point(0.150002, 0.00113657, 5e-05);
ramWritePoints[17]=point(0.150002, 0.0010082, 5e-05);
ramWritePoints[18]=point(0.150002, 0.000872052, 5e-05);
ramWritePoints[19]=point(0.150002, 0.00076039, 5e-05);
ramWritePoints[20]=point(0.1505, 0.00123931, 5e-05);
ramWritePoints[21]=point(0.1505, 0.00113657, 5e-05);
ramWritePoints[22]=point(0.1505, 0.0010082, 5e-05);
ramWritePoints[23]=point(0.1505, 0.000872052, 5e-05);
ramWritePoints[24]=point(0.1505, 0.00076039, 5e-05);                                                          
ramWritePoints[25]=point(0.150995, 0.00123931, 5e-05);                                                        
ramWritePoints[26]=point(0.150995, 0.00113657, 5e-05);                                                        
ramWritePoints[27]=point(0.150995, 0.0010082, 5e-05);                                                         
ramWritePoints[28]=point(0.150995, 0.000872052, 5e-05);                                                       
ramWritePoints[29]=point(0.150995, 0.00076039, 5e-05);                                                        
ramWritePoints[30]=point(0.151454, 0.00123931, 5e-05);                                                        
ramWritePoints[31]=point(0.151454, 0.00113657, 5e-05);                                                        
ramWritePoints[32]=point(0.151454, 0.0010082, 5e-05);                                                         
ramWritePoints[33]=point(0.151454, 0.000872052, 5e-05);                                                       
ramWritePoints[34]=point(0.151454, 0.00076039, 5e-05);         


//find cell ID
forAll(ramWritePoints,n){
	label tmp_cell = mesh.findCell(ramWritePoints[n]);
	if(tmp_cell > -1){
		Pout<<"find cell id: "<<n<<endl;
		//initialize the struct
		//ramWriteList[ramWriteList_size] = new(ramWriteCell);
		ramWriteCell* ramWrite_tmp = &ramWriteList[ramWriteList_size++];
		ramWrite_tmp->cellID = tmp_cell;
		ramWrite_tmp->name = n;
		ramWrite_tmp->location = ramWritePoints[n];
		
		fileName ramWriteDir("ramWrite" + runTime.timeName());
		if(!isDir(ramWriteDir)){
			mkDir(ramWriteDir);
			Pout<< "mkdir "<<ramWriteDir<<endl;
		}
		ramWrite_tmp->fileP = new OFstream(ramWriteDir + "/cell_" + name(n) + "_P");
                ramWrite_tmp->fileP2 = new OFstream(ramWriteDir + "/cell_" + name(n) + "_P2");
                ramWrite_tmp->filedeltaP = new OFstream(ramWriteDir + "/cell_" + name(n) + "_deltaP");
		ramWrite_tmp->fileU = new OFstream(ramWriteDir + "/cell_" + name(n) + "_U");
                ramWrite_tmp->fileU2 = new OFstream(ramWriteDir + "/cell_" + name(n) + "_U2");
                ramWrite_tmp->filedeltaU = new OFstream(ramWriteDir + "/cell_" + name(n) + "_deltaU");
  		ramWrite_tmp->fileVorticity = new OFstream(ramWriteDir + "/cell_" + name(n) + "_Vorticity");
                ramWrite_tmp->fileVorticity2 = new OFstream(ramWriteDir + "/cell_" + name(n) + "_Vorticity2");
		ramWrite_tmp->fileTensorA = new OFstream(ramWriteDir + "/cell_" + name(n) + "_TensorA");
		ramWrite_tmp->fileSigmaP = new OFstream(ramWriteDir + "/cell_" + name(n) + "_SigmaP");
		ramWrite_tmp->fileSigma = new OFstream(ramWriteDir + "/cell_" + name(n) + "_Sigma");
		ramWrite_tmp->fileSigmas = new OFstream(ramWriteDir + "/cell_" + name(n) + "_Sigmas");
                ramWrite_tmp->fileSigmas2 = new OFstream(ramWriteDir + "/cell_" + name(n) + "_Sigmas2");
                ramWrite_tmp->filedeltaSigmas = new OFstream(ramWriteDir + "/cell_" + name(n) + "_deltaSigmas");
                ramWrite_tmp->fileGradU = new OFstream(ramWriteDir + "/cell_" + name(n) + "_GradU");
                ramWrite_tmp->fileGradU2 = new OFstream(ramWriteDir + "/cell_" + name(n) + "_GradU2");
		*(ramWrite_tmp->fileP) << "the data of P from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileP2) << "the data of P2 from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->filedeltaP) << "the data of deltaP from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		*(ramWrite_tmp->fileU) << "the data of U from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileU2) << "the data of U2 from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->filedeltaU) << "the data of deltaU from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		*(ramWrite_tmp->fileVorticity) << "the data of Vorticity from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileVorticity2) << "the data of Vorticity2 from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		*(ramWrite_tmp->fileTensorA) << "the data of Tensor from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		*(ramWrite_tmp->fileSigmaP) << "the data of SigmaP from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		*(ramWrite_tmp->fileSigma) << "the data of Sigma from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileSigmas) << "the data of Sigmas from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileSigmas2) << "the data of Sigmas2 from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->filedeltaSigmas) << "the data of deltaSigmas from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;		
                *(ramWrite_tmp->fileGradU) << "the data of GradU from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
                *(ramWrite_tmp->fileGradU2) << "the data of GradU2 from cell_"<<n<<", at location "<<ramWritePoints[n]<<endl;
		ramWrite_tmp->ramP = new List<scalar>(outputStep);
                ramWrite_tmp->ramP2 = new List<scalar>(outputStep);
                ramWrite_tmp->ramdeltaP = new List<scalar>(outputStep);
		ramWrite_tmp->ramU = new List<vector>(outputStep);
                ramWrite_tmp->ramU2 = new List<vector>(outputStep);
                ramWrite_tmp->ramdeltaU = new List<vector>(outputStep);
		ramWrite_tmp->ramVorticity= new List<vector>(outputStep);
                ramWrite_tmp->ramVorticity2= new List<vector>(outputStep);
		ramWrite_tmp->ramTensorA = new List<tensor>(outputStep);
		ramWrite_tmp->ramSigmaP = new List<tensor>(outputStep);
		ramWrite_tmp->ramSigma = new List<tensor>(outputStep);
                ramWrite_tmp->ramSigmas = new List<tensor>(outputStep);		
                ramWrite_tmp->ramSigmas2 = new List<tensor>(outputStep);
                ramWrite_tmp->ramdeltaSigmas = new List<tensor>(outputStep);
                ramWrite_tmp->ramGradU = new List<tensor>(outputStep);
                ramWrite_tmp->ramGradU2 = new List<tensor>(outputStep);
	}
}
