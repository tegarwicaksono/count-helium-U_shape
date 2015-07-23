void StartValidity(void) {
  char nameDist[] = "./Appendix/Valid-Report/ValidRep.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.txt";
  char charDist[sizeof nameDist+200];

  sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD);
  validity = fopen(charDist,"w");

	fprintf(validity,"Timestep\tCorrelation Coefficient\n");
}

void CloseValidity(void) {
  fclose(validity);
}

void PerformValidity(void) {
  double SumAtomDist = 0.0;
	double SumCophDist = 0.0;
	double SumAtomDistSquare = 0.0;
	double SumCophDistSquare = 0.0;
	double SumAtomCophDist	 = 0.0;

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			SumAtomDist += atomDIST[ij][kl];
			SumCophDist += cophDIST[ij][kl];
			SumAtomDistSquare += (atomDIST[ij][kl]*atomDIST[ij][kl]);
			SumCophDistSquare += (cophDIST[ij][kl]*cophDIST[ij][kl]);
			SumAtomCophDist		+= (atomDIST[ij][kl]*cophDIST[ij][kl]);
		}
	}

	CorrCoefficient  = SumAtomCophDist - (SumAtomDist*SumCophDist)/NPAIRDIST;
	CorrCoefficient /= sqrt( (SumAtomDistSquare - (SumAtomDist*SumAtomDist)/NPAIRDIST));
	CorrCoefficient /= sqrt( (SumCophDistSquare - (SumCophDist*SumCophDist)/NPAIRDIST));
	fprintf(validity,"%d\t",timestep);
	fprintf(validity,"%.6lf\t",CorrCoefficient);
	fprintf(validity,"\n");
}

void PerformAgglomerateCluster(void) {
	double act, ref, dif;
	NCLUSTERS = 0;
	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		if (leadSORT[ij] == ij) {
			for (kl = (ij+1); kl < NHELIUM; kl++) {
				if (cophDIST[ij][kl] < CUTOFFRAD) {
					leadSORT[kl] = ij;
				}
			}
		}
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		if (leadSORT[ij] == ij) {
			revAtomName[ij] = NCLUSTERS;	//revAtomName range from 0 to NHELIUM
			clusSeq[NCLUSTERS] = ij;			//clusSeq range from 0 to NCLUSTERS
			for (kl = 0; kl < 3; kl++) {
				clusCOM[NCLUSTERS][kl] += atomCOOR[ij][kl];	//clusCentreofMass range from 0 to NCLUSTERS
			}
			
			NCLUSTERS++;
		} else {
			revAtomName[ij] = revAtomName[leadSORT[ij]];

			for (kl = 0; kl < 3; kl++) {
				act = atomCOOR[ij][kl];
				ref = atomCOOR[leadSORT[ij]][kl];
				dif = act - ref;
				if ( (revAtomName[ij] == 3)  && (kl == 0)) {
				} 
				if (dif > 0.5*LENGTH[kl])      	 	 act  = act - LENGTH[kl];
				else if (dif < (-0.5*LENGTH[kl]) ) act  = act + LENGTH[kl];

				clusCOM[revAtomName[leadSORT[ij]]][kl] += act;
			}
		}
	}

					
	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		if (clusSIZE[ij] > 0) {
			for (kl = (ij+1); kl < NHELIUM; kl++) {
				if (cophDIST[ij][kl] < CUTOFFRAD) {			
					clusSIZE[ij] += 1;
					clusSIZE[kl]  = 0;
				}
			}
		}
	}

	for (ij = 0; ij < NHELIUM; ij++) {
		if (clusSIZE[ij] > 0) {
			clusPOPU[clusSIZE[ij]-1] += 1;
			clusPOPU[0] -= clusSIZE[ij];
		}
	}

	for (ij = 0; ij < NCLUSTERS; ij++) {
		for (kl = 0; kl < 3; kl++) {
			clusCOM[ij][kl] = clusCOM[ij][kl] / (double)clusSIZE[clusSeq[ij]];
		}
	}
}

void PrintDistMatrix(void) {
  char nameDist[] = "./Appendix/Dist-Matrix/DistanceMatrix.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			fprintf(distmatrix,"%4d\t%4d\t%.8lf\n", ij, kl, atomDIST[ij][kl]);
		}
	}

	fclose(distmatrix);
}

void PrintCophMatrix(void) {
  char nameDist[] = "./Appendix/Coph-Matrix/CopheneticMatrix.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUMMINUSONE; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			fprintf(distmatrix,"%4d\t%4d\t%.8lf\t%.8lf\n", ij, kl, cophDIST[ij][kl], cophSIGN[ij][kl]);
		}
	}

	fclose(distmatrix);
}

void PrintClusterPopulation(void) {
  char nameDist[] = "./Clust-Distribution/ClustDistribution.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");
	fprintf(distmatrix,"Size\tTotalPopulation\tBULK\tBULKOUTSIDE\tEDGE\tLOOP\n");

	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t", ij+1);
		fprintf(distmatrix," %d\t", clusPOPU[ij]);
		fprintf(distmatrix," %d\t", labelBULK[ij]);
		fprintf(distmatrix," %d\t", labelBULKOUT[ij]);
		fprintf(distmatrix," %d\t", labelEDGE[ij]);
		fprintf(distmatrix," %d\t", labelLOOP[ij]);
		fprintf(distmatrix,"\n");
	}

	fclose(distmatrix);
}

void PrintAtomSize(void) {
  char nameDist[] = "./Appendix/Atom-Size/ATOMSIZE.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t%4d\n", ij,atomSIZE[ij]);
	}

	fclose(distmatrix);
}

void PrintAtomCoordinate(void) {
  char nameDist[] = "./Appendix/Atom-Coordinate/ATOMCoordinate.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");

	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t%.8lf\t%.8lf\t%.8lf\n", ij,atomCOOR[ij][0],atomCOOR[ij][1],atomCOOR[ij][2]);
	}

	fclose(distmatrix);
}

void PrintAtomLead(void) {
  char nameDist[] = "./Appendix/Atom-Lead/ATOMLEAD.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	int yes;
	FILE *distmatrix;

	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  distmatrix = fopen(charDist,"w");
	fprintf(distmatrix,"Index\tX\tY\tZ\tLeadIndx\tSIZE\tLeader?\n");
	for (ij = 0; ij < NHELIUM; ij++) {
		fprintf(distmatrix,"%4d\t",ij);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][0]);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][1]);
		fprintf(distmatrix,"%.8lf\t",atomCOOR[ij][2]);
		fprintf(distmatrix,"%4d\t",leadSORT[ij]);
		fprintf(distmatrix,"%4d\t",clusSIZE[leadSORT[ij]]);
		yes = (leadSORT[ij] == ij) ? 1 : 0;
		fprintf(distmatrix,"%4d\t",yes);
		fprintf(distmatrix,"\n");
	}

	fclose(distmatrix);
}

void PrintCFGClustCOMHeader(void) {
  fprintf(cfgnew,"Number of particles = %d\n", NCLUSTERS);
  fprintf(cfgnew,"A = 1 Angstrom (basic length-scale)\n");
  fprintf(cfgnew,"H0(1,1) = %lf A\n",LENGTH[0]);
  fprintf(cfgnew,"H0(1,2) = 0 A\n");
  fprintf(cfgnew,"H0(1,3) = 0 A\n");
  fprintf(cfgnew,"H0(2,1) = 0 A\n");
  fprintf(cfgnew,"H0(2,2) = %lf A\n",LENGTH[1]);
  fprintf(cfgnew,"H0(2,3) = 0 A\n");
  fprintf(cfgnew,"H0(3,1) = 0 A\n");
  fprintf(cfgnew,"H0(3,2) = 0 A\n");
  fprintf(cfgnew,"H0(3,3) = %lf A\n",LENGTH[2]);
  fprintf(cfgnew,".NO_VELOCITY.\n");
  fprintf(cfgnew,"entry_count = 6\n");
  fprintf(cfgnew,"auxiliary[0] = index\n");
  fprintf(cfgnew,"auxiliary[1] = leader\n");
  fprintf(cfgnew,"auxiliary[2] = size\n");
  fprintf(cfgnew,"4.008\n");
  fprintf(cfgnew,"Ga\n");
}

void PrintCFGClusterCentreOfMass(void) {
  char nameDist[] = "./Appendix/CFG-CentOfMass/CFG.COM.Clust.CURV.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.CUTOFFRAD%.0lf.%d.txt";
  char charDist[sizeof nameDist+200];
	sprintf(charDist,nameDist,TEMP,NHELIUM,NINIT,CUTOFFRAD,timestep);
  cfgnew = fopen(charDist,"w");

	PrintCFGClustCOMHeader();
	for (ij = 0; ij < NCLUSTERS; ij++) {
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][0]/LENGTH[0]);
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][1]/LENGTH[1]);
		fprintf(cfgnew,"%.8lf\t",clusCOM[ij][2]/LENGTH[2]);

		fprintf(cfgnew,"%d\t",ij);
		fprintf(cfgnew,"%d\t",clusSeq[ij]);
		fprintf(cfgnew,"%d\t",clusSIZE[clusSeq[ij]]);

		fprintf(cfgnew,"\n");	
	}	
}


