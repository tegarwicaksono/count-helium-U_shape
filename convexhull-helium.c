void BuildArrayHeAtoms(void) {
	int i, j;
	distHelium  = (		double **)malloc( NHELIUM*sizeof( double *) );
	for (i = 0; i < NHELIUM; i++) {
		distHelium[i] =  ( double *)malloc( NVERTEX*sizeof( double  ) );
		for (j = 0; j < NVERTEX; j++) 
			distHelium[i][j] = 0.0;
	}

	ClosestEdge  = (		double **)malloc( NHELIUM*sizeof( double *) );
	for (i = 0; i < NHELIUM; i++) {
		ClosestEdge[i] =  ( double *)malloc( 2*sizeof( double  ) );
		for (j = 0; j < 2; j++) 
			ClosestEdge[i][j] = 0.0;
	}

	LabelHeAtoms = ( int *)malloc( NHELIUM*sizeof( int ) );

	labelBULK 	 = ( int *)malloc( NHELIUM*sizeof( int ) );
	labelBULKOUT = ( int *)malloc( NHELIUM*sizeof( int ) );
	labelEDGE		 = ( int *)malloc( NHELIUM*sizeof( int ) );
	labelLOOP		 = ( int *)malloc( NHELIUM*sizeof( int ) );

	for (i = 0; i < NHELIUM; i++) {
		LabelHeAtoms[i] = 0;
		labelBULK[i] = 0;
		labelBULKOUT[i] = 0;
		labelEDGE[i] = 0;
		labelLOOP[i] = 0;
	}
}

void DistFromLineEdge(int index, double coor[]) {
	int 	i, j;
	double x = coor[0], y = coor[1];
	double P0[2], P1[2], P2[2];
	double label, slope, distance, dx, dy, dnewx, dnewy, t;
	double yinter;
	for (i = 0; i < NVERTEX; i++) {
		for (j = 0; j < 2; j++) {
			P0[j] = coor[j];
			P1[j] = VXPOINT[(i)%NVERTEX][j];
			P2[j] = VXPOINT[(i+1)%NVERTEX][j];
		}

		dx = P2[0] - P1[0];		dy = P2[1] - P1[1];

		t  = ( (P0[0]-P1[0])*dx + (P0[1]-P1[1])*dy) / (dx*dx + dy*dy);

		if (t < 0.0) {
			dx = P0[0] - P1[0];		dy = P0[1] - P1[1];
		} else if (t > 1.0) {
			dx = P0[0] - P2[0];		dy = P0[1] - P2[1];
		} else {
			dnewx = P1[0] + t*dx; dnewy = P1[1] + t*dy;
			dx = P0[0] - dnewx;	dy = P0[1] - dnewy;
		}

		distHelium[index][i] = sqrt(dx*dx + dy*dy);
	}

}

void FindMinDistance (int index) {
	double minvalue = distHelium[index][0];
	int i, indmin = 1;
	for (i = 1; i < NVERTEX; i++) {
		if (minvalue > distHelium[index][i]) {
			minvalue = distHelium[index][i];
			indmin = i+1;
		}
	}
	ClosestEdge[index][0] = minvalue;
	ClosestEdge[index][1] = (double)indmin;

}
	
void CalcDistHeAtoms(void) {
	int i,j;
	double coor[3];

  for (i = 0; i < NCLUSTERS; i++) {
		for (j = 0; j < 3; j++) coor[j] = clusCOM[i][j];
		DistFromLineEdge(i, coor);
		FindMinDistance(i);
  }

	fclose(heliuminput);
}

void PrintHeliumDistance(void) {
	int i, j;
  char nameHeCFG[]  = "./Appendix/Dist-He/He-Distance.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.txt";
  char charHeCFG[sizeof nameHeCFG+200];
	sprintf(charHeCFG,nameHeCFG,TEMP,NHELIUM,NINIT,timestep);
  heliumcfg = fopen(charHeCFG,"w"); 	\
		for(j = 0; j < NVERTEX; j++) {
			fprintf(heliumcfg, "edge-%d\t", j);
		}
		fprintf(heliumcfg, "min-dist\t");
		fprintf(heliumcfg, "clos-edge\t");
		fprintf(heliumcfg, "label\n");

	for (i = 0; i < NCLUSTERS; i++) {
		for(j = 0; j < NVERTEX; j++) {
			fprintf(heliumcfg, "%.5E\t",distHelium[i][j]);
		}
		fprintf(heliumcfg, "%.5E\t", ClosestEdge[i][0]);
		fprintf(heliumcfg, "%d\t", (int)ClosestEdge[i][1]);
		fprintf(heliumcfg, "%d\t", LabelHeAtoms[i]);
		fprintf(heliumcfg, "\n");
	}
	fclose(heliumcfg);
}

double whichmax(double a, double b) {
	return (a > b ? a : b);
}

double whichmin(double a, double b) {
	return (a < b ? a : b);
}

int CheckIfInsideShrinkingGrain(double P0[]){
	int counter = 0;
	int i, j;
	double xinters;
	double P1[2], P2[2];

	for (j = 0; j < 2; j++) {
		if (P0[j] < 0.0) P0[j] += LENGTH[j];
		else if (P0[j] > LENGTH[j]) P0[j] -= LENGTH[j];
	}

	for (j = 0; j < 2; j++) P1[j] = VXPOINT[0][j];
	for (i = 1; i <= NVERTEX; i++) {
		for (j = 0; j < 2; j++) P2[j] = VXPOINT[(i%NVERTEX)][j];
		if (P0[1] > whichmin(P1[1], P2[1])) {
			if (P0[1] <= whichmax(P1[1], P2[1])) {
				if (P0[0] <= whichmax(P1[0], P2[0])) {
					if (P1[1] != P2[1]) {
						xinters = (P0[1]-P1[1])*(P2[0] - P1[0])/(P2[1]-P1[1]) + P1[0];
						if (P1[0] == P2[0] || P0[0] <= xinters) {
							counter++;
						}
					}
				}
			}
		}
		for (j = 0; j < 2; j++) P1[j] = P2[j];
	}

	if (counter % 2 == 0) 
		return 0; //Outside
	else
		return 1;	//Inside
}
 

void BinHeAtoms(void) {
	int i, j, label, clustersize, inside;
	double coor[3];

	NHELIUMBULK = 0;
	NHELIUMEDGE = 0;
	NHELIUMLOOP	= 0;
	NHELIUMBULKOUT = 0;

	for (i = 0; i < NCLUSTERS; i++) {
		if (ClosestEdge[i][0] < CUTOFFDISTANCE) {
			label = (int)ClosestEdge[i][1] - 1;
			if (HULLEDGE[label][0] == -1.0) {
				LabelHeAtoms[i] = 0;
				NHELIUMBULK += clusSIZE[clusSeq[i]];
			} else if (HULLEDGE[label][0] == 1.0) {
				LabelHeAtoms[i] = 1;
				NHELIUMEDGE += clusSIZE[clusSeq[i]];
			} else {
				LabelHeAtoms[i] = 2;
				NHELIUMLOOP += clusSIZE[clusSeq[i]];
			}
		} else {
			for (j = 0; j < 3; j++) coor[j] = clusCOM[i][j];
			inside = CheckIfInsideShrinkingGrain(coor);
			if (inside == 0) {
				LabelHeAtoms[i] = -1;
				NHELIUMBULKOUT += clusSIZE[clusSeq[i]];
			} else {
				LabelHeAtoms[i] = 0;
				NHELIUMBULK += clusSIZE[clusSeq[i]];
			}
		}
	}

	for (i = 0; i < NCLUSTERS; i++) {
		clustersize = clusSIZE[clusSeq[i]] - 1;
		if (LabelHeAtoms[i] == -1) {
			labelBULKOUT[clustersize]++;
		} else if (LabelHeAtoms[i] == 0) {
			labelBULK[clustersize]++;
		} else if (LabelHeAtoms[i] == 1) {
			labelEDGE[clustersize]++;
		} else if (LabelHeAtoms[i] == 2) {
			labelLOOP[clustersize]++;
		}
	}
		
}

void PrintCFGHeliumHeader(void) {
  fprintf(heliumcfg,"Number of particles = %d\n", NHELIUM+NVERTEX);
  fprintf(heliumcfg,"A = 1 Angstrom (basic length-scale)\n");
  fprintf(heliumcfg,"H0(1,1) = %lf A\n",LENGTH[0]);
  fprintf(heliumcfg,"H0(1,2) = 0 A\n");
  fprintf(heliumcfg,"H0(1,3) = 0 A\n");
  fprintf(heliumcfg,"H0(2,1) = 0 A\n");
  fprintf(heliumcfg,"H0(2,2) = %lf A\n",LENGTH[1]);
  fprintf(heliumcfg,"H0(2,3) = 0 A\n");
  fprintf(heliumcfg,"H0(3,1) = 0 A\n");
  fprintf(heliumcfg,"H0(3,2) = 0 A\n");
  fprintf(heliumcfg,"H0(3,3) = %lf A\n",LENGTH[2]);
  fprintf(heliumcfg,".NO_VELOCITY.\n");
  fprintf(heliumcfg,"entry_count = 8\n");
	fprintf(heliumcfg, "auxiliary[0] = mindistance\n");
	fprintf(heliumcfg, "auxiliary[1] = whichedge\n");
	fprintf(heliumcfg, "auxiliary[2] = clustersize\n");
	fprintf(heliumcfg, "auxiliary[3] = clusterid\n");
	fprintf(heliumcfg, "auxiliary[4] = bulk_gb\n");
  fprintf(heliumcfg,"4.0026\n");
  fprintf(heliumcfg,"Ga\n");
}

void PrintLineHeAtom(int index) {
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][0] / LENGTH[0]);
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][1] / LENGTH[1]);
	fprintf(heliumcfg,"%.6lf\t", atomCOOR[index][2] / LENGTH[2]);
	fprintf(heliumcfg,"%.6lf\t",ClosestEdge[revAtomName[index]][0]);
	fprintf(heliumcfg,"%.0lf\t",ClosestEdge[revAtomName[index]][1]);
	fprintf(heliumcfg,"%d\t",	clusSIZE[leadSORT[index]]);
	fprintf(heliumcfg,"%d\t", revAtomName[leadSORT[index]]);
	fprintf(heliumcfg,"%d\n",	LabelHeAtoms[revAtomName[index]]);
}			

void PrintLineVertex(int index) {
	fprintf(heliumcfg,"%.6lf\t", VXPOINT[index][0] / LENGTH[0]);
	fprintf(heliumcfg,"%.6lf\t", VXPOINT[index][1] / LENGTH[1]);
	fprintf(heliumcfg,"%.6lf\t", 0.5);
	fprintf(heliumcfg,"%.6lf\t",0.0);
	fprintf(heliumcfg,"%d\t",index);
	fprintf(heliumcfg,"%d\t",0);
	fprintf(heliumcfg,"%d\t",-1);
	fprintf(heliumcfg,"%d\n",3);
}

void PrintCFGBinHelium(void) {
	int i;
  char nameHeCFG[]  = "./CFG-cluster/CFG.BinHe.COH.T%dK.DIA50.NHE.ALL%d.INIT%d-CUTOFFRAD_CLUS%.0lf_BULKGB%.1lf.%d.cfg";
  char charHeCFG[sizeof nameHeCFG+200];
	sprintf(charHeCFG,nameHeCFG,TEMP,NHELIUM,NINIT,CUTOFFRAD,CUTOFFDISTANCE,timestep);
  heliumcfg = fopen(charHeCFG,"w");

	PrintCFGHeliumHeader();


	for (i = 0; i < NHELIUM; i++) {

		PrintLineHeAtom(i);
	}

  fprintf(heliumcfg,"55.85\n");
  fprintf(heliumcfg,"Fe\n");

	for (i = 0; i < NVERTEX; i++) {
		PrintLineVertex(i);
	}

	fclose(heliumcfg);

}

void StartEvolSizeAll(void) {
	sprintf(charHeEvol,nameHeEvol,TEMP,NHELIUM,NINIT,(double)inistep/1.0e6, (double)finstep/1.0e6);
  evol = fopen(charHeEvol,"w");
	fprintf(evol, "time[ns]\tNumberEdge\tBulkInsideGrain\tBukOutsideShrinkingGrain\tEdgeHe\tLoopHe\n");
}

void StartEvolSizeBulkInside(void) {
	int i;
	sprintf(charHeEvolSizeBulkInside,nameHeEvolSizeBulkInside,TEMP,NHELIUM,NINIT,(double)inistep/1.0e6, (double)finstep/1.0e6);
  evolsizebulkin = fopen(charHeEvolSizeBulkInside,"w");
	fprintf(evolsizebulkin, "time[ns]\t");
	for (i = 0; i < NSIZEPRINT; i++) {
		fprintf(evolsizebulkin,"Size%d\t",i+1);
	}
	fprintf(evolsizebulkin,"\n");
}

void StartEvolSizeBulkOutside(void) {
	int i;
	sprintf(charHeEvolSizeBulkOutside,nameHeEvolSizeBulkOutside,TEMP,NHELIUM,NINIT,(double)inistep/1.0e6, (double)finstep/1.0e6);
  evolsizebulkout = fopen(charHeEvolSizeBulkOutside,"w");
	fprintf(evolsizebulkout, "time[ns]\t");
	for (i = 0; i < NSIZEPRINT; i++) {
		fprintf(evolsizebulkout,"Size%d\t",i+1);
	}
	fprintf(evolsizebulkout,"\n");
}

void StartEvolSizeEdge(void) {
	int i;
	sprintf(charHeEvolSizeEdge,nameHeEvolSizeEdge,TEMP,NHELIUM,NINIT,(double)inistep/1.0e6, (double)finstep/1.0e6);
  evolsizeedge = fopen(charHeEvolSizeEdge,"w");
	fprintf(evolsizeedge, "time[ns]\t");
	for (i = 0; i < NSIZEPRINT; i++) {
		fprintf(evolsizeedge,"Size%d\t",i+1);
	}
	fprintf(evolsizeedge,"\n");
}

void StartEvolSizeLoop(void) {
	int i;
	sprintf(charHeEvolSizeLoop,nameHeEvolSizeLoop,TEMP,NHELIUM,NINIT,(double)inistep/1.0e6, (double)finstep/1.0e6);
  evolsizeloop = fopen(charHeEvolSizeLoop,"w");
	fprintf(evolsizeloop, "time[ns]\t");
	for (i = 0; i < NSIZEPRINT; i++) {
		fprintf(evolsizeloop,"Size%d\t",i+1);
	}
	fprintf(evolsizeloop,"\n");
}

void StartEvolution(void) {
	StartEvolSizeAll();
	StartEvolSizeBulkInside();
	StartEvolSizeBulkOutside();
	StartEvolSizeEdge();
	StartEvolSizeLoop();
}

void PrintEvolution(void) {
	int i;
	fprintf(evol,"%.4lf\t", (double)timestep / 1.0e6);
	fprintf(evol,"%d\t",	NVERTEX);
	fprintf(evol,"%d\t",	NHELIUMBULK);
	fprintf(evol,"%d\t",	NHELIUMBULKOUT);
	fprintf(evol,"%d\t",  NHELIUMEDGE);
	fprintf(evol,"%d\t",  NHELIUMLOOP);
	fprintf(evol,"\n");

		fprintf(evolsizebulkin,	"%.4lf\t", (double)timestep / 1.0e6);
		fprintf(evolsizebulkout,"%.4lf\t", (double)timestep / 1.0e6);
		fprintf(evolsizeedge,		"%.4lf\t", (double)timestep / 1.0e6);
		fprintf(evolsizeloop,		"%.4lf\t", (double)timestep / 1.0e6);


	for (i = 0; i < NSIZEPRINT; i++) {
		fprintf(evolsizebulkin,	"%d\t",labelBULK[i]);
		fprintf(evolsizebulkout,"%d\t",labelBULKOUT[i]);
		fprintf(evolsizeedge,		"%d\t",labelEDGE[i]);
		fprintf(evolsizeloop,		"%d\t",labelLOOP[i]);
	}
		fprintf(evolsizebulkin,	"\n");
		fprintf(evolsizebulkout,"\n");
		fprintf(evolsizeedge,		"\n");
		fprintf(evolsizeloop,		"\n");
}

void CloseEvolution(void) {
	fclose(evol);
	fclose(evolsizebulkin);
	fclose(evolsizebulkout);
	fclose(evolsizeedge);
	fclose(evolsizeloop);
}

