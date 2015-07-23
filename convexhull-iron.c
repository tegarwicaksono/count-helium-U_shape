
void OpenFeAtomicFile(void) {
	int i;
		sprintf(charFeInput,nameFeInput,TEMP,NHELIUM,NINIT,timestep);
    ironinput = fopen(charFeInput,"r");

		fscanf(ironinput,"%s %s %s %s %d", &w1, &w2, &w3,&w4, &NIRON);
		fclose(ironinput);

		sprintf(charFeInput,nameFeInput,TEMP,NHELIUM,NINIT,timestep);
    ironinput = fopen(charFeInput,"r");
    for (i = 0; i < 2; i++) {
      do
        skip = fgetc(ironinput);
      while (skip != '\n');
    }
    fscanf(ironinput,"%s %s", &ig1, &ig2);
    fscanf(ironinput,"%lf", &LENGTH[0]);
    fscanf(ironinput,"%s", &ig3);
    fclose(ironinput);

		sprintf(charFeInput,nameFeInput,TEMP,NHELIUM,NINIT,timestep);
    ironinput = fopen(charFeInput,"r");
    for (i = 0; i < 6; i++) {
      do
        skip = fgetc(ironinput);
      while (skip != '\n');
    }
    fscanf(ironinput,"%s %s", &ig1, &ig2);
    fscanf(ironinput,"%lf", &LENGTH[1]);
    fscanf(ironinput,"%s", &ig3);
    fclose(ironinput);

		sprintf(charFeInput,nameFeInput,TEMP,NHELIUM,NINIT,timestep);
    ironinput = fopen(charFeInput,"r");
    for (i = 0; i < 10; i++) {
      do
        skip = fgetc(ironinput);
      while (skip != '\n');
    }
    fscanf(ironinput,"%s %s", &ig1, &ig2);
    fscanf(ironinput,"%lf", &LENGTH[2]);
    fscanf(ironinput,"%s", &ig3);
    fclose(ironinput);
}

void BuildArrayFeAtoms(void) {
	int i;
	LEFTCAND 	= (	  double * )malloc( NIRON*sizeof( double  ) );
	ZLEFT			= (	  double * )malloc( NIRON*sizeof( double  ) );
	RIGHTCAND = (	  double * )malloc( NIRON*sizeof( double  ) );
	ZRIGHT		= (	  double * )malloc( NIRON*sizeof( double  ) );
	LOOPCAND  = (		double **)malloc( NIRON*sizeof( double *) );
	for (i = 0; i < NIRON; i++) {
		LOOPCAND[i] =  ( double *)malloc( 3*sizeof( double  ) );
	}
}

void ReadFeAtoms(void) {
	int i;
	sprintf(charFeInput,nameFeInput,TEMP,NHELIUM,NINIT,timestep);
  ironinput = fopen(charFeInput,"r");

  for (i = 0; i < LINEFE; i++) {
    do
      skip = fgetc(ironinput);
    while (skip != '\n');
  }
}

void ScanForIronCoordinate(void) {
	int op;
	fscanf(ironinput,"%lf", &ATOMCOOR[0]);
  fscanf(ironinput,"%lf", &ATOMCOOR[1]);
  fscanf(ironinput,"%lf", &ATOMCOOR[2]);

	for (op = 0; op < ENTRYFE; op++) {
	  fscanf(ironinput,"%lf", &ENTRYCFG);
	}
}

void ScanBinIronAtoms(void) {
	int i;
	ReadFeAtoms();
	double maxleft = -LENGTH[1], minleft = LENGTH[1];
	double maxright = -LENGTH[1], minright = LENGTH[1];	

	NLCAND = 0; NRCAND = 0; NMCAND = 0;
	for (i = 0; i < NIRON; i++) {
		ScanForIronCoordinate();

		if ((ATOMCOOR[0] >= LEFTEDGE[0]) && (ATOMCOOR[0] <= LEFTEDGE[1])) {
			LEFTCAND[NLCAND] = ATOMCOOR[1]*LENGTH[1];
			ZLEFT[NLCAND]    = ATOMCOOR[2]*LENGTH[2];
			if (LEFTCAND[NLCAND] > maxleft)	{
				maxleft = LEFTCAND[NLCAND];
				idMAXLEFT = NLCAND;
			}

			if (LEFTCAND[NLCAND] < minleft) {
				minleft = LEFTCAND[NLCAND];
				idMINLEFT = NLCAND;
			}
			NLCAND++;
		} else if ((ATOMCOOR[0] >= RIGHTEDGE[0]) && (ATOMCOOR[0] <= RIGHTEDGE[1])) {

			RIGHTCAND[NRCAND] = ATOMCOOR[1]*LENGTH[1];
			ZRIGHT[NRCAND]    = ATOMCOOR[2]*LENGTH[2];
			if (RIGHTCAND[NRCAND] > maxright)	{
				maxright = RIGHTCAND[NLCAND];
				idMAXRIGHT = NRCAND;
			}

			if (RIGHTCAND[NRCAND] < minright) {
				minright = RIGHTCAND[NLCAND];
				idMINRIGHT = NRCAND;
			}
			NRCAND++;
		} else {
			LOOPCAND[NMCAND][0] = ATOMCOOR[0]*LENGTH[0];
			LOOPCAND[NMCAND][1] = ATOMCOOR[1]*LENGTH[1];
			LOOPCAND[NMCAND][2] = ATOMCOOR[2]*LENGTH[2];
			NMCAND++;
		}
	}

	fclose(ironinput);
}

void SortLeftRightCand(void) {
	LEFTCAND[idMINLEFT]  = upperBorder;
	RIGHTCAND[idMINRIGHT] = upperBorder;
}

void AssignVertexCandidate(void) {
	int i, j;
	NCAND = NMCAND + 4;
	VXCAND  = (		double **)malloc( NCAND*sizeof( double *) );
	for (i = 0; i < NCAND; i++) {
		VXCAND[i] =  ( double *)malloc( 2*sizeof( double  ) );
	}

	VXCAND[0][0] = LEFTLOC*LENGTH[0];	VXCAND[0][1] = LEFTCAND[idMINLEFT];
	VXCAND[1][0] = LEFTLOC*LENGTH[0]; VXCAND[1][1] = LEFTCAND[idMAXLEFT];

	for (i = 0; i < NMCAND; i++) {
		for (j = 0; j < 2; j++) {
			VXCAND[i+2][j] = LOOPCAND[i][j];
		}
	}
	VXCAND[NCAND-2][0] = RIGHTLOC*LENGTH[0]; VXCAND[NCAND-2][1] = RIGHTCAND[idMAXRIGHT];
	VXCAND[NCAND-1][0] = RIGHTLOC*LENGTH[0]; VXCAND[NCAND-1][1] = RIGHTCAND[idMINRIGHT];

}

void FindVertexCandidate(void) {
	ScanBinIronAtoms();
	SortLeftRightCand();
	AssignVertexCandidate();
}



void PrintCFGHeader(void) {
  fprintf(ironCand,"Number of particles = %d\n", NMCAND + 4);
  fprintf(ironCand,"A = 1 Angstrom (basic length-scale)\n");
  fprintf(ironCand,"H0(1,1) = %lf A\n",LENGTH[0]);
  fprintf(ironCand,"H0(1,2) = 0 A\n");
  fprintf(ironCand,"H0(1,3) = 0 A\n");
  fprintf(ironCand,"H0(2,1) = 0 A\n");
  fprintf(ironCand,"H0(2,2) = %lf A\n",LENGTH[1]);
  fprintf(ironCand,"H0(2,3) = 0 A\n");
  fprintf(ironCand,"H0(3,1) = 0 A\n");
  fprintf(ironCand,"H0(3,2) = 0 A\n");
  fprintf(ironCand,"H0(3,3) = %lf A\n",LENGTH[2]);
  fprintf(ironCand,".NO_VELOCITY.\n");
  fprintf(ironCand,"entry_count = 3\n");
  fprintf(ironCand,"55.85\n");
  fprintf(ironCand,"Fe\n");
}

void PrintVertexCandidate(void) {
	int i;
	sprintf(charFeCand,nameFeCand,TEMP,NHELIUM,NINIT,timestep);
  ironCand = fopen(charFeCand,"w");
	fprintf(ironCand,"%.4E\t%.4E\n", LEFTLOC*LENGTH[0], LEFTCAND[idMINLEFT]);
	fprintf(ironCand,"%.4E\t%.4E\n", LEFTLOC*LENGTH[0], LEFTCAND[idMAXLEFT]);
	for (i = 0; i < NMCAND; i++) {
		fprintf(ironCand,"%.4E\t%.4E\n", LOOPCAND[i][0],  LOOPCAND[i][1]);
	}
	fprintf(ironCand,"%.4E\t%.4E\n", RIGHTLOC*LENGTH[0], RIGHTCAND[idMAXRIGHT]);
	fprintf(ironCand,"%.4E\t%.4E\n", RIGHTLOC*LENGTH[0], RIGHTCAND[idMINRIGHT]);	
	fclose(ironCand);
}

void PrintCFGCandidate(void) {
	int i, j;

  char nameFeCandCFG[]  = "./Cand-CFG-Fe/CFG.FeCand.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.cfg";
  char charFeCandCFG[sizeof nameFeCand+200];
	sprintf(charFeCandCFG,nameFeCandCFG,TEMP,NHELIUM,NINIT,timestep);
  ironCand = fopen(charFeCandCFG,"w");

	PrintCFGHeader();
	
	fprintf(ironCand,"%.4E\t%.4E\t%.4E\n", LEFTLOC, LEFTCAND[idMINLEFT] / LENGTH[1], ZLEFT[idMINLEFT] / LENGTH[2]);
	fprintf(ironCand,"%.4E\t%.4E\t%.4E\n", LEFTLOC, LEFTCAND[idMAXLEFT] / LENGTH[1], ZLEFT[idMAXLEFT] / LENGTH[2]);
	for (i = 0; i < NMCAND; i++) {
		for (j = 0; j < 3; j++) {
			fprintf(ironCand,"%.4E\t", LOOPCAND[i][j] / LENGTH[j]);
		}
		fprintf(ironCand,"\n");
	}
	fprintf(ironCand,"%.4E\t%.4E\t%.4E\n", RIGHTLOC, RIGHTCAND[idMAXRIGHT] / LENGTH[1], ZRIGHT[idMAXRIGHT]/ LENGTH[2]);
	fprintf(ironCand,"%.4E\t%.4E\t%.4E\n", RIGHTLOC, RIGHTCAND[idMINRIGHT] / LENGTH[1], ZRIGHT[idMINRIGHT]/ LENGTH[2]);	
	fclose(ironCand);
}	



