void OpenHeAtomicFile(void) {
	int i;
	sprintf(charHeInput,nameHeInput,TEMP,NHELIUM,NINIT,timestep);
  heliuminput = fopen(charHeInput,"r");

  for (i = 0; i < LINEHE; i++) {
    do
      skip = fgetc(heliuminput);
    while (skip != '\n');
  }
}

void	BuildArrayClusters(void) {
	NPAIRDIST = NHELIUM*(NHELIUM - 1) / 2;
	CUTOFFRADSQ = CUTOFFRAD*CUTOFFRAD;

	for (ij = 0; ij < 3; ij++)
		HALFLENGTH[ij] = 0.5*LENGTH[ij];

  atomNAME = (		int  	 *)malloc(	NHELIUM*sizeof(	int	)	);
	for (ij = 0; ij < NHELIUM; ij++)
		atomNAME[ij] = 0;

	leadSORT = (		int		 *)malloc(	NHELIUM*sizeof(	int ) );
	for (ij = 0; ij < NHELIUM; ij++)
		leadSORT[ij] = ij;
	
	atomCOOR = (		double **)malloc(	NHELIUM*sizeof(	double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		atomCOOR[ij] = ( double *)malloc(3*sizeof(double));
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < 3; kl++)
			atomCOOR[ij][kl] = 0.0;

	atomDIST = (		double **)malloc(	NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		atomDIST[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			atomDIST[ij][kl] = 0.0;

	clusCOM = (		double **)malloc(	NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++) 
		clusCOM[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			clusCOM[ij][kl] = 0.0;

	cophDIST = (	  double **)malloc( NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++)
		cophDIST[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );

	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			cophDIST[ij][kl] = 0.0;	

	cophSIGN = (	  double **)malloc( NHELIUM*sizeof( double *) );
	for (ij = 0; ij < NHELIUM; ij++)
		cophSIGN[ij] = ( double *)malloc( NHELIUM*sizeof( double  ) );
	for (ij = 0; ij < NHELIUM; ij++)
		for (kl = 0; kl < NHELIUM; kl++)
			cophSIGN[ij][kl] = 0.0;

  clusPOPU 		= (		int	*) malloc( NHELIUM*sizeof( int ) );
  clusGB 			= (		int	*) malloc( NHELIUM*sizeof( int ) );
  clusBULK 		= (		int	*) malloc( NHELIUM*sizeof( int ) );
	clusSeq			= (		int *) malloc( NHELIUM*sizeof( int ) );
	revAtomName = (		int	*) malloc( NHELIUM*sizeof( int ) );

	for (ij = 0; ij < NHELIUM; ij++) {
		clusPOPU[ij] = 0;
		clusGB[ij] = 0;
		clusBULK[ij] = 0;
	}
	clusPOPU[0] = NHELIUM;

	atomSIZE = (		int *) malloc( NHELIUM*sizeof( int ) );
	for (ij = 0; ij < NHELIUM; ij++)
		atomSIZE[ij] = 1;


	clusSIZE = (		int *) malloc( NHELIUM*sizeof( int ) );
	for (ij = 0; ij < NHELIUM; ij++)
		clusSIZE[ij] = 1;

	entryCFG = (double *) malloc( NENTRY*sizeof( double));

	poteng	 = (double *) malloc( NHELIUM*sizeof( double));

}

void BuildAtomCoordinate(void) {
  InitiateScan();
	for (ij = 0; ij < NHELIUM; ij++) {
		ScanForCoordinate();
	}
	CloseScan();
}

void BuildDistanceMatrix(void) {
	double dist = 0.0, del[3];
  for (ij = 0; ij < NHELIUM; ij++) {
		for (kl = (ij+1); kl < NHELIUM; kl++) {
			dist = 0.0;
			for (mn = 0; mn < 3; mn++) {
				del[mn] = fabs(atomCOOR[ij][mn] - atomCOOR[kl][mn]);
				if (del[mn] > HALFLENGTH[mn])
					del[mn] = LENGTH[mn] - del[mn];
				dist += (del[mn]*del[mn]);
			}
			atomDIST[ij][kl] = dist;
			cophDIST[ij][kl] = atomDIST[ij][kl];
		}
	}			
}

void FreeMalloc(void) {
	free(atomNAME);
	free(atomSIZE);
	free(atomLEAD);
	free(atomCOOR);
	free(atomDIST);
	free(leadSORT);
	free(cophDIST);
	free(cophSIGN);
	free(clusPOPU);
	free(clusSIZE);
	free(clusGB);
	free(atomBulkOrGB);
	free(atomBulkOrGBlocal);
	free(entryCFG);
	free(poteng);
}
