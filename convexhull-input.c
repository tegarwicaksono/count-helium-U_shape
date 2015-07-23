  char nameFeInput[] = "./Filter-Ovito/CFG.OVITO.CURV.GB.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.cfg";
  char charFeInput[sizeof nameFeInput+200];

  char nameFeCand[]  = "./Candidate-Fe/IronCandidate.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.txt";
  char charFeCand[sizeof nameFeCand+200];

  char nameFeVert[]  = "./Vertex-Fe/ConvexHullVertex.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.txt";
  char charFeVert[sizeof nameFeVert+200];

  char nameHeInput[] = "./CFG-He/CFG.CURV.HELIUMONLY.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.cfg";
  char charHeInput[sizeof nameHeInput+200];

  char nameHeEvol[] = "./He-evolution/Evol.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvol[sizeof nameHeEvol+200];

  char skip;

	int NHELIUM = 500; //In absolute number
	int	NINIT		= 5;
  int TEMP		= 1000;
  int SIGMA   = 3;
	int NIRON		= 1;
	int NLCAND	= 0;
	int NRCAND  = 0;
	int NMCAND	= 0;
	int	NCAND		= 0;
	int NVERTEX = 1;
  int LINEFE  = 15;
	int LINEHE	= 20;
	int ENTRYFE = 0;
	int ENTRYHE = 5;
	int NENTRY  = 1;
	int NSIZEPRINT = 15;


  int inistep  =   100000;
  int finstep  =  5000000;
	int intstep	 =   100000;
	int	NPAIRDIST, NCLUSTERS;
	int	NHELIUMBULK, NHELIUMBULKOUT, NHELIUMEDGE, NHELIUMLOOP, NHELIUMMINUSONE;

  int    	step, atom, timestep, c, lines = 0;
	int idMAXLEFT, idMINLEFT, idMAXRIGHT, idMINRIGHT;

	double	int1lo = 0.31;
	double	int1hi = 0.19;
	double	int2lo = 0.81;
	double	int2hi = 0.69;
	double	cutoffPOTEN = 1.00;

//	int			NATOMS, NATOMSM;

  double 	x, y, z, ix, iy, iz;
  double 	lx = 0.0, ly = 0.0, lz = 0.0;

	double	CUTOFFDISTANCE = 7.5;	//in Angstrom
	double 	CorrCoefficient;
	double 	CUTOFFRAD = 150.0; //In Angstrom
	double 	CUTOFFRADSQ;

	double	TOPLEFT[3], TOPRIGHT[3], BOTLEFT[3], BOTRIGHT[3];
	double	LEFTEDGE[2] = {0.138, 0.152};		//For the left edge of U-loop
	double	RIGHTEDGE[2] = {0.850, 0.864};	//For the right edge of U-loop
	double	LEFTLOC  = 0.140;
	double	RIGHTLOC = 0.857;
	double	upperBorder  = 1.0;

  double 	LENGTH[3], HALFLENGTH[3], ATOMCOOR[3], ACTCOOR[3], ENTRYCFG;

	int 		target1, target2;

	int 		*atomNAME, *atomLEAD, *atomSIZE, *leadSORT, *revAtomName;
	int 		*labelBULK, *labelEDGE, *labelLOOP, *labelBULKOUT;
	int 		*clusPOPU, *clusSIZE,  *clusSeq;
	double 	**atomCOOR, **atomDIST, **cophDIST, **cophSIGN, **clusCOM;
	double	*entryCFG, *poteng;



	double	*LEFTCAND, *RIGHTCAND, **LOOPCAND, *ZLEFT, *ZRIGHT, **VXCAND;
	double	**VXPOINT, **HULLEDGE, **distHelium, **ClosestEdge;
	int			*LabelHeAtoms;

  char 		xa, ya, za;
  char 		ig1[7], ig2[1], ig3[1];
	char		w1[6], w2[2], w3[9], w4[1];
  int  		ch, ij, kl, mn;

  char nameHeEvolSizeBulkInside[] = "./He-evolution/BulkInside/Evol.BulkIn.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeBulkInside[sizeof nameHeEvolSizeBulkInside+500];

  char nameHeEvolSizeBulkOutside[] = "./He-evolution/BulkOutside/Evol.BulkOut.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeBulkOutside[sizeof nameHeEvolSizeBulkOutside+500];

  char nameHeEvolSizeEdge[] = "./He-evolution/Edge/Evol.Edge.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeEdge[sizeof nameHeEvolSizeEdge+500];

  char nameHeEvolSizeLoop[] = "./He-evolution/Loop/Evol.Loop.COH.T%dK.NHE.ALL%d.INIT%d.time_%.1lfns_to_%.1lf_ns.txt";
  char charHeEvolSizeLoop[sizeof nameHeEvolSizeLoop+500];


	FILE *input, *validity, *cfgnew;
	FILE *ironinput, *ironCand, *heliuminput, *ironvert, *evol, *heliumcfg;
	FILE *evolsizebulkin, *evolsizebulkout, *evolsizeedge, *evolsizeloop;


void FindInit(void) {
	if (NHELIUM == 100) NINIT = 31;
	if (NHELIUM == 200) NINIT = 21;
	if (NHELIUM == 300) NINIT = 1;
	if (NHELIUM == 400) NINIT = 4;
	if (NHELIUM == 500) NINIT = 5;
	if (NHELIUM == 600) NINIT = 6;
	if (NHELIUM == 700) NINIT = 7;
	if (NHELIUM == 1000) NINIT = 0;
	if (NHELIUM == 1500) NINIT = 0;
	if (NHELIUM == 2000) NINIT = 0;
}

void AskForInput(void) {
	double ini, fin;
	printf("Enter NHELIUM: ");
	scanf("%d", &NHELIUM);
	FindInit();
	printf("Enter the timestep of the first index (in ns): ");
	scanf("%lf", &ini);
	ini = ini*1.0e6;
	inistep = (int)ini;
	
	printf("Enter the timestep of the last index (in ns): ");
	scanf("%lf", &fin);
	fin = fin*1.0e6;
	finstep = (int)fin;


}
