void ComputeConvexHull(void) {
	int i;
	FILE *pipe = popen("python","w");

	fprintf(pipe, "from convexhull_python import convex_hull\n");
	fprintf(pipe, "from numpy import *\n");
	fprintf(pipe, "data = [");
		for (i = 0; i < NCAND-1; i++) {
			fprintf(pipe, "(%.4lf,%.4lf),", VXCAND[i][0], VXCAND[i][1]);
		}
		fprintf(pipe, "(%.4lf,%.4lf)]\n", VXCAND[i][0], VXCAND[i][1]);

	fprintf(pipe, "res  = convex_hull(data)\n");
	fprintf(pipe, "savetxt(\"./Vertex-Fe/ConvexHullVertex.COH.T%dK.DIA50.NHE.ALL%d.INIT%d.%d.txt\",res,fmt=\'%%.4e\',delimiter=\'\\t\')\n",TEMP,NHELIUM,NINIT,timestep);

  pclose(pipe);
}

void FindNumberVertex(void) {
	int i;
	char ch;

	NVERTEX = 0;
	sprintf(charFeVert,nameFeVert,TEMP,NHELIUM,NINIT,timestep);
  ironvert = fopen(charFeVert,"r");

	while(!feof(ironvert))	{
  	ch = fgetc(ironvert);
  	if(ch == '\n')
  	{
   		NVERTEX++;
  	}
	}

	fclose(ironvert);

	VXPOINT  = (		double **)malloc( (NVERTEX+1)*sizeof( double *) );
	for (i = 0; i <= NVERTEX; i++) {
		VXPOINT[i] =  ( double *)malloc( 2*sizeof( double  ) );
	}

	HULLEDGE  = (		double **)malloc( NVERTEX*sizeof( double *) );
	for (i = 0; i < NVERTEX; i++) {
		HULLEDGE[i] =  ( double *)malloc( 3*sizeof( double  ) );
	}

}


void FindVertexConvexHull(void) {
	int i;
	double label, slope, intercept, dx, dy;
	FindNumberVertex();
	sprintf(charFeVert,nameFeVert,TEMP,NHELIUM,NINIT,timestep);
  ironvert = fopen(charFeVert,"r");

	for (i = 0; i < NVERTEX; i++) {
		fscanf(ironvert, "%lf %lf", &VXPOINT[i][0], &VXPOINT[i][1]);
	}

	VXPOINT[NVERTEX][0] = VXPOINT[0][0];
	VXPOINT[NVERTEX][1] = VXPOINT[0][1];

	for (i = 0; i < NVERTEX; i++) {
		dx = VXPOINT[(i+1)%NVERTEX][0] - VXPOINT[(i)%NVERTEX][0];
		dy = VXPOINT[(i+1)%NVERTEX][1] - VXPOINT[(i)%NVERTEX][1];
		if (dx == 0.0) {		//Vertical line
			slope = 0.0;
			label = 1.0;
			intercept = VXPOINT[(i+1)%NVERTEX][0];
		} else if (dy == 0.0) {
			slope = 0.0;
			label = -1.0;
			intercept = VXPOINT[(i)%NVERTEX][1];
		} else {
			slope = dy/dx;
			label = 0.0;
			intercept = VXPOINT[(i)%NVERTEX][1] - slope*VXPOINT[(i)%NVERTEX][0];
		}
		HULLEDGE[i][0] = label;
		HULLEDGE[i][1] = slope;
		HULLEDGE[i][2] = intercept;
	}

	fclose(ironvert);

	sprintf(charFeVert,nameFeVert,TEMP,NHELIUM,NINIT,timestep);
  ironvert = fopen(charFeVert,"w");
	fprintf(ironvert,"Indx\tX\tY\tLabel\tSlope\tIntercept\n");
	for (i = 0; i < NVERTEX; i++) {
		fprintf(ironvert, "%d\t",i+1);
		fprintf(ironvert, "%.6lf\t", VXPOINT[i][0]);
		fprintf(ironvert, "%.6lf\t", VXPOINT[i][1]);
		fprintf(ironvert, "%.1lf\t", HULLEDGE[i][0]);
		fprintf(ironvert, "%.6lf\t", HULLEDGE[i][1]);
		fprintf(ironvert, "%.6lf\t", HULLEDGE[i][2]);
		fprintf(ironvert, "\n", HULLEDGE[i][2]);
	}
	fclose(ironvert);

}

	

