#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "convexhull-input.c"
#include "convexhull-iron.c"
#include "convexhull-convexhull.c"
#include "convexhull-cluster-buildarray.c"
#include "convexhull-cluster-readatoms.c"
#include "convexhull-cluster-performanalysis.c"
#include "convexhull-cluster-performvalidity.c"
#include "convexhull-helium.c"


int main(void) {
	AskForInput();
	StartEvolution();

  for (timestep = inistep; timestep <= finstep; timestep += intstep) {
		printf("timestep = %d\n", timestep);
		//1.Find points to construct Convex Hull
		OpenFeAtomicFile();
		BuildArrayFeAtoms();
		FindVertexCandidate();
		PrintVertexCandidate();	
//		PrintCFGCandidate();
		printf("\tStep1: Find candidate for convex hull done\n");
		FreeArrayVertex();

		//2.Perform Convex Hull
		ComputeConvexHull();
		FindVertexConvexHull();
		printf("\tStep2: Identify vertices for convex hull done\n");


		//3.Characterize Helium clusters 
		OpenHeAtomicFile();
	  BuildArrayClusters();
		ExtractHeliumCoordinate();
		BuildDistanceMatrix();
		PerformAnalysis();
		PerformAgglomerateCluster();
		PrintAtomLead();
		PrintCFGClusterCentreOfMass();
		printf("\tStep3: Characterize clusters done\n");

		//4.Calculate the distance of He atoms to the convex hull
		BuildArrayHeAtoms();

		CalcDistHeAtoms();
		BinHeAtoms();	
		PrintClusterPopulation();
		PrintCFGBinHelium();
		PrintHeliumDistance();
		PrintEvolution();

		FreeMalloc();
		printf("\tStep4: Categorize clusters done\n");

	}

	CloseEvolution();
	return 0;
}

