#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <mpi.h>


#define ARRAY_SIZE 1000

void shuffle(int* solution);
int cost(const int* solution, const uint8_t* distanceMatrix);
void pairwiseExchange(int* solution);


int main(int argc, char* argv[])
{
	long long int itrs = 0;
	long long int totalItrs = 0;
	int myRank, commSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	time_t start, end, now;
	uint8_t* distanceMatrix = calloc(ARRAY_SIZE * ARRAY_SIZE, sizeof(uint8_t));

	int i, j;
	if(myRank == 0)
	{
		char fileName[50];
		// Create filename
		sprintf(fileName, "../DistanceMatrix1000_v2.csv");
		printf("Reading file %s of size %dx%d\n", fileName, ARRAY_SIZE, ARRAY_SIZE);
		// Open file
		FILE* dataFile = fopen(fileName, "r");
		// Read elements
		char data[5000];
		for(i = 0; i < ARRAY_SIZE; ++i)
		{
			fgets(data, 5000, dataFile);
			const char* token = strtok(data, ",");
			distanceMatrix[i * ARRAY_SIZE] = (uint8_t) atoi(token);
			for(j = 1; j < ARRAY_SIZE; ++j)
			{
				token = strtok(NULL, ",");
				distanceMatrix[i * ARRAY_SIZE + j] = (uint8_t) atoi(token);
			}
		}

		printf("Matrix has been read.\n");
	}

	MPI_Bcast(distanceMatrix, ARRAY_SIZE * ARRAY_SIZE, MPI_UINT8_T, 0, MPI_COMM_WORLD);

	//Start the SA
	time(&start);
	//GENERATE INITIAL SOLUTION
	//Generate distanceMatrix random permutation of the cities
	int* myCurrentSolution = calloc(ARRAY_SIZE, sizeof(int));

	for(i = 0; i < ARRAY_SIZE; ++i)
	{
		myCurrentSolution[i] = i;
	}
	shuffle(myCurrentSolution);

	int* pertSolution = calloc(ARRAY_SIZE, sizeof(int));

	//Calculate the cost of the initial solution
	int myCurrentCost = cost(myCurrentSolution, distanceMatrix);

	//Set the initial temperature
	double temperature = 1000;
	//Set the cooling rate
	double coolingRate = 0.003;
	//Set the number of iterations without improvement
	int iterationsWithoutImprovement = 0;
	//Set the maximum number of iterations without improvement
	int maxIterationsWithoutImprovement = 50;
	//Set how many perturbations to do when kicking out of distanceMatrix local minimum
	int kickNumber = 20;

	//Set the best solution
	int* myBestSolution = calloc(ARRAY_SIZE, sizeof(int));
	memcpy(myBestSolution, myCurrentSolution, ARRAY_SIZE * sizeof(int));
	int myBestCost = myCurrentCost;

	int stop = 0;
	//Loop until the stopping condition is met
	while(!stop)
	{
		itrs++;
		//create perturbation
		memcpy(pertSolution, myCurrentSolution, ARRAY_SIZE * sizeof(int));
		pairwiseExchange(pertSolution);

		//Calculate the cost of the perturbed solution
		int pertCost = cost(pertSolution, distanceMatrix);

		//decide weather to accept the perturbation
		if(pertCost < myCurrentCost)
		{
			//accept the perturbation
			memcpy(myCurrentSolution, pertSolution, ARRAY_SIZE * sizeof(int));
			myCurrentCost = pertCost;
			//reset the number of iterations without improvement
			iterationsWithoutImprovement = 0;
			//check if the new solution is the best solution
			if(myCurrentCost < myBestCost)
			{
				memcpy(myBestSolution, myCurrentSolution, ARRAY_SIZE * sizeof(int));
				myBestCost = myCurrentCost;
			}
		}
		else
		{
			//increase the number of iterations without improvement
			iterationsWithoutImprovement++;
		}

		//kick if it is stuck
		if(iterationsWithoutImprovement > maxIterationsWithoutImprovement)
		{
			for(i = 0; i < kickNumber; ++i)
			{
				pairwiseExchange(myCurrentSolution);
			}
		}

		//detect if the stopping condition is met
		time(&now);
		if(difftime(now, start) > 57)
		{
			stop = 1;
		}
	}

	int globalBestCost;
	MPI_Allreduce(&myBestCost, &globalBestCost, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	if(myBestCost == globalBestCost)
	{
		if(myRank != 0)
		{
			MPI_Send(myBestSolution, ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	if(myRank == 0)
	{
		int* globalBestSolution = calloc(ARRAY_SIZE, sizeof(int));
		if(myBestCost == globalBestCost)
		{
			memcpy(globalBestSolution, myBestSolution, ARRAY_SIZE * sizeof(int));
		}
		else
		{
			MPI_Recv(globalBestSolution, ARRAY_SIZE, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		time(&end);
		//Print the best solution
		printf("Best solution: \n");
		for(i = 0; i < ARRAY_SIZE; ++i)
		{
			printf("%d ", globalBestSolution[i]);
		}
		printf("\n");
		printf("Best cost: %d\n", globalBestCost);
		printf("Time: %f\n", difftime(end, start));
	}
	MPI_Reduce(&itrs, &totalItrs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	if(myRank == 0)
	{
		printf("Total iterations: %lld\n", totalItrs);
	}

	MPI_Finalize();
}


void shuffle(int* solution)
{
	size_t i;
	for(i = 0; i < ARRAY_SIZE - 1; i++)
	{
		size_t j = i + rand() / (RAND_MAX / (ARRAY_SIZE - i) + 1);
		int t = solution[j];
		solution[j] = solution[i];
		solution[i] = t;
	}
}


int cost(const int* solution, const uint8_t* distanceMatrix)
{
	int i;
	int cost = 0;
	for(i = 0; i < ARRAY_SIZE - 1; ++i)
	{
		cost += distanceMatrix[solution[i] * ARRAY_SIZE + solution[i + 1]];
	}
	cost += distanceMatrix[solution[ARRAY_SIZE - 1] * ARRAY_SIZE + solution[0]];
	return cost;
}


void pairwiseExchange(int* solution)
{
	int i = rand() % ARRAY_SIZE;
	int j = rand() % ARRAY_SIZE;
	int temp = solution[i];
	solution[i] = solution[j];
	solution[j] = temp;
}

void cycleOfThree(int* solution)
{
	int i = rand() % ARRAY_SIZE;
	int j = rand() % ARRAY_SIZE;
	int k = rand() % ARRAY_SIZE;
	int temp = solution[i];
	solution[i] = solution[j];
	solution[j] = solution[k];
	solution[k] = temp;
}

void inversionPerturbation(int* solution)
{
	int i = rand() % ARRAY_SIZE;
	int j = rand() % ARRAY_SIZE;
	int temp;
	while(i < j)
	{
		temp = solution[i];
		solution[i] = solution[j];
		solution[j] = temp;
		i++;
		j--;
	}
}

void swapBlockPair(int* solution)
{
	int i = rand() % ARRAY_SIZE;
	int j = rand() % ARRAY_SIZE;
	int k = rand() % ARRAY_SIZE;
	int temp;
	while(i < j)
	{
		temp = solution[i];
		solution[i] = solution[k];
		solution[k] = temp;
		i++;
		k++;
	}
}