#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <time.h>


#define ARRAY_SIZE 1000

void shuffle(int* solution);
int cost(const int* solution, const uint8_t* distanceMatrix);
void pairwiseExchange(int* solution);


int main(int argc, char* argv[])
{
	time_t start, end, now;
	char fileName[50];
	uint8_t* a = calloc(ARRAY_SIZE * ARRAY_SIZE, sizeof(uint8_t));
	int i, j;

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
		a[i * ARRAY_SIZE] = (uint8_t) atoi(token);
		for(j = 1; j < ARRAY_SIZE; ++j)
		{
			token = strtok(NULL, ",");
			a[i * ARRAY_SIZE + j] = (uint8_t) atoi(token);
		}
	}

	printf("Matrix has been read.");

	//Start the SA
	time(&start);
	//GENERATE INITIAL SOLUTION
	//Generate a random permutation of the cities
	int* currentSolution = calloc(ARRAY_SIZE, sizeof(int));
	for(i = 0; i < ARRAY_SIZE; ++i)
	{
		currentSolution[i] = i;
	}
	shuffle(currentSolution);
//	printf("Initial solution: ");
//	for(i = 0; i < ARRAY_SIZE; ++i)
//	{
//		printf("%d ", currentSolution[i]);
//	}
//	printf("\n");

	int* pertSolution = calloc(ARRAY_SIZE, sizeof(int));

	//Calculate the cost of the initial solution
	int currCost = cost(currentSolution, a);

	//Set the initial temperature
	double temperature = 1000;
	//Set the cooling rate
	double coolingRate = 0.003;
	//Set the number of iterations without improvement
	int iterationsWithoutImprovement = 0;
	//Set the maximum number of iterations without improvement
	int maxIterationsWithoutImprovement = 50;
	//Set how many perturbations to do when kicking out of a local minimum
	int kickNumber = 20;

	//Set the best solution
	int* bestSolution = calloc(ARRAY_SIZE, sizeof(int));
	memcpy(bestSolution, currentSolution, ARRAY_SIZE * sizeof(int));
	int bestCost = currCost;

	int stop = 0;
	//Loop until the stopping condition is met
	while(!stop)
	{
		//create perturbation
		memcpy(pertSolution, currentSolution, ARRAY_SIZE * sizeof(int));
		pairwiseExchange(pertSolution);

		//Calculate the cost of the perturbed solution
		int pertCost = cost(pertSolution, a);

		//decide weather to accept the perturbation
		//if(pertCost < currCost || (rand() / (double) RAND_MAX) < exp((currCost - pertCost) / temperature))
		if(pertCost < currCost)
		{
			//accept the perturbation
			memcpy(currentSolution, pertSolution, ARRAY_SIZE * sizeof(int));
			currCost = pertCost;
			//reset the number of iterations without improvement
			iterationsWithoutImprovement = 0;
			//check if the new solution is the best solution
			if(currCost < bestCost)
			{
				memcpy(bestSolution, currentSolution, ARRAY_SIZE * sizeof(int));
				bestCost = currCost;
			}
		}
		else
		{
			//increase the number of iterations without improvement
			iterationsWithoutImprovement++;
		}

		//cool the system
		temperature *= (1 - coolingRate);

		//kick if it is stuck
		if(iterationsWithoutImprovement > maxIterationsWithoutImprovement)
		{
			for(i = 0; i < kickNumber; ++i)
			{
				pairwiseExchange(currentSolution);
			}
		}

		//detect if the stopping condition is met
		time(&now);
		if(difftime(now,start) > 57)
		{
			stop = 1;
		}
	}
	time(&end);
	//Print the best solution
	printf("Best solution: ");
	for(i = 0; i < ARRAY_SIZE; ++i)
	{
		printf("%d ", bestSolution[i]);
	}
	printf("\n");
	printf("Best cost: %d\n", bestCost);
	printf("Time: %f\n", difftime(end, start));
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

