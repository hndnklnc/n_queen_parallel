#include "Genetic.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <omp.h>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;



struct Choromosome{
	float fitness;
	int* geneCode;
};

Genetic genetic;
const int TOURNAMENT_SIZE = 5;
const float MUTATION_RATE = 0.2;
const int NEW_POPULATION_RATE = 0.5;
const int START_SIZE = 250;
int popSIZE = 0; //population size
std::vector<Choromosome> population;
int SIZE;
int MAX_POPULATION;
int tour = 0;
int* bestChoro;
float bestFitness;
double old_clock;
double tour_begin = 0;
double tour_end = 0;
float total_tour_time = 0;

Genetic::Genetic(){
}
Genetic::Genetic(int size){
    SIZE = size;
    MAX_POPULATION =  SIZE + SIZE * NEW_POPULATION_RATE;
    bestChoro = new int[SIZE];
}

int* swapArray(int indx1, int indx2, int* array){
	int * code = new int[SIZE];

	std::swap(array[indx1], array[indx2]);
	for(int i = 0; i < SIZE; i++){
		code[i] = array[i];
	}
	return code;
}

bool controlChoromosome(int* gene){
	bool isTrue = true;
	for(int i = 0; i < SIZE; i++){
		if(gene[i] > SIZE - 1){
			isTrue = false;
		}
	}
	return isTrue;
}
void copyArray(int * a1, int * a2){

		for(int i = 0; i < SIZE; i++){
			a1[i] = a2[i];
		}


}
//Chooses 7 random chromosome and select two best fitness choromosome
int* makeTournament(){
	int *matesIndex = new int[2];
	int temp, best,second;
	best = rand() & (popSIZE - 1);
	second = rand() & (popSIZE - 1);
	if(population[best].fitness < population[second].fitness){
		temp = best;
		best = second;
		best = temp;
	}
	for(int i = 0; i < TOURNAMENT_SIZE; i++){
		temp = rand()% (popSIZE - 1);
		if(population[temp].fitness < population[best].fitness){
			best = temp;
		}else if(population[temp].fitness < population[second].fitness){
			second = temp;
		}
	}
	matesIndex[0] = best;
	matesIndex[1] = second;
	return matesIndex;
}


void clearPopulation(){
	popSIZE = 0;
	for(int i = 0; i < population.size(); i++){
		delete population[i].geneCode;
	}
	population.clear();
}


//finds the difference in two arrays
vector<int> findDifference(int* piece1, int* piece2, int pieceSize){
	vector <int> difElements;
	bool isSame = false;

		for(int i = 0; i < pieceSize; i++){
			isSame = false;
			for(int j = 0; j < pieceSize; j++){
				if(piece1[i] == piece2[j]){
					isSame = true;
				}
			}
			if(!isSame){
				difElements.push_back(piece1[i]);
			}
		}

	return difElements;
}
int findPosition(int i1, int i2, int* array, int elem){
	int index = -1;

	for(int i = i1; i < i2; i++){
		if(array[i] == elem){
			index = i;
			break;
		}
	}
	return index;
}
void nQueenControl(Choromosome c){
	if(c.fitness == 0){
		double current_clock = omp_get_wtime();
		cout<<"!!!!!!!!!!! n queen is found !!!!!!!!"<<endl;
		std::ostream_iterator< int > output( cout, " " );
		std::copy( c.geneCode, c.geneCode+ SIZE, output );
		cout<<tour<<". tur"<<endl;
		cout<<"elapsed time: "<<current_clock - old_clock<<endl;
		if(tour > 1){
			cout<<"tour time:"<< total_tour_time / float(tour - 1)<<endl;
		}else{
			cout<<"tour didnt finish"<<endl;
		}
		exit(EXIT_SUCCESS);
	}
}

float findFitness(int* code){
	int diagSize = 2 * SIZE;
	int *rightDiagonal = new int[diagSize];
	int *leftDiagonal = new int[diagSize];
	for(int i = 0; i < diagSize; i++){
		rightDiagonal[i] = 0;
		leftDiagonal[i] = 0;
	}
	for(int i = 0; i < SIZE; i++){
		rightDiagonal[i + code[i]]++;
		leftDiagonal[SIZE - i + code[i]]++;
	}
	float sum = 0;
	int counter;

		for(int j = 0; j < diagSize; j++){
			counter = 0;
			if(leftDiagonal[j] >1)
				counter  += leftDiagonal[j] - 1;
			if(rightDiagonal[j] >1)
				counter  += rightDiagonal[j] - 1;

			sum += (float) counter / (float) ((SIZE - abs(SIZE - j)) + 1); //normalization with size of corresponding diagonal

		}
		delete rightDiagonal;
		delete leftDiagonal;

	return sum;
}
void addNewIndividual(Choromosome indv){
	if(!controlChoromosome(indv.geneCode)){
		cout<<"chromosome hatasi "<<endl;
		std::ostream_iterator< int > output( cout, "_" );
		std::copy( indv.geneCode, indv.geneCode+ SIZE, output );
		cout<<"\n tour: "<<tour;
		exit(EXIT_FAILURE);
	}
	float fitness = findFitness(indv.geneCode);
	indv.fitness = fitness;
	if(bestFitness > fitness){
		copyArray(bestChoro, indv.geneCode);
		bestFitness = fitness;
	}
	nQueenControl(indv);
	population.push_back(indv);
	popSIZE++;
}
void repairChoromosome(vector<int> sameElements, Choromosome parent, Choromosome child, int index1, int index2, vector<int> exElements){

	int parentIndex;
	int exIndex = 0;
	for(int i = 0; i < sameElements.size(); i++){
		parentIndex =findPosition(0, index1, parent.geneCode, sameElements[i]);
		if(parentIndex == -1)
			parentIndex = findPosition(index2, SIZE, parent.geneCode, sameElements[i]);
		child.geneCode[parentIndex] = exElements[exIndex];
		exIndex++;
	}
	addNewIndividual(child);
}



void Crossover(Choromosome parent1, Choromosome parent2){
	Choromosome child1, child2;
	child1.geneCode = new int[SIZE];
	child2.geneCode = new int[SIZE];
	int randpos1 = rand() % SIZE;
	int randpos2 = randpos1;

	while(randpos1 == randpos2){
		randpos2 = rand() % SIZE;
	}
	if(randpos1 > randpos2){
		int tempPos = randpos1;
		randpos1 = randpos2;
		randpos2 = tempPos;
	}
	int pieceSize = randpos2 - randpos1;
	int* pieceCode2 = new int[pieceSize];
	int* pieceCode1 = new int[pieceSize];
	int l = 0;
	for(int i = randpos1; i <randpos2; i++){
		pieceCode1[l] = parent1.geneCode[i];
		pieceCode2[l] = parent2.geneCode[i];
		child1.geneCode[i] = parent2.geneCode[i];
		child2.geneCode[i] = parent1.geneCode[i];
		l++;
	}
	for(int j = 0; j < randpos1; j++){
		child1.geneCode[j] = parent1.geneCode[j];
		child2.geneCode[j] = parent2.geneCode[j];
	}

	for(int k = randpos2; k < SIZE; k++){
		child1.geneCode[k] = parent1.geneCode[k];
		child2.geneCode[k] = parent2.geneCode[k];
	}

	vector<int> difP2C1Elements = findDifference(pieceCode2, pieceCode1, pieceSize);
	vector<int> difP1C2Elements = findDifference(pieceCode1, pieceCode2, pieceSize);
	repairChoromosome(difP2C1Elements, parent1, child1, randpos1, randpos2, difP1C2Elements);
	repairChoromosome(difP1C2Elements, parent2, child2, randpos1, randpos2, difP2C1Elements);
	delete pieceCode1;
	delete pieceCode2;
	difP2C1Elements.clear();
	difP1C2Elements.clear();

}




void mutation(Choromosome chro){
	int rand1 = rand() % (SIZE - 1);
	int rand2 = rand() & (SIZE - 1);

	std::swap(chro.geneCode[rand1], chro.geneCode[rand2]);
	int fitness = chro.fitness;
	int mutFitness = findFitness(chro.geneCode);
	if(mutFitness > fitness){
		std::swap(chro.geneCode[rand1], chro.geneCode[rand2]);
		chro.fitness = findFitness(chro.geneCode);
	}
}
void createNextGeneration(){
	int randIndx = 0;
	float mutRound ;
	int count = 0;
	while(popSIZE < MAX_POPULATION){
		int* mates = makeTournament();
		for(int i = 0 ; i < SIZE; i++){
			if((population[mates[0]].geneCode[i] > SIZE - 1) || (population[mates[1]].geneCode[i] > SIZE - 1))
				exit(EXIT_FAILURE);
		}
		Crossover(population[mates[0]], population[mates[1]]);
		if(count % 20 == 0){
             mutRound = MUTATION_RATE * popSIZE;
            for(int i = 0; i < (int) mutRound; i++){
                randIndx = rand() % (popSIZE - 1);
                mutation(population[randIndx]);
            }
        }
        count++;
	}
}


void createMoreIndv(int* randomCode){
	Choromosome c;

	for(int i = 0 ; i < SIZE; i++){
		if((randomCode[i] > SIZE - 1) )
			exit(EXIT_FAILURE);
	}
	int randomIndex1 = rand() % (SIZE );
	int randomIndex2 = randomIndex1;
	while(randomIndex1 == randomIndex2){
		randomIndex2 = rand() % (SIZE );
	}
	c.geneCode = swapArray(randomIndex1, randomIndex2, randomCode);
	addNewIndividual(c);
}

void makeNewInitialization(){
	tour_begin = omp_get_wtime() ;
	clearPopulation();
	int * randCode = new int[SIZE];
	copyArray(randCode, bestChoro);
	tour++;
	for(int i = 0 ; i < SIZE; i++){
		if((bestChoro[i] > SIZE - 1) ){
			cout<<"best chro error";
			exit(EXIT_FAILURE);
		}
	}
	for(int i = 1; i < START_SIZE; i++){
		createMoreIndv(randCode);
	}
	createNextGeneration();
	tour_end = omp_get_wtime();
	total_tour_time += (tour_end - tour_begin);
	if(tour % 500 == 0){
	    cout<<"tour time:"<< total_tour_time / float( tour)<<" sec"<<endl;
	}

	makeNewInitialization();
}
void Genetic:: initializePopulation(){
	old_clock = omp_get_wtime();
	srand(time(NULL));
	int* randomCode = new int[SIZE];
	int SIZEPlusOne = SIZE + 1; //fill randomCode with it to understand that all idexes of it has different value
	bool isFilled = false; //to control idexes of randoCode are filled
	for(int j = 0; j < SIZE; j++){
		randomCode[j] = SIZEPlusOne;

	}
	for(int i = 0; i < SIZE; i++){
		int randomGene = rand() % SIZE;
		isFilled = false;
		while(!isFilled){
			if(randomCode[randomGene] == SIZEPlusOne){
				randomCode[randomGene] = i;
				isFilled = true;
			}else{
				randomGene = rand() % SIZE;
				isFilled =false;
			}
		}
	}
	copyArray(bestChoro, randomCode);
	bestFitness = findFitness(randomCode);
	Choromosome chro;
	chro.geneCode = randomCode;
	chro.fitness = bestFitness;
	addNewIndividual(chro);
	for(int k = 1; k < START_SIZE ; k++){
		createMoreIndv(randomCode);
	}
	createNextGeneration();
	makeNewInitialization();

}
