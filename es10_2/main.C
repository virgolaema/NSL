#include "TSP.h"
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[]){
    int n_cities;  //number of cities
    int n_indiv;  //number of individuals in each gen
    int n_gens;   //number of generations
    int disp, loss, p;

    ifstream in ("input.dat");
    in >> n_cities;
    cout << "Number of cities " << n_cities << endl;
    in >> n_indiv;
    in >> n_gens;
    cout << "Number of generations " << n_gens << endl;
    in >> disp;
    in >> loss;
    cout << "Loss function " << loss << endl;
    in >> p;
    in.close();

    //multicore
    int nstep = 2000; 
    int n_migr= 20; //nodes after which exchange is made
    int size, rank;
    int exchange, exchange2, exchange3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat[4];
    MPI_Request req, req2;
    int itag=1;
    cout << "\n RANK E' " << rank << endl << endl;
    Genetic tsp (n_cities);

    Random rnd_appo;
    int seed[4];
    int p1, p2;
    ifstream Primes("../Primes");
    if (Primes.is_open()){
       for (int r = 0; r <= rank; r++) Primes >> p1 >> p2; //differentiate the primes between ranks
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    cout << "\n generatori P1 E P2 SONO " << p1 << "," << p2 << endl << endl;

    ifstream input("../seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd_appo.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    tsp.rnd = rnd_appo;

    cout << "random: " << tsp.rnd.Rannyu() << endl;

    tsp.dispose(disp);        //dispose the n cities on a circ (0) or a square (1)
    tsp.whichLoss(loss);      //1 for L1, 2 for L2
    tsp.whichP(p);            //exp of selection operator
    tsp.whichProb(0.1);       //probabilty of mutations

    ofstream pos ("cities.out");
    for (int i = 0; i < n_cities; i++) pos << tsp.cities[i].getX() << "," << tsp.cities[i].getY() << endl;
    pos.close();

    //Create first generation
    vector <Individual> generation; 
    generation.resize(n_indiv);
    for (int i = 0; i < n_indiv; i++){
        generation[i].setN(n_cities);
        tsp.newIndividual(generation[i].sequence); //start creating random individuals
        tsp.check(generation[i].sequence);
        generation[i].fitness = tsp.loss(generation[i].sequence);
    }
    tsp.order(generation);

    int best_ind [32];

    //Create new generations
    ofstream fit ("fitness.out");
    ofstream meanL ("meanL.out");

    for (int i = 1; i <= n_gens; i++){ //should be n_gens
        if(i%10 == 1) cout << "Creating the " << i << " generation" << endl;
        //NEW GEN
        tsp.crossover(generation);
        //MUTATIONS. Every mutation also puts element in fitness order
        tsp.randomSwap(generation); 
        tsp.contSwap(generation);
        tsp.shift(generation);
        for (int y = 0; y < n_cities; y++) best_ind[y] = generation[0].sequence[y];
        if(i%10 ==0) cout << "Loss of best element at this gen: " << generation[0].fitness << endl;
        fit << generation[0].fitness << endl;
        /*double mean_L = 0.;
        for (int h = 0; h < n_indiv/2.; h++)  mean_L += generation[h].fitness;
        meanL << mean_L/(n_indiv/2.) << endl;
        */
        // EXCHANGE best individual
        if(i%n_migr==0){ 
            exchange = tsp.rnd.Rannyu(1,4);
            if(exchange == 1) {exchange2=2; exchange3=3;}
            else if(exchange == 2) {exchange2=1; exchange3=3;}
            else {exchange2=1; exchange3=2;}

            if(rank == 0){ //exchange between rank 0 and exchange
                MPI_Isend(best_ind, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &req);
                MPI_Recv(best_ind, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &stat[1]);
            }
            else if(rank == exchange){
                MPI_Send(best_ind, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
                MPI_Recv(best_ind, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &stat[0]);
            }
            else if(rank == exchange2){ //faccio scambio tra rank==exchange2 e rank==exchange3
                MPI_Isend(best_ind, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &req2);
                MPI_Recv(best_ind, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &stat[3]);
            }
            else if(rank == exchange3){
                MPI_Send(best_ind, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD);
                MPI_Recv(best_ind, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD, &stat[2]);
            }
            for (int y = 0; y < n_cities; y++) generation[0].sequence[y] = best_ind[y];
            generation[0].fitness = tsp.loss(generation[0].sequence);
            tsp.order(generation);
        }
    }
    meanL.close();
    fit.close();

    ofstream final("final.out");
    cout << "Printing best configuration...\n"; 
    for (int h = 0; h < n_cities; h++) final << tsp.cities[generation[0].sequence[h]].getX() << "," << tsp.cities[generation[0].sequence[h]].getY() << endl;
    final << tsp.cities[generation[0].sequence[0]].getX() << "," << tsp.cities[generation[0].sequence[0]].getY() << endl;
    final.close();
    cout << "Done\n";
    MPI::Finalize();
    cout << "finalized\n";
    return 0;
}