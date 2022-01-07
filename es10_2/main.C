#include "TSP.h"
#include "mpi.h"

using namespace std;

int Pbc (int max, int el){
    if(el >= max) el -= max;
    else if(el < 0) el += max;
    return el;
}

int main(int argc, char *argv[]){
    int n_cities;  //number of cities
    int n_indiv;  //number of individuals in each gen
    int n_gens;   //number of generations
    int disp, loss, p;
    double prob;

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
    cout << "Exp of selection operator: " << p << endl;
    in >> prob;
    cout << "Prob of mutations: " << prob << endl;
    in.close();

    //multicore
    int nstep = 2000; 
    int n_migr= 20; //nodes after which exchange is made
    int size, rank, max_rank = atoi(argv[1]);

    cout << "Cores used: " << max_rank << endl;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << "\n RANK is " << rank << endl << endl;
    Genetic tsp (n_cities);

    Random rnd_appo;
    int seed[4];
    int p1, p2;
    ifstream Primes("../Primes");
    if (Primes.is_open()){
       for (int r = 0; r <= rank; r++) Primes >> p1 >> p2; //differentiate the primes between ranks
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    // cout << "\n generatori P1 E P2 SONO " << p1 << "," << p2 << endl << endl;

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

    if (rank == 0) tsp.dispose(disp);        //dispose the n cities on a circ (0) or a square (1)
    tsp.whichLoss(loss);      //1 for L1, 2 for L2
    tsp.whichP(p);            //exp of selection operator
    tsp.whichProb(prob);      //probabilty of mutations

    double cities_x [n_cities] = {};
    double cities_y [n_cities] = {};

    if (rank == 0){
        ofstream pos ("cities.out");
        for (int i = 0; i < n_cities; i++){
            pos << tsp.cities[i].getX() << "," << tsp.cities[i].getY() << endl;
            cities_x[i] = tsp.cities[i].getX();
            cities_y[i] = tsp.cities[i].getY();
        }
        pos.close();
    }

    // Broadcast the disposition of cities to other ranks
    MPI_Bcast(cities_x,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);		
	MPI_Bcast(cities_y,n_cities,MPI_DOUBLE,0, MPI_COMM_WORLD);
    for (int i = 0; i < n_cities; i++){
        tsp.cities[i].setX(cities_x[i]);
        tsp.cities[i].setY(cities_y[i]);
    }

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

    int best_ind [n_cities] = {};

    //Create new generations
    ofstream fit, meanL;
    if (rank == 0){
        fit.open("fitness.out");
        meanL.open("meanL.out");
    }
     
    for (int i = 1; i <= n_gens; i++){ 
        //NEW GEN
        tsp.crossover(generation);
        //MUTATIONS. Every mutation also puts element in fitness order
        tsp.randomSwap(generation); 
        tsp.contSwap(generation);
        tsp.shift(generation);
        for (int y = 0; y < n_cities; y++) best_ind[y] = generation[0].sequence[y];
        if(i%10 ==0) cout << "RANK " << rank <<  ", Loss of best element at gen " << i << ": " << generation[0].fitness << endl;
        if (rank == 0) fit << generation[0].fitness << endl;
        double mean_L = 0.;
        if (rank == 0){
            for (int h = 0; h < n_indiv/2.; h++)  mean_L += generation[h].fitness;
            meanL << mean_L/(n_indiv/2.) << endl;
        }			
        //exchange best individuals every n_migr generations
        if (i % n_migr == 0 and max_rank != 1){
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1,req2;
			int itag[max_rank];
            for (int k = 0; k < max_rank; k++) itag[k] = k+1;
            int exchange [max_rank]; //just because can't use a vector
            if (rank == 0) {
			    vector<int> vec_exc;
                for (int k = 0; k < max_rank; k++) vec_exc.push_back(k);
                // auto rng = default_random_engine {};
                // random_shuffle(begin(vec_exc),end(vec_exc), rng);
                if(max_rank != 2) tsp.shuffle(vec_exc); 
                for (int k = 0; k < max_rank; k++) {exchange[k] = vec_exc[k]; cout << exchange[k];}
                cout << endl;
            }

            MPI_Bcast(exchange,max_rank,MPI_INTEGER,0, MPI_COMM_WORLD);
            
            //now every continent sends the best_ind to the next, and receives from the precious (in pbc)
            for (int h = 0; h < max_rank; h++){
                int this_el = h;
                int prev_el = Pbc(max_rank, h-1);
                int next_el = Pbc(max_rank,h+1);
                cout << "Elems in order: " << prev_el << " " << this_el << " " << next_el << endl;
                if (rank == exchange[this_el]){		
                    MPI_Isend(best_ind,n_cities,MPI_INTEGER,exchange[next_el],itag[this_el],MPI_COMM_WORLD,&req1);
                    MPI_Recv(best_ind,n_cities,MPI_INTEGER,exchange[prev_el],itag[prev_el], MPI_COMM_WORLD,&stat2);
                }
            }
        }
    }

    if (rank == 0){
        meanL.close();
        fit.close();
        ofstream final("final.out");
        cout << "Printing best configuration...\n"; 
        for (int h = 0; h < n_cities; h++) final << tsp.cities[generation[0].sequence[h]].getX() << "," << tsp.cities[generation[0].sequence[h]].getY() << endl;
        final << tsp.cities[generation[0].sequence[0]].getX() << "," << tsp.cities[generation[0].sequence[0]].getY() << endl;
        final.close();
    }

    cout << "RANK " << rank << " Done\n";
    MPI::Finalize();
    cout << "RANK " << rank << "finalized\n";
    return 0;
}