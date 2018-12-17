//
//  main.cpp
//  Inverse-statistical-methods-and-pseudolikelihood-approximation
//
//  Created by Samuel Bosch on 11/10/18.
//

#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

using namespace std;


int main(int argc, const char * argv[]){
    cout << "\n\n\nProgram started. This is a small version of the simulation which will be completed in less than 30s.\n\n";
    random_device rd;  // only used once to initialise (seed) engine
    mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    int min_spin = 0;  // min value spin can take
    int max_spin = 20; // max value spin can take
    int N = 25;        //Number of atoms in molecule
    double T = 1.0;  // This is the temperature of the experiment
    double K_b = 1.0; // Don't change this
    double beta = 1/(K_b*T);
    string T_str = to_string(T);
    uniform_int_distribution<int> random_spin(min_spin,max_spin); // definition of random function
    uniform_int_distribution<int> random_atom(0,N-1); // definition of random function
    // (max_spin-min_spin+1)*N
    // array<array<double, 525>, 525> J;
    double J[(max_spin-min_spin+1)*N][(max_spin-min_spin+1)*N];
    ifstream myfile;
    myfile.open("J.txt");
    if(!myfile){ //Testing if file is actually open.
        cout << "Error opening file" << endl;
        return -1;
    }else{
        cout << "File 'J.txt' is open" << "\n\n";
    }
    for (int i = 0; i < (max_spin-min_spin+1)*N; i++){
        for(int j = 0; j < (max_spin-min_spin+1)*N; ++j){
            myfile >> J[i][j]; // Here we read the file number by number
        }
    }

    // Random initialization of spins
    cout << "Initial spin configuration:\n";
    vector<int> v(N);    // declares a vector of integers
    for(int i=0; i<N; i++){
        auto random_integer = random_spin(rng);
        cout << random_integer << ' ';
        v[i] = random_integer;
    }

    // Calculation of the energy using Pott's model (H = -J*sum(Kronecker_delta(i,j))
    double E = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            E = E + J[21*i+v[i]][21*j+v[j]]; // Is this the right formula for calculating the energy?
        }
    }
    cout << "\n\nEnergy(initial)" << " = " << E << "\n";

    cout << "Starting 1st part of simulation for thermalization and Blocking method analysis...\n";
    long max_number_of_interations = 200000;
    vector<double> Energy(max_number_of_interations);
    long iteration_number = 0;
    int t = 1;
    for(; iteration_number<max_number_of_interations; iteration_number++){
        if (iteration_number % int(max_number_of_interations/10) == 0){
            cout << "Progress: " << t*10 << "%\n";
            t++;
        }
        Energy[iteration_number] = E;
        auto atom_number = random_atom(rng); //Pick random atom for changing the spin
        int old_spin = v[atom_number]; //Saving old spin in case we still want to use it
        double E_old = E;
        v[atom_number] = random_spin(rng); //Pick random new spin for selected atom
        if(v[atom_number] == old_spin){
            continue; // If the spin didn't change, we do nothing
        }
        // Calculation of the new energy (new, more efficient algorithm)
        for(int i=0; i<N; i++){
            if (atom_number==i){
                continue;
            }
            E = E + J[21*i+v[i]][21*atom_number+v[atom_number]]; // Here we add the new energy
            E = E - J[21*i+v[i]][21*atom_number+old_spin]; // Here we substract the old energy
        }

        if(E > E_old){
            double uniform_random_0_1 = (double)rand()/(double)RAND_MAX;
            double prob = exp(-(E-E_old)/(K_b*T));
            if (prob < uniform_random_0_1){
                // With some small probability, we revert the change (Metropolis algorithm)
                E = E_old;
                v[atom_number] = old_spin;
            }
        }
    }
    cout << "...done with 1st part of the simulation\n";
    cout << "Final spin configuration:\n";
    for(int i=0; i<N; i++){
        cout << v[i] << ' ';
    }
    for(int iteration_number=0; iteration_number<max_number_of_interations; iteration_number++){
        // cout << Energy[iteration_number] << '\n';
    }

    cout << endl;
    cout << "Energy(" << iteration_number << ")" << " = " << E << "\n";
    cout << '\n';

    double mean = accumulate(begin(Energy), end(Energy), 0.0) / Energy.size();
    //cout << "Average Energy = " << mean << '\n';


    // Writing the Energy vs step number into a .txt file
    // The specific path was need, as it is otherwise saved in the xcode hidden folder
    ofstream energy_file;
    energy_file.open("Energy_vs_time.txt");
    if (myfile.is_open()) { cout << "File 'Energy_vs_time.txt' is open and ready for writing...\n\n"; }
    energy_file << N << '\n';
    for(int i=0; i<max_number_of_interations; i++){
        energy_file << Energy[i] << '\n';
    }
    energy_file.close();



// Autocorrelation function
    cout << "Starting autocorrelation function calculation....\n";
    double autocorrelation_fraction = 0.03; // Through what fraction of the data do you want the autocorr. function to go?
    vector<double> autocorrelation((int)(autocorrelation_fraction*max_number_of_interations));
    t = 1;
    for(int i=0; i<int(autocorrelation_fraction*max_number_of_interations); i++){
        if (i % int(int(autocorrelation_fraction*max_number_of_interations)/10) == 0){
            cout << "Progress: " << t*10 << "%\n";
            t++;
        }
        autocorrelation[i] = 0;
        for(int j=0; j<max_number_of_interations-i; j++){
            autocorrelation[i] += Energy[j]*Energy[j+i] - mean*mean;
        }
        if (max_number_of_interations-i>0){
            autocorrelation[i] /= (max_number_of_interations-i);
        }
    }
    double normalisation_factor = autocorrelation[0];
    for(int i=0; i<int(autocorrelation_fraction*max_number_of_interations); i++){
        autocorrelation[i] /= normalisation_factor;
    }


    // Writing the autocorrelation function into a .txt file
    // The specific path was need, as it is otherwise saved in the xcode hidden folder
    int relaxation_time = -1;
    int check_var = 1;
    ofstream autocorrelation_file;
    autocorrelation_file.open("Autocorrelation.txt");
    cout << "...done with autocorrelation function calculation.\n";
    if (autocorrelation_file.is_open()) {
        cout << "File 'Autocorrelation.txt' is open and ready for writing\n\n";
    }
    for(int i=0; i<int(autocorrelation_fraction*max_number_of_interations); i++){
        autocorrelation_file << autocorrelation[i] << '\n';
        if (autocorrelation[i]<0.3679 && check_var){
            check_var = 0;
            relaxation_time = i;
        }
    }
    autocorrelation_file.close();


    // The blocking method analysis
    int n = 10; //number of blocks
    n++;
    vector<double> block_averages(n);
    vector<double> block_std(n);
    double sum = 0;
    int k = int(max_number_of_interations/n);
    int j = 0;
    for (int i=0; i<max_number_of_interations; i++){
        sum += Energy[i];
        if (i==k){
            block_averages[j] = sum/(int(max_number_of_interations/n));
            sum = 0;
            k += int(max_number_of_interations/n);
            j++;
        }
    }
    sum = 0;
    j = 0;
    k = int(max_number_of_interations/n);
    for (int i=0; i<max_number_of_interations; i++){
        sum += pow((Energy[i]-block_averages[j]),2);
        if (i==k){
            block_std[j] = sqrt(sum/(int(max_number_of_interations/n)-1));
            sum = 0;
            k += int(max_number_of_interations/n);
            j++;
        }
    }
    cout << "Blocking method analysis with " << n-1 << " blocks:\n";
    for (int i=0; i<n-1; i++){
        cout << "E[block " << i+1 << "] = " << block_averages[i] << " +/- " << block_std[i] << '\n';
    }
    
    cout << "Relaxation time = " << 1.0*relaxation_time/N << '\n';
    
    
    // (Re)calculation of the energy using Pott's model (H = -J*sum(Kronecker_delta(i,j))
    E = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            E = E + J[21*i+v[i]][21*j+v[j]];
        }
    }

    // Defining the f_{i,j}(A,B) and f_{i}(A) functions
    double f_2D[N][N][max_spin+1][max_spin+1];
    double C[N][N][max_spin+1][max_spin+1];
    for (int i=0; i<max_spin+1; i++){
        for (int j=0; j<max_spin+1; j++){
            for (int ii=0; ii<N; ii++){
                for (int jj=0; jj<N; jj++){
                    f_2D[ii][jj][i][j] = 0.0;
                    C[ii][jj][i][j] = 0.0;
                }
            }
        }
    }
    double f_1D[N][max_spin+1];
    for (int i=0; i<max_spin+1; i++){
        for (int ii=0; ii<N; ii++){
            f_1D[ii][i] = 0.0;
        }
    }

    // New round of simulations
    max_number_of_interations = 20000000;
    cout << "\n\n\nStarting main simulation with " << max_number_of_interations << " steps...\n";
    long counter = 0;
    iteration_number = 0;
    t = 1;
    for(; iteration_number<max_number_of_interations; iteration_number++){
        if (iteration_number % int(max_number_of_interations/10) == 0){
            cout << "Progress: " << t*10 << "%\n";
            t++;
        }
        auto atom_number = random_atom(rng); //Pick random atom for changing the spin
        int old_spin = v[atom_number]; //Saving old spin in case we still want to use it
        double E_old = E;
        v[atom_number] = random_spin(rng); //Pick random new spin for selected atom
        if(v[atom_number] == old_spin){
            continue; // If the spin didn't change, we do nothing
        }
        // Calculation of the new energy (new, more efficient algorithm)
        for(int i=0; i<N; i++){
            if (atom_number==i){
                continue;
            }
            E += J[21*i+v[i]][21*atom_number+v[atom_number]]; // Here we add the new energy
            E -= J[21*i+v[i]][21*atom_number+old_spin]; // Here we substract the old energy
        }

        if(E > E_old){
            double uniform_random_0_1 = (double)rand()/(double)RAND_MAX;
            double prob = exp(-(E-E_old)/(K_b*T));
            if (prob < uniform_random_0_1){
                // With some small probability, we revert the change (Metropolis algorithm)
                E = E_old;
                v[atom_number] = old_spin;
            }
        }
        if (iteration_number % (3*relaxation_time) == 0){ // We wait for 3 times the relaxation time before we take another measurement
            counter++;
            for (int j=0; j<N; j++){
                f_1D[j][v[j]]++;
            }
            for (int j=0; j<N; j++){
                for (int i=0; i<N; i++){
                    f_2D[j][i][v[j]][v[i]]++;
                }
            }
        }
    }


    // converting f_1D and f_2D to probabilities
    for (int i=0; i<max_spin+1; i++){
        for (int j=0; j<max_spin+1; j++){
            for (int ii=0; ii<N; ii++){
                for (int jj=0; jj<N; jj++){
                    f_2D[ii][jj][i][j] /= counter;
                    // cout << f_2D[ii][jj][i][j] << '\n';
                }
            }
        }
    }
    for (int i=0; i<max_spin+1; i++){
        for (int j=0; j<N; j++){
            f_1D[j][i] /= counter;
            // cout << f_1D[j][i] << '\n';
        }
    }

    for (int i=0; i<max_spin+1; i++){
        for (int j=0; j<max_spin+1; j++){
            for (int ii=0; ii<N; ii++){
                for (int jj=0; jj<N; jj++){
                    C[ii][jj][i][j] = f_2D[ii][jj][i][j] - f_1D[ii][i]*f_1D[jj][j];
                    // cout << C[ii][jj][i][j] << '\n';
                }
            }
        }
    }


    // Writing the autocorrelation function into a .txt file
    // The specific path was need, as it is otherwise saved in the xcode hidden folder
    cout << "...main simulation successfully completed.\n\nPreparing for writing data to file.\n";
    ofstream C_file;
    string s1 = "C_T_";
    string s2 = "_counter_";
    string s3 = ".txt";
    string counter_str = to_string(counter);
    string filename = s1 + T_str + s2 + counter_str + s3;
    // For practicality, in this small simulation, we are just going to use C.txt
    filename = "C_test_simulation.txt";
    cout << filename << "\n";
    C_file.open(filename);
    if (C_file.is_open()){
        cout << "File 'C.txt' is open\n\n";
        C_file << max_spin << '\n';
        C_file << N << '\n';
        C_file << beta << '\n';
        C_file << counter << '\n';
        C_file << 1.0*relaxation_time/N << '\n';
        C_file << max_number_of_interations << '\n';
        for (int i=0; i<max_spin+1; i++){
            for (int j=0; j<max_spin+1; j++){
                for (int ii=0; ii<N; ii++){
                    for (int jj=0; jj<N; jj++){
                        C_file << C[ii][jj][i][j] << '\n';
                    }
                }
            }
        }
    }
    C_file.close();
    cout << "Program completed successfully. All data is written to files.\n\n\n\n";

    return 0;
}



