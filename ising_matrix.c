//sample C file to add 2 numbers - int and floats

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

// static variables
static bool  *spin;
static short *pairs;  // for each site, how many pairs anti-align with neighbors
static int   N;       // size of array
static int   NN;      // = N*N; number of elements in the array
static int   Npairs;  // sum of pairs
static int   Nspin;   // sum of spin
static float J = 1.0; // pair interactions strength
static float flip_prop = 0.1; // percent of sites allowed to flip
static float E;
static float M;
static int   HALF_RAND_MAX;
static int   FLIP_PROP_RAND_MAX;

// function signatures
int   allocate (int matrix_size);
float get_E   ();
float get_M ();
int   set_flip_prop(float);
int   set_J (float);
int   print_spins();
int   print_E();
float step(float T, float B); // returns E
int   free_mem ();
int   rand_spins();
int   get_spin(int i, int j);
int   get_N();


// function definitions
int get_spin(int i, int j){
    return (spin[i*N+j]) ? 1 : -1;
}
int free_mem() { free(spin); free(pairs); return 0; }
float get_E() { return E; }
float get_M() { return M; }
int   set_flip_prop(float flip_prop_in){ 
    flip_prop = flip_prop_in;  return 0; 
    FLIP_PROP_RAND_MAX = (int)(flip_prop * RAND_MAX);
}

int   set_J(float J_in) { J = J_in; return 0; }
int   print_spins() {
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){ printf(" %i", spin[i*N+j]); }
        printf("\n");
    }
    return 0;
}
int print_E(){ printf("\nE: %f\n",E); return 0;}
int get_N(){ return N; }

int allocate(int matrix_size){
    srand(time(NULL));
    HALF_RAND_MAX = RAND_MAX / 2;
    FLIP_PROP_RAND_MAX = (int)(flip_prop * RAND_MAX);
    N = matrix_size;
    NN = matrix_size * matrix_size;
    spin  = malloc( NN * sizeof(bool) );
    pairs = malloc( NN * sizeof(short) );
    
    // initialize spin matrix
    for (int i = 0; i < NN; ++i) { spin[i] = rand() > HALF_RAND_MAX; }

    return 0;
}


int rand_spins(){
    for (int i = 0; i < NN; ++i) { spin[i] = rand() > HALF_RAND_MAX; }
    return 0;
}

float step(float T, float B){
    //#Calculating the total spin of neighbouring cells (in this case,
    // the number of pairs pointing opposite).
    // Also sum the number pointed "up" in the direction of the magnetic field.

    //   look through the cells pairing each with the one to the right, and below
    //   it, recording the result for both.
    //
    //    0     1       2     ....  j  ...   N - 1
    //    2*N   2*N+1   ....  2*N+j    ...  2*N-1
    //    .
    //    .
    //    .
    //    i*N   ...      i*N+j         ... (i+1)*N-1
    //    .
    //    .
    //    .
    //    (N-1)*N        ...               N*N-1
    //

    //zero out pairs
    for (int i = 0; i < NN; ++i){ pairs[i] = 0; }

    for (int i = 0; i<N-1; ++i){
        // get the first and last in the row
        for (int j = 0; j<N-1; ++j){
            if (spin[i*N+j]^spin[i*N+j+1]){ // O -> O
                pairs[i*N+j] += 1;
                pairs[i*N+j+1] += 1;
            }
            if (spin[i*N+j]^spin[(i+1)*N+j]){ // O
               pairs[i*N+j] += 1;             // |
               pairs[(i+1)*N+j] += 1;         // O
            }
        }
        if (spin[(i+1)*N-1]^spin[i*N]){ // last to first column
            pairs[(i+1)*N-1] += 1;
            pairs[i*N]       += 1;
        }
        if (spin[(i+1)*N-1]^spin[(i+2)*N-1]){
           pairs[(i+1)*N-1] += 1;
           pairs[(i+2)*N-1] += 1;
        }
    }
    // Same math, but for the final row (i = N -1), and the 
    // for (int i = N-1; i < N; ++i){ // <- as if this conditon was true
    int final_row = NN-N;
    for (int j = 0; j < N - 1; ++j){
        if (spin[final_row+j] ^ spin[j]){
            pairs[final_row+j] += 1;
            pairs[j]           += 1;
        }
        if (spin[final_row+j]^spin[final_row+j+1]){
            pairs[final_row+j]   += 1;
            pairs[final_row+j+1] += 1;
        }
    }
    if (spin[NN-1]^spin[N-1]){
        pairs[NN-1] += 1;
        pairs[N-1]  += 1;
    }
    if (spin[NN-1]^spin[final_row]){
        pairs[NN-1]      += 1;
        pairs[final_row] += 1;
    }

    // #Sum up our variables of interest, normalize by N^2
    Nspin = 0;
    for (int i = 0; i < NN; ++i) Nspin += spin[i];
    M = (float)(2*Nspin - NN)/NN; // this value is now available to get taken.

    Npairs = 0;
    for (int i = 0; i < NN; ++i) Npairs += pairs[i];
    // from above, Npairs = Npairs of opposite spin
    //  -spin*neighbors  =    Npairs_opposite - Npairs_same;
    //                   =    N_opp - N_same;   
    // N_opp + N_same = 4*NN
    // Therefore, N_opp - N_same = N_opp - (4*NN - N_opp) = 2*N_opp - 4*NN
    // Divide by two to account for double counting, 
    // therefore: -spin*neighbors/2 = N_opp - 2*NN;
    // therefore: energy_spin = J*(N_opp-2*NN)
    E = (J*(Npairs-2*NN)-B*(2*Nspin-NN))/(float)NN;
    // can now get E with get_E

    // now allow sites to fip
    for (int i = 0; i < NN; ++i){
        /* if ((float)rand()/(float)RAND_MAX > flip_prop) continue; */
        if (rand() > FLIP_PROP_RAND_MAX) continue;
        float DeltaE = 2.0*(J*(4.0 - 2.0*pairs[i]) + B*(2.0*spin[i]-1.0));
        if ( DeltaE < 0.0 || exp(-1.0*DeltaE/T) > (float)rand()/(float)RAND_MAX ) spin[i] = !spin[i];
    }
    /* print_spins(); */

    return E;
}
