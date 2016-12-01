#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    I decided to load all the values except for N as doubles since, as much as some of them might intuitively
    be integers (initial temperature for example), it is better to load an int as double than a double as int.
    The only exception being N which must be an int for it to make sense.
*/


struct Var {
    int N;
    double x_L;
    double x_R;
    double t_f;
    double t_d;
    double X_s;
    double X_p;
    double T_c;
    double T_w;
    double T_0;
    double T_1;
    double x_T;
    double x_w;
    double P;
};
typedef struct Var Var;


/*
    Loads variables form "input.txt" and saves them in @param v, struct which contains all inputed double values
    (all except N, number of grid points)

    @param v: Var struct
*/
void load_struct( Var *v ){
    int i;
    double temp;
    FILE *input = fopen("input.txt", "r");

    for(i=0; i<14; i++){
        fscanf( input, "%lf", &temp );
        switch( i ){
            case 0: v->x_L  = temp; break;
            case 1: v->x_R  = temp; break;
            case 2: v->N    = (int)temp; break;
            case 3: v->t_f  = temp; break;
            case 4: v->t_d  = temp; break;
            case 5: v->X_s  = temp; break;
            case 6: v->X_p  = temp; break;
            case 7: v->T_c  = temp; break;
            case 8: v->T_w  = temp; break;
            case 9: v->T_0  = temp; break;
            case 10: v->T_1 = temp; break;
            case 11: v->x_T = temp; break;
            case 12: v->x_w = temp; break;
            case 13: v->P   = temp; break;

            default: exit(1); break;
        }
    }
}

/*
    Function used in testing to make sure all variables were
    correctly loaded
*/
void print_struct( Var *v ){
    int i;
    double temp;

    for(i=0; i<14; i++){
        switch( i ){
            case 0: temp = v->x_L; break;
            case 1: temp = v->x_R; break;
            case 2: temp = (double)v->N; break;
            case 3: temp = v->t_f; break;
            case 4: temp = v->t_d; break;
            case 5: temp = v->X_s; break;
            case 6: temp = v->X_p; break;
            case 7: temp = v->T_c;  break;
            case 8: temp = v->T_w;  break;
            case 9: temp = v->T_0; break;
            case 10: temp = v->T_1 ; break;
            case 11: temp = v->x_T; break;
            case 12: temp = v->x_w; break;
            case 13: temp = v->P;   break;

            default: exit(1); break;
        }
        printf("%g ", temp);
    }
    printf("\n");
}

double thermal_conductivity( double T, Var *v ){

    double X = 0.5*(v->X_p - v->X_s);
    double DT = T - v->T_c;
    X *= (1 + tanh( DT / v->T_w ) );

    return v->X_s + X;
}

double X_forward( double T0, double T1, Var *v ){
    double X1 = thermal_conductivity( T0, v );
    double X2 = thermal_conductivity( T1, v );

    return 0.5*(X2+X1);
}


/*
    Initialises the two grid arrays, T and T_next, including allocation
    Note: I decided to ignore the bdd condition for the initialisation

    @param N:           number of grid points
    @param dx:          delta x between grid points
    @param v_address:   address of array of variables (so that the stack of the function is minimal, instead of "sending" the 14 vairables)
    @param T:           array to be initialised
    @param T_next:      -
*/
void initialise_grid( double dx, Var *v, double** T, double** T_next ){
    int N = v->N;
    *T = malloc( sizeof(double) * N );
    *T_next = malloc( sizeof(double) * N );

    int i;
    for( i=0; i<N; i++){
        double param = dx*i - v->x_T;
        param /= v->x_w;

        (*T)[i] = v->T_0 + 0.5* ( v->T_1 - v->T_0 ) * (1 + tanh(param));
    }
}

/*
    Swap pointers for efficient value exchange

    @param array1:  pointer to be swapped
    @param array2:  -
*/
void swap_mem(double **array1, double **array2) {
    double *tmp;
    tmp     = *array1;
    *array1 = *array2;
    *array2 = tmp;
}

void print_grid( FILE* output, double *T, int len, double time, double dx ){
    int i;
    for(i=0; i<len; i++){
        fprintf(output, "%g %g %g\n", time, i*dx, T[i]);
    }
    printf("\n");
}

double bdd_calculate_a( double T_0, double T_1, double dx, double dT){
    double result = T_0 - T_1 - dT*dx;
    result /= dx*dx;

    return result;
}

/*
    Sets boundary conditions
    @param T:       Array of "old" grid points
    @param T_next:  Array of new grid
    @param v:       Structure containing system's constants
    @param dx:      Grid spacing value
    @param dt:      Time step
*/
void set_bdd_conditions( double *T, double *T_next, Var *v, double dx, double dt ){

    double a;

    // Left side boundary condition (solve a*x^2+c=0)
    //a = bdd_calculate_left_a();
    //c = bdd_calculate_left_c();

     a = bdd_calculate_a( T[ 1 ], T[ 0 ], dx, 0.0 );
     T_next[ 0 ] = 2*thermal_conductivity( T[0], v )*a*dt+T[0];
     //(*T_next)[0] = (*T_next)[1];
    //(*T_next)[0]            = (*T_next)[1];
    //(*T_next)[ v->N - 1 ]   = (*T_next)[ v->N - 2 ];

    // Right side boundary condition (solve for a in a*x^2 + b*x + c = 0 )
    // Then dT/dt = chi( T(t-1, N-1) ) * 2 * a
    double dT = v->P / thermal_conductivity( T[ v->N-1 ], v );
    a = bdd_calculate_a( T[ v->N-2 ], T[ v->N-1 ], -dx, dT );
    T_next[ v->N - 1 ] = 2*thermal_conductivity( T[v->N-1], v )*a*dt+T[v->N-1];

    //double x = dx * (v->N);
    //(*T_next)[v->N-1] = a*(dx*dx) + b*dx + c;

}

/*
    Calculates next time step temperature for all points but boundaries
    @param T:   Temperature grid
    @param v:   Structure containing system's constants
    @param i:   Grid point index
*/
double difference_T( double* T, Var *v, int i){

    double X1 = X_forward(T[i+1], T[i], v);
    double X2 = X_forward(T[i], T[i-1], v);

    double T1 = T[i+1]-T[i];
    double T2 = T[i]-T[i-1];

    return X1*T1-X2*T2;

}

/*
    Returns maximum of the two given parameters
*/
double max( double a, double b ){
    if( a>b )
        return a;
    else
        return b;
}

void open_file( FILE** output ){
    (*output) = fopen("output.txt", "w+");
    if( (*output) == NULL ){
        printf("fopen failed to create/open file\n");
        exit(1);
    }
}
int main()
{
    int i;
    int N;
    int print_loop = 0;
    double dx;
    double d2x;
    double ts;
    double temp_ts;
    double cur_time = 0.0;
    double next_output_time = 0.0;

    double* T = NULL;
    double* T_next = NULL;

    FILE *output = NULL;

    Var v;

    /* INITIALISATION */
    load_struct( &v );
    N = v.N;
    dx = ( v.x_R - v.x_L )/(N-1);
    d2x = dx*dx;
    next_output_time = v.t_d;
    initialise_grid( dx, &v, &T, &T_next );

    // Stability condition
    ts = 0.5*d2x/max(abs(v.X_p), (v.X_s) );

    // Printing set up
    open_file( &output );
    print_grid( output, T, v.N, 0, dx );
    cur_time += ts;


    // MAIN LOOP
    while( cur_time < v.t_f ){
        temp_ts = ts;
        //Makes sure you calculate and print the correct T
        if( cur_time > next_output_time ){
            temp_ts = ts - cur_time + next_output_time;
            cur_time = next_output_time;
            print_loop = 1;
        }
         //loop over points
        for (i=1; i<(N-1); i++) {
            // Centred finite difference evaluation of gradient
            T_next[i] = difference_T( T, &v, i )*temp_ts/d2x+T[i];
        }
        set_bdd_conditions( T, T_next, &v, dx, temp_ts );

        if( print_loop == 1 ){
            //printf("Time: %lf\n", cur_time);
            print_grid( output, T_next, N, cur_time, dx );
            next_output_time += v.t_d;
            print_loop = 0;
        }
        // Efficiently copy next values at time step to y array.
        swap_mem(&T, &T_next);

        cur_time += ts;
    }

    free( T );
    free( T_next );
    fclose( output );
}

