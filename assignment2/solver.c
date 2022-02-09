/* Code for the Jacbi equation solver. 
 * Author: Naga Kandasamy
 * Date modified: January 25, 2022
 *
 * Compile as follows:
 * gcc -o solver solver.c solver_gold.c -O3 -Wall -std=c99 -lm -lpthread
 *
 * If you wish to see debug info add the -D DEBUG option when compiling the code.
 */

// System includes
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#include "grid.h" 

extern int compute_gold(grid_t *);
int compute_using_pthreads_jacobi(grid_t *, int);
void compute_grid_differences(grid_t *, grid_t *);
grid_t *create_grid(int, float, float);
grid_t *copy_grid(grid_t *);
void print_grid(grid_t *);
void print_stats(grid_t *);
double grid_mse(grid_t *, grid_t *);


int main(int argc, char **argv)
{	
	if (argc < 5) {
        fprintf(stderr, "Usage: %s grid-dimension num-threads min-temp max-temp\n", argv[0]);
        fprintf(stderr, "grid-dimension: The dimension of the grid\n");
        fprintf(stderr, "num-threads: Number of threads\n"); 
        fprintf(stderr, "min-temp, max-temp: Heat applied to the north side of the plate is uniformly distributed between min-temp and max-temp\n");
        exit(EXIT_FAILURE);
    }
    
    /* Parse command-line arguments. */
    int dim = atoi(argv[1]);
    int num_threads = atoi(argv[2]);
    float min_temp = atof(argv[3]);
    float max_temp = atof(argv[4]);
    
    /* Generate grids and populate them with initial conditions. */
 	grid_t *grid_1 = create_grid(dim, min_temp, max_temp);
    /* Grid 2 should have the same initial conditions as Grid 1. */
    grid_t *grid_2 = copy_grid(grid_1); 

	/* Compute reference solution using the single-threaded version. */
	fprintf(stderr, "\nUsing the single threaded version to solve the grid\n");
	int num_iter = compute_gold(grid_1);
	fprintf(stderr, "Convergence achieved after %d iterations\n", num_iter);
    /* Print key statistics for the converged values. */
	fprintf(stderr, "Printing statistics for the interior grid points\n");
    print_stats(grid_1);
#ifdef DEBUG
    print_grid(grid_1);
#endif
	
	/* Use pthreads to solve the equation using the jacobi method. */
	fprintf(stderr, "\nUsing pthreads to solve the grid using the jacobi method\n");
	num_iter = compute_using_pthreads_jacobi(grid_2, num_threads);
	fprintf(stderr, "Convergence achieved after %d iterations\n", num_iter);			
    fprintf(stderr, "Printing statistics for the interior grid points\n");
	print_stats (grid_2);
#ifdef DEBUG
    print_grid (grid_2);
#endif
    
    /* Compute grid differences. */
    double mse = grid_mse(grid_1, grid_2);
    fprintf (stderr, "MSE between the two grids: %f\n", mse);

	/* Free up the grid data structures. */
	free((void *) grid_1->element);	
	free((void *) grid_1); 
	free((void *) grid_2->element);	
	free((void *) grid_2);

	exit(EXIT_SUCCESS);
}

/* FIXME: Edit this function to use the jacobi method of solving the equation. The final result should be placed in the grid data structure. */

typedef struct pthread_args {
	int pthread_id;
	int pthread_num;
	grid_t *grid1;
	grid_t *grid2;

	double *pthread_diff;
} PTHREAD_ARGS;

int pthread_iter = 0;
int pthread_done = 0;

pthread_barrier_t barrier1;
pthread_barrier_t barrier2;
pthread_barrier_t barrier3;

int compute_using_pthreads_jacobi(grid_t *grid, int num_threads)
{
	int i;
	
	PTHREAD_ARGS *pthread_args;
	double pthread_diffs[num_threads];

	grid_t *grid2 = copy_grid(grid);

	pthread_t *worker = (pthread_t *) malloc (num_threads *sizeof (pthread_t));

	pthread_barrier_init(&barrier1, 0, num_threads);
	pthread_barrier_init(&barrier2, 0, num_threads);
	pthread_barrier_init(&barrier3, 0, num_threads);

	for(i = 0; i < num_threads; i++){
		pthread_args = (PTHREAD_ARGS *) malloc (sizeof (ARGS_FOR_THREAD));

		pthread_args->pthread_id = i;
		pthread_args->pthread_num = num_threads;
		pthread_args->grid1 = grid;
		pthread_args->grid2 = grid2;
		pthread_args->pthread_diff = pthread_diffs + i;

		if ((pthread_create (&worker[i], NULL, jacobi, (void *) pthread_args)) != 0){
			exit(EXIT_FAILURE);
		}
	}

	for(i = 0; i < num_threads; i++){
		pthread_join(worker[i], NULL);
	}

	free((void *) grid2->element);
	free((void *) grid2);
	pthread_barrier_destroy(&barrier1);
	pthread_barrier_destroy(&barrier2);
	pthread_barrier_destroy(&barrier3);

    return pthread_iters;
}

void *jacobi(void *pthread_args)
{
	
}

/* Create a grid with the specified initial conditions. */
grid_t* create_grid(int dim, float min, float max)
{
    grid_t *grid = (grid_t *)malloc (sizeof(grid_t));
    if (grid == NULL)
        return NULL;

    grid->dim = dim;
	fprintf(stderr, "Creating a grid of dimension %d x %d\n", grid->dim, grid->dim);
	grid->element = (float *) malloc(sizeof(float) * grid->dim * grid->dim);
    if (grid->element == NULL)
        return NULL;

    int i, j;
	for (i = 0; i < grid->dim; i++) {
		for (j = 0; j < grid->dim; j++) {
            grid->element[i * grid->dim + j] = 0.0; 			
		}
    }

    /* Initialize the north side, that is row 0, with temperature values. */ 
    srand((unsigned)time(NULL));
	float val;		
    for (j = 1; j < (grid->dim - 1); j++) {
        val =  min + (max - min) * rand ()/(float)RAND_MAX;
        grid->element[j] = val; 	
    }

    return grid;
}

/* Creates a new grid and copies over the contents of an existing grid into it. */
grid_t* copy_grid(grid_t *grid) 
{
    grid_t *new_grid = (grid_t *)malloc(sizeof(grid_t));
    if (new_grid == NULL)
        return NULL;

    new_grid->dim = grid->dim;
	new_grid->element = (float *)malloc(sizeof(float) * new_grid->dim * new_grid->dim);
    if (new_grid->element == NULL)
        return NULL;

    int i, j;
	for (i = 0; i < new_grid->dim; i++) {
		for (j = 0; j < new_grid->dim; j++) {
            new_grid->element[i * new_grid->dim + j] = grid->element[i * new_grid->dim + j] ; 			
		}
    }

    return new_grid;
}

/* Print grid to screen. */
void print_grid(grid_t *grid)
{
    int i, j;
    for (i = 0; i < grid->dim; i++) {
        for (j = 0; j < grid->dim; j++) {
            printf("%f\t", grid->element[i * grid->dim + j]);
        }
        printf("\n");
    }
    printf("\n");
}


/* Print out statistics for the converged values of the interior grid points, including min, max, and average. */
void print_stats(grid_t *grid)
{
    float min = INFINITY;
    float max = 0.0;
    double sum = 0.0;
    int num_elem = 0;
    int i, j;

    for (i = 1; i < (grid->dim - 1); i++) {
        for (j = 1; j < (grid->dim - 1); j++) {
            sum += grid->element[i * grid->dim + j];

            if (grid->element[i * grid->dim + j] > max) 
                max = grid->element[i * grid->dim + j];

             if(grid->element[i * grid->dim + j] < min) 
                min = grid->element[i * grid->dim + j];
             
             num_elem++;
        }
    }
                    
    printf("AVG: %f\n", sum/num_elem);
	printf("MIN: %f\n", min);
	printf("MAX: %f\n", max);
	printf("\n");
}

/* Calculate the mean squared error between elements of two grids. */
double grid_mse (grid_t *grid_1, grid_t *grid_2)
{
    double mse = 0.0;
    int num_elem = grid_1->dim * grid_1->dim;
    int i;

    for (i = 0; i < num_elem; i++) 
        mse += (grid_1->element[i] - grid_2->element[i]) * (grid_1->element[i] - grid_2->element[i]);
                   
    return mse/num_elem; 
}



		

