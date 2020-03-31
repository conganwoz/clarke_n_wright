#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

using namespace std;

struct Route
{
  int start;     // start_customer
  int end;       // end_customer
  double length; // sumary length distances
  double load;
  bool is_begin_with_station;
  bool is_merged;
};

// -- data problem --
int dimention = 0;
int num_stations = 0;
int num_customers = 0;
int num_saving_distances = 0;
double max_capacity_vh = 0.0;
double max_energy_vh = 0.0;
double eng_consumtion = 0.0;
double optimal_val = 0.0;

double **distances;
double saving_distances[1000][3];
double *demands;
double **coords;
bool *visited;

// -- data in process --
struct Route *routes = NULL;
int *route_ids;
int *next_arr, *pred_arr;
int **best_station;
double **best_stat_distance;
bool *is_throught_station;
bool *is_interior;
double *available_energy;

// sumary data
double *F_CAP;
double *B_CAP;
double *F_ENERGY;
double *B_ENERGY;

void read_file(char *file_src)
{
  int i;
  FILE *infile;
  infile = fopen(file_src, "r");

  if (infile == NULL)
  {
    printf("READ FILE ERROR\n");
    exit(-1);
  }
  else
  {
    printf("READ FILE SUCCESS\n");
  }

  fscanf(infile, "%lf\n", &optimal_val);
  fscanf(infile, "%d\n", &dimention);
  fscanf(infile, "%d\n", &num_stations);
  fscanf(infile, "%lf\n", &max_capacity_vh);
  fscanf(infile, "%lf\n", &max_energy_vh);
  fscanf(infile, "%lf\n", &eng_consumtion);

  num_customers = dimention + 1;
  dimention = num_customers + num_stations;

  printf("data: dimention: %d - num_customer: %d - capacity_vh: %lf - energy_vh: %lf - energy_consumtion: %lf\n", dimention, num_customers, max_capacity_vh, max_energy_vh, eng_consumtion);
  create_spaces_mem();
  double index, x, y;
  double demand = 0.0;

  for (i = 0; i < num_customers; i++)
  {
    fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
    coords[i][0] = x;
    coords[i][1] = y;
  }

  for (i = num_customers; i < num_stations; i++)
  {
    fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
    coords[i][0] = x;
    coords[i][1] = y;
  }

  for (i = 0; i < num_customers; i++)
  {
    fscanf(infile, "%lf %lf\n", &index, &demand);
    demands[i] = demand;
  }

  fclose(infile);
  printf("\n");
  for (i = 0; i < dimention; i++)
  {
    printf("%d: %lf - %lf\n", i, coords[i][0], coords[i][1]);
  }
  printf("\n");
  for (i = 0; i < num_customers; i++)
  {
    printf("%d:  %lf", i, demands[i]);
  }
}

void create_spaces_mem()
{
  int i = 0;
  coords = (double **)malloc(dimention * sizeof(double*));
  demands = (double *)malloc(num_customers * sizeof(double));
  for (i = 0; i < dimention; i++)
  {
    coords[i] = (double *)malloc(2 * sizeof(double));
  }
}

int main()
{
  read_file("E-n22-k4.evrp");
  return -1;
}
