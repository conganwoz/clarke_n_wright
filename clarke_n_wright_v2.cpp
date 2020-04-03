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

double **distances;               // done
double saving_distances[1000][3]; // done
// double *demands;                  // done
double demands[2000];
double **coords; // done
bool *visited;   // done

// -- data in process --
struct Route *routes = NULL;
int *route_ids;               // done
int *next_arr, *pred_arr;     // done
int **best_stations;          // done
double **best_stat_distances; // done
bool *is_throught_station;    // done
bool *is_interior;            // done
double *available_energy;     // done

// sumary data
double *F_CAP;
double *B_CAP;
double *F_ENERGY;
double *B_ENERGY;

void create_spaces_mem()
{
  int i = 0;
  coords = (double **)malloc(dimention * sizeof(double *));
  distances = (double **)malloc(dimention * sizeof(double *));
  best_stations = (int **)malloc(num_customers * sizeof(int *));
  best_stat_distances = (double **)malloc(num_customers * sizeof(int *));

  // demands = (double *)malloc(num_customers * sizeof(double));
  next_arr = (int *)malloc(num_customers * sizeof(int));
  pred_arr = (int *)malloc(num_customers * sizeof(int));
  is_throught_station = (bool *)malloc(num_customers * sizeof(bool));
  is_interior = (bool *)malloc(num_customers * sizeof(bool));
  visited = (bool *)malloc(num_customers * sizeof(bool));
  route_ids = (int *)malloc(num_customers * sizeof(int));
  available_energy = (double *)malloc(num_customers * sizeof(double));

  routes = (struct Route *)malloc(num_customers * sizeof(struct Route));

  F_CAP = (double *)malloc(num_customers * sizeof(double));
  B_CAP = (double *)malloc(num_customers * sizeof(double));
  F_ENERGY = (double *)malloc(num_customers * sizeof(double));
  B_ENERGY = (double *)malloc(num_customers * sizeof(double));

  for (i = 0; i < dimention; i++)
  {
    coords[i] = (double *)malloc(2 * sizeof(double));
    distances[i] = (double *)malloc(dimention * sizeof(double));
    if (i < num_customers)
    {
      best_stations[i] = (int *)malloc(num_customers * sizeof(int));
      best_stat_distances[i] = (double *)malloc(num_customers * sizeof(double));
    }
  }
}

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

  num_customers = dimention;
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

  for (i = num_customers; i < dimention; i++)
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
}

void find_best_stations(int i, int j)
{
  double min_distance = 100000000;
  int k;
  int best_stat = -1;
  for (k = num_customers; k < dimention; k++)
  {
    if (distances[i][k] + distances[j][k] < min_distance)
    {
      min_distance = distances[i][k] + distances[j][k];
      best_stat = k;
    }
  }

  best_stat_distances[i][j] = min_distance;
  best_stat_distances[j][i] = min_distance;

  best_stations[i][j] = best_stat;
  best_stations[j][i] = best_stat;
}

void bubble_sort(int k)
{
  int i, j;
  for (i = 0; i < k; i++)
  {
    for (j = i + 1; j < k; j++)
    {
      if (saving_distances[j][2] > saving_distances[i][2])
      {
        double tempI = saving_distances[i][0];
        double tempJ = saving_distances[i][1];
        double tempS = saving_distances[i][2];

        saving_distances[i][0] = saving_distances[j][0];
        saving_distances[i][1] = saving_distances[j][1];
        saving_distances[i][2] = saving_distances[j][2];

        saving_distances[j][0] = tempI;
        saving_distances[j][1] = tempJ;
        saving_distances[j][2] = tempS;
      }
    }
  }
}

void prepare_data()
{
  int i, j, k = 0;
  for (i = 0; i < dimention; i++)
  {
    for (j = i + 1; j < dimention; j++)
    {
      double distance = sqrt((coords[i][0] - coords[j][0]) * (coords[i][0] - coords[j][0]) + (coords[i][1] - coords[j][1]) * (coords[i][1] - coords[j][1]));
      distances[i][j] = distance;
      distances[j][i] = distance;
    }
  }

  for (i = 0; i < num_customers; i++)
  {
    is_throught_station[i] = false;
    is_interior[i] = false;
    visited[i] = false;
  }

  // find best station of couple node(i,j)
  for (i = 0; i < num_customers; i++)
  {
    for (j = i + 1; j < num_customers; j++)
      find_best_stations(i, j);
  }

  for (i = 1; i < num_customers; i++)
  {
    for (j = i + 1; j < num_customers; j++)
    {
      double saving = distances[0][i] + distances[0][j] - distances[i][j];
      saving_distances[k][0] = i;
      saving_distances[k][1] = j;
      saving_distances[k][2] = saving;
      k++;
    }
  }
  num_saving_distances = k;
  bubble_sort(num_saving_distances);
}

double dist_consum(double distance)
{
  return distance * eng_consumtion;
}

void init()
{
  int i = 0;
  for (i = 1; i < num_customers; i++)
  {
    int x = demands[1];
    if (demands[i] <= max_capacity_vh)
    {
      int best_stat = best_stations[i][0];
      double best_stat_dis = best_stat_distances[i][0];

      if (dist_consum(distances[i][0] * 2) < max_energy_vh)
      {
        visited[i] = true;
        next_arr[i] = 0;
        pred_arr[i] = 0;
        is_throught_station[i] = false;
        is_interior[i] = false;
        available_energy[i] = max_energy_vh - distances[i][0] * eng_consumtion;
        route_ids[i] = i;
        routes[i].start = i;
        routes[i].end = i;
        routes[i].length = 2 * distances[0][i];
        routes[i].load = demands[i];
        routes[i].is_begin_with_station = false;
        routes[i].is_merged = false;
        F_ENERGY[i] = available_energy[i];
        B_ENERGY[i] = available_energy[i];
      }
      else if (
          (dist_consum(distances[0][best_stat]) < max_energy_vh) && (max_energy_vh > dist_consum(distances[best_stat][i] + distances[i][0])))
      {
        visited[i] = true;
        next_arr[i] = 0;
        pred_arr[i] = 0;
        is_throught_station[i] = false;
        is_interior[i] = false;
        available_energy[i] = max_energy_vh - distances[best_stat][i] * eng_consumtion;
        route_ids[i] = i;
        routes[i].start = i;
        routes[i].end = i;
        routes[i].length = best_stat_distances[0][i] + distances[0][i];
        routes[i].load = demands[i];
        routes[i].is_begin_with_station = true;
        routes[i].is_merged = false;
        F_ENERGY[i] = available_energy[i];
        B_ENERGY[i] = max_energy_vh - dist_consum(distances[i][0]);
      }
      else if (
          dist_consum(distances[0][i] + distances[i][best_stat]) < max_energy_vh && max_energy_vh >= dist_consum(distances[best_stat][0]))
      {
        visited[i] = true;
        next_arr[i] = 0;
        pred_arr[i] = 0;
        is_throught_station[i] = true;
        is_interior[i] = false;
        available_energy[i] = max_energy_vh - dist_consum(distances[i][0]);
        route_ids[i] = i;
        routes[i].start = i;
        routes[i].end = i;
        routes[i].length = best_stat_distances[0][i] + distances[0][i];
        routes[i].load = demands[i];
        routes[i].is_begin_with_station = false;
        routes[i].is_merged = false;
        F_ENERGY[i] = available_energy[i];
        // B_ENERGY[i] =
      }
      else if (
          // (dist_consum(distances[0][i]) > max_energy_vh) && (dist_consum(distances[0][best_stat]) < max_energy_vh) && (max_energy_vh > dist_consum(distances[i][best_stat])) && (max_energy_vh - dist_consum(distances[i][best_stat]) < dist_consum(distances[0][i])) && (max_energy_vh - dist_consum(distances[i][best_stat]) > dist_consum(distances[i][best_stat])) && (max_energy_vh > dist_consum(distances[best_stat][0]))
          // dist_consum(distances[0][i] + distances[i][best_stat] < max_energy_vh) && max_energy_vh >= dist_consum(distances[best_stat][0]))
          (max_energy_vh > dist_consum(distances[0][best_stat])) && (max_energy_vh > dist_consum(distances[best_stat][i] * 2)))
      {
        visited[i] = true;
        next_arr[i] = 0;
        pred_arr[i] = 0;
        is_throught_station[i] = true;
        is_interior[i] = false;
        available_energy[i] = max_energy_vh - dist_consum(distances[best_stat][i]);
        route_ids[i] = i;
        routes[i].start = i;
        routes[i].end = i;
        routes[i].length = 2 * best_stat_distances[0][i];
        routes[i].load = demands[i];
        routes[i].is_begin_with_station = true;
        routes[i].is_merged = false;
      }
      else
      {
        printf("can NOT find route for all Node!");
        exit(-1);
      }
    }
    else
    {
      printf("capacity required of custom exceed max capacity of VH!");
      exit(-1);
    }
  }
}

bool check_is_same_route(int i, int j)
{
  return route_ids[i] == route_ids[j];
}

bool check_is_valid_merge_route_capacity(int i, int j)
{
  int route_i = route_ids[i];
  int route_j = route_ids[j];
  return (routes[route_i].load + routes[route_j].load) < max_capacity_vh;
}

// temp data for validate
double temp_avail_eng[1000];
int temp_is_throught_station[1000];

void reset_temp_data()
{
  int i = 0;
  for (i = 0; i < 1000; i++)
  {
    temp_avail_eng[i] = -1.0;
    temp_is_throught_station[i] = 0;
  }
}

bool check_valid_merge_route_energy(int i, double available_eng, bool is_forward)
{
  reset_temp_data();
  temp_avail_eng[i] = available_eng;
  if (is_forward)
  {
    int x = i;
    while (x != 0)
    {
      int next_customer = next_arr[x];
      int best_stat = best_stations[x][next_customer];
      double best_dis_stat = best_stat_distances[x][next_customer];
      if (dist_consum(distances[x][next_customer]) < temp_avail_eng[x])
      {
        temp_avail_eng[next_customer] = temp_avail_eng[x] - dist_consum(distances[x][next_customer]);
        if (temp_avail_eng[next_customer] <= 0)
          return false;
        temp_is_throught_station[x] = -1;
      }
      else if (temp_avail_eng[x] > dist_consum(distances[x][best_stat]) && (max_energy_vh > dist_consum(distances[best_stat][next_customer])))
      {
        temp_avail_eng[next_customer] = max_energy_vh - dist_consum(distances[best_stat][next_customer]);
        if (temp_avail_eng[next_customer] <= 0)
          return false;
        temp_is_throught_station[x] = 1;
      }
      else
      {
        return false;
      }
      x = next_arr[x];
    }
  }
  else
  {
    int x = i;
    while (x != 0)
    {
      int pred_customer = pred_arr[x];
      int best_stat = best_stations[x][pred_customer];
      double best_dis_stat = best_stat_distances[x][best_stat];
      if (dist_consum(distances[x][pred_customer]) < temp_avail_eng[x])
      {
        temp_avail_eng[pred_customer] = temp_avail_eng[x] - dist_consum(distances[x][pred_customer]);
        if (temp_avail_eng[pred_customer] <= 0)
          return false;
        temp_is_throught_station[x] = -1;
      }
      else if (temp_avail_eng[x] > dist_consum(distances[x][best_stat]) && (max_energy_vh > dist_consum(distances[best_stat][pred_customer])))
      {
        temp_avail_eng[pred_customer] = max_energy_vh - dist_consum(distances[best_stat][pred_customer]);
        if (temp_avail_eng[pred_customer] <= 0)
          return false;
        temp_is_throught_station[x] = 1;
      }
      else
      {
        return false;
      }
      x = pred_arr[x];
    }
  }

  for (int m = 0; m < 1000; m++)
  {
    if (temp_avail_eng[m] > 0)
      available_energy[m] = temp_avail_eng[m];
    if (temp_is_throught_station[m] == 1)
    {
      is_throught_station[m] = true;
    }
    else if (temp_is_throught_station[m] == -1)
    {
      is_throught_station[m] = false;
    }
  }

  return true;
}

double revert_route(int i, bool is_forward)
{
  double interior_length = 0;
  if (is_forward)
  {
    int x = i;
    while (x != 0)
    {
      if (x == i)
      {
        pred_arr[x] = next_arr[x];
        if (next_arr[x] != 0)
          interior_length += distances[x][next_arr[x]];
        x = next_arr[x];
      }
      else
      {
        int temp = pred_arr[x];
        pred_arr[x] = next_arr[x];
        next_arr[x] = temp;
        if (pred_arr[x] != 0)
          interior_length += distances[x][pred_arr[x]];
        x = pred_arr[x];
      }
      //interior_length += distances[x][]
      // x = next_arr[x];
    }
  }
  else
  {
    int x = i;
    while (x != 0)
    {
      if (x == i)
      {
        next_arr[x] = pred_arr[x];
        if (pred_arr[x] != 0)
          interior_length += distances[x][pred_arr[x]];
        x = pred_arr[x];
      }
      else
      {
        int temp = pred_arr[x];
        pred_arr[x] = next_arr[x];
        next_arr[x] = temp;
        if (next_arr[x] != 0)
          interior_length += distances[x][next_arr[x]];
        x = next_arr[x];
      }
    }
  }

  int start = routes[route_ids[i]].start;
  int end = routes[route_ids[i]].end;

  routes[route_ids[i]].start = end;
  routes[route_ids[i]].end = start;
  return interior_length;
}

void set_interior(int i)
{
  if (next_arr[i] == 0 && pred_arr[i] == 0)
    is_interior[i] = false;
  else
  {
    is_interior[i] = true;
  }
}

double available_after_j(int j)
{
  int start = routes[route_ids[j]].end;
  int x = start;
  double last_avaible = 0.0;
  if (max_energy_vh > dist_consum(distances[x][0]))
  {
    last_avaible = max_energy_vh - dist_consum(distances[x][0]);
  }
  else if (max_energy_vh > dist_consum(distances[0][best_stations[0][x]]) && max_energy_vh > dist_consum(distances[x][best_stations[0][x]]))
  {
    last_avaible = max_energy_vh - dist_consum(distances[x][best_stations[0][x]]);
  }
  else
    return -1.0;
  x = pred_arr[x];

  while (x != 0)
  {
    int pred_customer = next_arr[x];
    if (last_avaible > dist_consum(distances[pred_customer][x]))
    {
      last_avaible -= dist_consum(distances[pred_customer][x]);
    }
    else if (last_avaible > dist_consum(distances[best_stations[pred_customer][x]][pred_customer]) && max_capacity_vh > dist_consum(distances[best_stations[pred_customer][x]][x]))
    {
      last_avaible = max_energy_vh - dist_consum(distances[best_stations[x][pred_customer]][x]);
    }
    else
      return -1.0;

    x = pred_arr[x];
  }
  return last_avaible;
}

double get_interior_length(int i)
{
  int route_i = route_ids[i];
  int start = routes[route_i].start;
  int end = routes[route_i].end;

  double begin_length = 0.0;
  double end_length = 0.0;

  if (routes[i].is_begin_with_station)
    begin_length = best_stat_distances[start][0];
  else
    begin_length = distances[start][0];

  if (is_throught_station[end])
    end_length = best_stat_distances[end][0];
  else
    end_length = distances[end][0];

  return routes[route_i].length - begin_length - end_length;
}

void set_route_id(int i, int new_route_id)
{
  int start = routes[route_ids[i]].start;
  int end = routes[route_ids[i]].end;
  int x = start;

  while (x != 0)
  {
    if (x != i)
    {
      route_ids[x] = new_route_id;
    }
    x = next_arr[x];
  }
}

void temp_process(int i, int j)
{
  int best_stat = best_stations[i][j];
  double interior_length_j = get_interior_length(j);
  if (available_energy[i] > dist_consum(distances[i][j]))
  {
    if (check_valid_merge_route_energy(j, available_energy[i] - dist_consum(distances[i][j]), true))
    {
      set_route_id(j, route_ids[i]);
      int route_i = route_ids[i];
      int route_j = route_ids[j];
      set_interior(i);
      set_interior(j);
      next_arr[i] = j;
      pred_arr[j] = i;

      is_throught_station[i] = false;

      routes[route_ids[i]].end = routes[route_ids[j]].end;
      routes[route_ids[i]].length += (interior_length_j + distances[i][j]);
      routes[route_ids[i]].load += routes[route_ids[j]].load;
      routes[route_ids[j]].is_merged = true;
      route_ids[j] = route_ids[i];
    }
    else if ((available_energy[i] > dist_consum(distances[i][best_stat])) && (max_energy_vh > dist_consum(distances[best_stat][j])))
    {
      set_route_id(j, route_ids[i]);
      pred_arr[j] = i;
      next_arr[i] = j;
      is_throught_station[i] = true;

      set_interior(i);
      set_interior(j);

      routes[route_ids[i]].end = routes[route_ids[j]].end;
      routes[route_ids[i]].length += (interior_length_j + best_stat_distances[i][j]);
      routes[route_ids[i]].load += routes[route_ids[j]].load;
      routes[route_ids[j]].is_merged = true;
      route_ids[j] = route_ids[i];
    }
  }
}

void merge_route(int i, int j)
{
  int x = -1;
  if (!visited[i] && !visited[j])
    x = 0;
  if (visited[i] && !visited[j])
    x = 1;
  if (!visited[i] && visited[j])
    x = 2;
  if (visited[i] && visited[j])
    x = 3;
  if (is_interior[i] || is_interior[j])
    x = -1;

  if (x == -1)
    return;

  switch (x)
  {
  case 0:
    break; // 2 node i va j khong the toi tu 0 --> noi vs nhau cung khong the toi 0
  case 1:
    break; // tuong tu nhu tren
  case 2:
    break; //
  case 3:
    if (check_is_valid_merge_route_capacity(i, j) && !is_interior[i] && !is_interior[j])
    {
      int best_stat = best_stations[i][j];
      double best_distance_stat = best_stat_distances[i][j];

      if (next_arr[i] == 0 && next_arr[j] == 0)
      {
        if (available_energy[i] > dist_consum(distances[i][j]))
        {
          if (check_valid_merge_route_energy(j, available_energy[i] - dist_consum(distances[i][j]), false))
          {
            set_route_id(j, route_ids[i]);
            double interior_length = revert_route(j, false);
            set_interior(i);
            set_interior(j);
            pred_arr[j] = i;
            next_arr[i] = j;
            is_throught_station[i] = false;

            // setup route info
            routes[route_ids[i]].end = routes[route_ids[j]].end;
            routes[route_ids[i]].length += (interior_length + distances[i][j]);
            routes[route_ids[i]].load += routes[route_ids[j]].load;
            routes[route_ids[j]].is_merged = true;
            route_ids[j] = route_ids[i];
          }
        }
        else if ((available_energy[i] > dist_consum(distances[i][best_stat])) && (max_energy_vh > dist_consum(distances[best_stat][j])))
        {
          set_route_id(j, route_ids[i]);
          double interior_length = revert_route(j, false);
          set_interior(i);
          set_interior(j);
          pred_arr[j] = i;
          next_arr[i] = j;
          is_throught_station[i] = true;

          // setup route info
          routes[route_ids[i]].end = routes[route_ids[j]].end;
          routes[route_ids[i]].length += (interior_length + best_stat_distances[i][j]);
          routes[route_ids[i]].load += routes[route_ids[j]].load;
          routes[route_ids[j]].is_merged = true;
          route_ids[j] = route_ids[i];
        }
      }
      else if (pred_arr[i] == 0 && pred_arr[j] == 0)
      {
        double temp_avail_eng = available_after_j(j);
        if (temp_avail_eng > dist_consum(distances[i][j]))
        {
          if (check_valid_merge_route_energy(i, temp_avail_eng - dist_consum(distances[i][j]), true))
          {
            set_route_id(i, route_ids[j]);
            double interior_length = revert_route(j, true);
            set_interior(i);
            set_interior(j);
            pred_arr[i] = j;
            next_arr[j] = i;

            is_throught_station[j] = false;

            // setup route info
            routes[route_ids[j]].end = routes[route_ids[i]].end;
            routes[route_ids[j]].length += (interior_length + distances[i][j]);
            routes[route_ids[j]].load += routes[route_ids[i]].load;
            routes[route_ids[i]].is_merged = true;
            route_ids[i] = route_ids[j];
          }
        }
        else if (temp_avail_eng > dist_consum(distances[best_stat][i]) && max_energy_vh > dist_consum(distances[best_stat][j]))
        {
          if (check_valid_merge_route_energy(i, max_energy_vh - dist_consum(distances[best_stat][j]), true))
          {
            set_route_id(i, route_ids[j]);
            double interior_length = revert_route(j, true);
            set_interior(i);
            set_interior(j);
            pred_arr[i] = j;
            next_arr[j] = i;

            is_throught_station[j] = true;

            routes[route_ids[j]].end = routes[route_ids[i]].end;
            routes[route_ids[j]].length += (interior_length + best_stat_distances[i][j]);
            routes[route_ids[j]].load += routes[route_ids[i]].load;
            routes[route_ids[i]].is_merged = true;
            route_ids[i] = route_ids[j];
          }
        }
      }
      else if (next_arr[i] == 0 && pred_arr[j] == 0)
      {
        temp_process(i, j);
      }
      else if (pred_arr[i] == 0 && next_arr[j] == 0)
      {
        temp_process(j, i);
      }
    }
  default:
    break;
  }
}

void show_result()
{
  int i;
  for (i = 1; i < dimention; i++)
  {
    if (!routes[i].is_merged)
    {
      int start = routes[i].start;
      int end = routes[i].end;

      int x = start;
      while (x != 0)
      {
        printf("%d ", x);
        x = next_arr[x];
      }
      printf("\n");
    }
  }
}

void clarke_wright()
{
  int m;
  for (m = 0; m < num_saving_distances; m++)
  {
    int i, j, saving;
    i = (int)saving_distances[m][0];
    j = (int)saving_distances[m][1];
    saving = (int)saving_distances[m][2];
    if (!check_is_same_route(i, j))
    {
      merge_route(i, j);
      // printf("\n\n\n=================: %d\n", m);
      // show_result();
      // printf("\n=================\n\n\n");
    }
  }
}

int main()
{
  read_file("E-n22-k4.evrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
  printf("===================================\n\n");

  read_file("./data_set/E-n23-k3.evrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
  printf("===================================\n\n");

  read_file("./data_set/E-n30-k3.evrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
  printf("===================================\n\n");

  read_file("./data_set/E-n33-k4.evrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
  printf("===================================\n\n");

  read_file("./data_set/E-n51-k5.evrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
  printf("===================================\n\n");

  // read_file("./data_set/E-n76-k7.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n101-k8.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n143-k7.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n214-k11.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n351-k40.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n459-k26.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n573-k30.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n685-k75.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n749-k98.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n819-k171.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n916-k207.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  // read_file("./data_set/E-n1001-k43.evrp");
  // prepare_data();
  // init();
  // clarke_wright();
  // show_result();
  // printf("===================================\n\n");

  return -1;
}
