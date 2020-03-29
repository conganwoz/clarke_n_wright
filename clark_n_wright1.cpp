#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

using namespace std;

struct route
{
  int start;
  int end;
  double length;
  double capacity;
  struct route *next;
  bool is_begin_with_station;
  bool is_merged;
};

const int ADDED = 1;
const int UNUSED = 2;
const int INTERIOR = 3;

// Thong so khoi tao CT
struct route *ROUTE = NULL;
struct route *END_ROUTE = NULL;
int *next_arr, *pred_arr; // tracking position
int **station;            // best station between i and j
int **dis_station;        // cost distance through station
bool *is_throught_station;
bool *is_interior;
double *availble_energy;

// data problem
double distances[1000][1000];
double saving_distances[1000][3];
double demands[1000];
double coords[1000][2];

int dimention = 0;    // so luong node cua bai toan
int num_customer = 0; // so luong customer

double capacity = 0;         // capacity cua xe
double energy = 0;           // max energy of xe
int num_saving_distance = 0; // do dai mang saving_distance
double consum_coeff = 0.2;
double max_length = 1000;

// Thong so trong qua trinh chay CT
bool visited[1000];
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
  }
  else
    printf("READ FILE SUCCESS\n");

  fscanf(infile, "%d\n", &dimention);
  fscanf(infile, "%d\n", &num_customer);
  fscanf(infile, "%lf\n", &capacity);
  fscanf(infile, "%lf\n", &energy);

  printf("data: dimension: %d - num_customer: %d - capacity: %lf - energy: %lf\n", dimention, num_customer, capacity, energy);

  double index, x, y;
  double demand = 0;

  //int num_station = dimention - num_customer;
  for (i = 0; i < num_customer; i++)
  {
    fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
    coords[i][0] = x;
    coords[i][1] = y;
  }

  for (i = num_customer; i < dimention; i++)
  {
    fscanf(infile, "%lf %lf %lf\n", &index, &x, &y);
    coords[i][0] = x;
    coords[i][1] = y;
  }

  for (i = 0; i < num_customer; i++)
  {
    fscanf(infile, "%lf %lf\n", &index, &demand);
    demands[i] = demand;
  }

  // set init data
  for (i = 0; i < num_customer; i++)
  {
    visited[i] = false;
  }

  fclose(infile);
}

int find_best_station(int i, int j)
{
  double min_distance = 10000000;
  int k;
  int best_station = -1;
  for (k = num_customer; k < dimention; k++)
  {
    if ((distances[i][k] + distances[j][k]) < min_distance)
    {
      min_distance = distances[i][k] + distances[j][k];
      best_station = k;
    }
  }

  dis_station[i][j] = min_distance;
  dis_station[j][i] = min_distance;

  return best_station;
}

void bubble_sort(int k)
{
  int i, j;
  for (int i = 0; i < k; i++)
  {
    for (int j = i + 1; j < k; j++)
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
  // create distance matrix
  for (i = 0; i < dimention; i++)
  {
    for (j = 0; j < dimention; j++)
    {
      double temp = sqrt((coords[i][0] - coords[j][0]) * (coords[i][0] - coords[j][0]) + (coords[i][1] - coords[j][1]) * (coords[i][1] - coords[j][1]));
      distances[i][j] = temp;
      distances[j][i] = temp;
    }
  }

  // malloc space for store data
  next_arr = (int *)malloc(num_customer * sizeof(int));
  pred_arr = (int *)malloc(num_customer * sizeof(int));

  is_throught_station = (bool *)malloc(num_customer * sizeof(bool));
  is_interior = (bool *)malloc(num_customer * sizeof(bool));
  availble_energy = (double *)malloc(num_customer * sizeof(double));
  B_CAP = (double *)malloc(num_customer * sizeof(double));
  F_CAP = (double *)malloc(num_customer * sizeof(double));
  F_ENERGY = (double *)malloc(num_customer * sizeof(double));
  B_ENERGY = (double *)malloc(num_customer * sizeof(double));
  for (i = 0; i < num_customer; i++)
  {
    is_throught_station[i] = false;
    is_interior[i] = false;
  }

  station = (int **)malloc(sizeof(int *) * num_customer);
  dis_station = (int **)malloc(sizeof(int *) * num_customer);

  for (i = 0; i < num_customer; i++)
  {
    station[i] = (int *)malloc(sizeof(int) * num_customer);
    dis_station[i] = (int *)malloc(sizeof(int) * num_customer);
  }

  // prepare meta data
  for (i = 0; i < num_customer; i++)
  {
    for (j = i + 1; j < num_customer; j++)
    {
      int best = find_best_station(i, j);
      station[i][j] = best;
      station[j][i] = best;
    }
  }

  for (i = 1; i < num_customer; i++)
  {
    for (j = i + 1; j < num_customer; j++)
    {
      double s = distances[0][i] + distances[0][j] - distances[i][j];
      saving_distances[k][0] = i;
      saving_distances[k][1] = j;
      saving_distances[k][2] = s;
      k++;
    }
  }

  num_saving_distance = k;
  bubble_sort(num_saving_distance);
  double sum = 0;
  for (i = 0; i < num_customer; i++)
  {
    for (j = i + 1; j < num_customer; j++)
    {
      sum += distances[i][j];
    }
  }
  printf("avg: %lf", sum / num_customer);
}

bool check_in_the_same_route(int i, int j)
{
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    bool check_i = false;
    bool check_j = false;
    int start = temp->start;
    int end = temp->end;
    int x = start;
    while (x != end && x != -1)
    {
      if (i == x)
        check_i = true;
      if (j == x)
        check_j = true;
      if (check_i && check_j)
        return true;
      x = next_arr[x];
    }
    if (i == end)
      check_i = true;
    if (j == end)
      check_j = true;
    if (check_i && check_j)
      return true;
    temp = temp->next;
  }

  return false;
}

void init()
{
  ROUTE = (struct route *)malloc(sizeof(struct route));
  END_ROUTE = ROUTE;
  for (int i = 1; i < num_customer; i++)
  {
    if (demands[i] < capacity)
    {
      if (max_length > distances[i][0])
      {
        availble_energy[i] = max_length - distances[0][i];
        next_arr[i] = -1;
        pred_arr[i] = -1;
        is_interior[i] = false;
        visited[i] = true;
        END_ROUTE->next = (struct route *)malloc(sizeof(struct route));
        END_ROUTE = END_ROUTE->next;
        END_ROUTE->start = i;
        END_ROUTE->end = i;
        END_ROUTE->length = 2 * distances[0][i];
        END_ROUTE->capacity = demands[i];
        END_ROUTE->is_merged = false;
        END_ROUTE->next = NULL;
        if (availble_energy[i] < distances[i][0])
          is_throught_station[i] = true;
        continue;
      }
    }
    else
    {
      printf("error_demand");
      exit(-1);
    }
  }
}

bool is_validate_new_node(int i, int j)
{
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    int start = temp->start;
    int end = temp->end;
    int x = start;
    while (x != end)
    {
      if (x == i)
        return temp->capacity + demands[j] < capacity;
      x = next_arr[x];
    }
    if (x == i)
      return temp->capacity + demands[j] < capacity;

    temp = temp->next;
  }
  return false;
}

struct route *update_bound(int i, int j, bool is_tail, double distance, bool is_init_throught_station, int k)
{
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    int start = temp->start;
    int end = temp->end;
    int x = start;
    if (!is_tail && x == i)
    {
      if (j == -1)
      {
        temp->start = i;
        temp->length = temp->length - distance;
        temp->capacity = temp->capacity - demands[k];
      }
      temp->start = j;
      temp->length = temp->length + distance;
      temp->capacity = temp->capacity + demands[j];
      temp->is_begin_with_station = is_init_throught_station;
      return temp;
    }
    else
    {
      while (x != end)
      {
        x = next_arr[x];
      }
      if (x == i)
      {
        temp->end = j;
        temp->length = temp->length + distance;
        temp->capacity = temp->capacity + demands[j];
        temp->is_begin_with_station = is_init_throught_station;
        return temp;
      }
    }
    temp = temp->next;
  }
  return NULL;
}

void add_new_node_to_route(int i, int j)
{
  if (next_arr[i] == -1)
  {
    if (availble_energy[i] > distances[i][j])
    {
      availble_energy[j] = availble_energy[i] - distances[i][j];
      next_arr[i] = j;
      pred_arr[j] = i;
      next_arr[j] = -1;
      is_interior[i] = true;
      is_interior[j] = false;
      is_throught_station[i] = false;
      update_bound(i, j, true, distances[i][j] + distances[j][0] - distances[i][0], false, 0);
    }
    else
    {
      availble_energy[j] = max_length - distances[j][station[i][j]];
      next_arr[i] = j;
      pred_arr[j] = i;
      next_arr[j] = -1;
      is_interior[i] = true;
      is_interior[j] = false;
      is_throught_station[i] = true;
      update_bound(i, j, true, dis_station[i][j] + distances[j][0] - distances[i][0], false, 0);
    }
  }
  else
  {
    if (max_length > distances[0][j])
    {
      availble_energy[j] = max_length - distances[0][j];
      if (availble_energy[j] > distances[i][j])
      {
        availble_energy[i] = availble_energy[j] - distances[i][j];
        pred_arr[j] = -1;
        next_arr[j] = i;
        pred_arr[i] = j;
        is_interior[i] = true;
        is_interior[j] = false;
        update_bound(i, j, false, distances[i][j] + distances[j][0] - distances[i][0], false, 0);
      }
      else
      {
        availble_energy[i] = max_length - distances[station[i][j]][i];
        pred_arr[j] = -1;
        next_arr[j] = i;
        pred_arr[i] = j;
        is_interior[i] = true;
        is_interior[j] = false;
        is_throught_station[j] = true;
        update_bound(i, j, false, dis_station[0][j] + distances[j][i] - distances[i][0], false, 0);
      }
    }
    else
    {
      availble_energy[j] = max_length - distances[station[0][j]][j];
      if (availble_energy[j] > distances[i][j])
      {
        availble_energy[i] = availble_energy[j] - distances[i][j];
        pred_arr[j] = -1;
        next_arr[j] = i;
        pred_arr[i] = j;
        is_interior[i] = true;
        is_interior[j] = false;
        update_bound(i, j, false, dis_station[0][j] + distances[j][i] - distances[i][0], true, 0);
      }
      else
      {
        availble_energy[i] = max_length - distances[station[i][j]][i];
        pred_arr[j] = -1;
        next_arr[j] = i;
        pred_arr[i] = j;
        is_throught_station[j] = true;
        is_interior[i] = true;
        is_interior[j] = false;
        update_bound(i, j, false, dis_station[0][j] + distances[j][i] - distances[i][0], true, 0);
      }
    }
  }
}

bool is_valid_merge_route(int i, int j)
{
  struct route *tempi = ROUTE->next;
  struct route *tempj = ROUTE->next;

  while (tempi != NULL)
  {
    int start = tempi->start;
    int end = tempi->end;
    int x = start;
    bool is_found = false;
    while (x != end)
    {
      if (x == i)
      {
        is_found = true;
        break;
      }
      x = next_arr[x];
    }
    if (x == i)
      is_found = true;
    if (is_found)
      break;

    tempi = tempi->next;
  }

  while (tempj != NULL)
  {
    int start = tempj->start;
    int end = tempj->end;
    int x = start;
    bool is_found = false;
    while (x != end)
    {
      if (x == j)
      {
        is_found = true;
        break;
      }
      x = next_arr[x];
    }
    if (x == j)
      is_found = true;
    if (is_found)
      break;

    tempj = tempj->next;
  }

  if (tempi != NULL && tempj != NULL)
  {
    return tempi->capacity + tempj->capacity < capacity;
  }
  else
    return false;
}

struct route *find_route(int i)
{
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    if (!(temp->is_merged))
    {
      int start = temp->start;
      int end = temp->end;
      int x = start;

      while (x != end)
      {
        if (x == i)
          return temp;
        x = next_arr[x];
      }
      if (x == i)
        return temp;
    }
    temp = temp->next;
  }

  return NULL;
}

void revert_route(int i)
{
  struct route *routei = find_route(i);
  int start = routei->start;
  int end = routei->end;

  int x = end;
  while (x != start && x != -1)
  {
    int true_pred_x = pred_arr[x];
    if (x == end)
    {
      next_arr[x] = pred_arr[x];
      pred_arr[x] = -1;
      availble_energy[x] = max_length - distances[end][0];
    }
    else
    {
      int temp = pred_arr[x];
      pred_arr[x] = next_arr[x];
      next_arr[x] = temp;
      availble_energy[x] = availble_energy[pred_arr[x]] - distances[x][pred_arr[x]];
    }
    x = pred_arr[true_pred_x];
  }
  if (x != -1)
  {
    int temp_v = pred_arr[x];
    pred_arr[x] = next_arr[x];
    next_arr[x] = temp_v;
    availble_energy[x] = availble_energy[pred_arr[x]] - distances[x][pred_arr[x]];

    routei->start = end;
    routei->end = start;
    pred_arr[routei->start] = -1;
    next_arr[routei->end] = -1;
  }
}

int update_route(int i)
{
  struct route *routei = find_route(i);
  int start = routei->start;
  int end = routei->end;
  int x = start;
  while (x != end)
  {
    int pred_x = pred_arr[x];
    if (availble_energy[pred_x] > distances[pred_x][x])
    {
      availble_energy[x] = availble_energy[pred_x] - distances[pred_x][x];
    }
    else
    {
      availble_energy[x] = max_length - distances[station[pred_x][x]][x];
      is_throught_station[pred_x] = true;
    }
    x = next_arr[x];
  }
  int pred_x = pred_arr[x];
  if (pred_x == -1)
  {
    routei->is_merged = true;
    return x;
  }
  else
  {
    if (availble_energy[pred_x] > distances[pred_x][x])
    {
      availble_energy[x] = availble_energy[pred_x] - distances[pred_x][x];
    }
    else
    {
      availble_energy[x] = max_length - distances[station[pred_x][x]][x];
      is_throught_station[pred_x] = true;
    }
  }

  routei->is_merged = true;

  return end;
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
    break;
  case 1:
    if (!is_interior[i] && is_validate_new_node(i, j))
    {
      add_new_node_to_route(i, j);
    }
    break;
  case 2:
    if (!is_interior[j] && is_validate_new_node(i, j))
    {
      add_new_node_to_route(j, i);
    }
    break;
  case 3:
    if (is_valid_merge_route(i, j) && !is_interior[i] && !is_interior[j])
    {
      if (pred_arr[i] == -1 && pred_arr[j] == -1)
      {
        struct route *routei = find_route(i);
        struct route *routej = find_route(j);
        if (routei->start == routei->end)
          is_interior[i] = true;
        if (routej->start == routej->end)
          is_interior[j] = true;

        revert_route(j);
        if (availble_energy[j] > distances[i][j])
        {
          availble_energy[i] = availble_energy[j] - distances[i][j];
        }
        else
        {
          availble_energy[i] = max_length - distances[station[i][j]][i];
          is_throught_station[j] = true;
        }
        int end = update_route(i);
        struct route *r = find_route(j);
        next_arr[j] = i;
        pred_arr[i] = j;
        r->end = end;
      }
      else if (next_arr[i] == -1 && next_arr[j] == -1)
      {
        revert_route(j);
        if (availble_energy[i] > distances[i][j])
        {
          availble_energy[j] = availble_energy[i] - distances[i][j];
        }
        else
        {
          availble_energy[j] = max_length - distances[station[i][j]][j];
          is_throught_station[i] = true;
        }
        int end = update_route(j);
        struct route *r = find_route(i);
        next_arr[i] = j;
        pred_arr[j] = i;
        r->end = end;
      }
      else if (next_arr[i] == -1 && pred_arr[j] == -1)
      {
        if (availble_energy[i] > distances[i][j])
        {
          availble_energy[j] = availble_energy[i] - distances[i][j];
        }
        else
        {
          availble_energy[j] = max_length - distances[station[i][j]][j];
          is_throught_station[i] = true;
        }
        int end = update_route(j);
        struct route *r = find_route(i);
        next_arr[i] = j;
        pred_arr[j] = i;
        r->end = end;
      }
      else if (next_arr[j] == -1 && pred_arr[i] == -1)
      {
        if (availble_energy[j] > distances[i][j])
        {
          availble_energy[i] = availble_energy[j] - distances[i][j];
        }
        else
        {
          availble_energy[i] = max_length - distances[station[i][j]][i];
          is_throught_station[j] = true;
        }
        int end = update_route(i);
        struct route *r = find_route(j);
        next_arr[j] = i;
        pred_arr[i] = j;
        r->end = end;
      }
    }
    break;
  default:
    return;
  }
}

void clarke_wright()
{
  int m;
  for (m = 0; m < num_saving_distance; m++)
  {
    int i, j, s;
    i = saving_distances[m][0];
    j = saving_distances[m][1];
    s = saving_distances[m][2];
    if (!check_in_the_same_route(i, j))
      merge_route(i, j);
  }
}

void create_capacity_seq(int start, int end)
{
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    if (!temp->is_merged)
    {
      int start = temp->start;
      int end = temp->end;
      int x = start;
      while (x != end)
      {
        if (x == start)
        {
          F_CAP[x] = demands[x];
        }
        else
        {
          F_CAP[x] = F_CAP[pred_arr[x]] + demands[x];
        }
        F_ENERGY[x] = availble_energy[x];
        x = next_arr[x];
      }

      if (x != -1)
      {
        F_CAP[x] = F_CAP[pred_arr[x]] + demands[x];
        F_ENERGY[x] = availble_energy[x];
      }

      x = end;
      while (x != start)
      {
        if (x == end)
        {
          B_CAP[x] = demands[x];
          if (is_throught_station[x])
          {
            B_ENERGY[x] = max_length - distances[station[x][0]][x];
          }
          else
          {
            B_ENERGY[x] = max_length - distances[0][x];
          }
        }
        else
        {
          B_CAP[x] = B_CAP[next_arr[x]] + demands[x];
          if (is_throught_station[x])
          {
            B_ENERGY[x] = max_length - distances[station[x][next_arr[x]]][x];
          }
          else
          {
            B_ENERGY[x] = B_ENERGY[next_arr[x]] - distances[x][next_arr[x]];
          }
        }
        x = pred_arr[x];
      }

      if (x != -1)
      {
        B_CAP[x] = B_CAP[next_arr[x]] + demands[x];
        if (is_throught_station[x])
        {
          B_ENERGY[x] = max_length - distances[station[x][next_arr[x]]][x];
        }
        else
        {
          B_ENERGY[x] = B_ENERGY[next_arr[x]] - distances[x][next_arr[x]];
        }
      }
    }

    temp = temp->next;
  }
}

void show_result()
{
  if (ROUTE == NULL)
  {
    printf("no data");
    return;
  }
  printf("\n");
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    if (!temp->is_merged)
    {
      int start = temp->start;
      int end = temp->end;
      create_capacity_seq(start, end);
      int x = start;
      printf("start: %d - end: %d\n", start, end);
      while (x != end)
      {
        printf("%d ", x);
        x = next_arr[x];
      }
      printf("%d \n", end);

      // show F_CAP
      x = start;
      while (x != end)
      {
        printf("%.1lf ", F_CAP[x]);
        x = next_arr[x];
      }
      printf("%.1lf ", F_CAP[x]);
      printf("\n");
      x = start;
      while (x != end)
      {
        printf("%.1lf ", B_CAP[x]);
        x = next_arr[x];
      }
      printf("%.1lf ", B_CAP[x]);
      printf("\n");
      x = start;
      while (x != end)
      {
        printf("%.1lf ", F_ENERGY[x]);
        x = next_arr[x];
      }
      printf("%.1lf ", F_ENERGY[x]);
      printf("\n");
      x = start;
      while (x != end)
      {
        printf("%.1lf ", B_ENERGY[x]);
        x = next_arr[x];
      }
      printf("%.1lf ", B_ENERGY[x]);
      printf("\n");
    }
    temp = temp->next;
  }
}

int main()
{
  read_file("./A-n32-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n33-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n33-k6.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n34-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n36-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n37-k6.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n38-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n39-k5.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n39-k6.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();

  printf("=============================");
  read_file("./A-n45-k6.vrp");
  prepare_data();
  init();
  clarke_wright();
  show_result();
}
