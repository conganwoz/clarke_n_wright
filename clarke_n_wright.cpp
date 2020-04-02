#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

using namespace std;

// class Route
struct route
{
  int start;
  int end;
  double length;
  double capacity;
  struct route *next;
  bool is_begin_with_station;
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

  printf("data: dimension: %d - num_customer: %d - capacity: %lf - energy: %lf", dimention, num_customer, capacity, energy);

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

bool check_all_visited()
{
  int i;
  for (i = 0; i < num_customer; i++)
  {
    if (!visited[i])
      return false;
  }
  return true;
}

void init()
{
  ROUTE = (struct route *)malloc(sizeof(struct route));
  END_ROUTE = ROUTE;
  for (int i = 1; i < num_customer; i++)
  {
    if (demands[i] <= capacity)
    {
      int best_station = station[0][i];
      double best_dis_station = dis_station[0][i];
      double distant_i = distances[best_station][i];

      if (max_length > 2 * distances[0][i])
      {
        //struct route * temp = (struct route *)malloc(sizeof(struct route));
        availble_energy[i] = max_length - distances[0][i];
        next_arr[i] = -1;
        pred_arr[i] = -1;
        is_interior[i] = false;
        visited[i] = true;
        END_ROUTE->next = (struct route *)malloc(sizeof(struct route));
        ;
        END_ROUTE = END_ROUTE->next;
        END_ROUTE->start = i;
        END_ROUTE->end = i;
        END_ROUTE->length = 2 * distances[0][i];
        END_ROUTE->capacity = demands[i];
        END_ROUTE->next = NULL;
        continue;
      }
      else if (
          (max_length < distances[0][i]) && (distances[best_station][0] < max_length) && (max_length > distances[best_station][i]) && (max_length - distances[best_station][i] > distances[0][i]))
      {
        //struct route * temp = (struct route *)malloc(sizeof(struct route));
        availble_energy[i] = max_length - distances[best_station][i];
        next_arr[i] = -1;
        pred_arr[i] = -1;
        is_interior[i] = false;
        visited[i] = true;
        // temp = temp->next;
        END_ROUTE->next = (struct route *)malloc(sizeof(struct route));
        END_ROUTE = END_ROUTE->next;

        END_ROUTE->start = i;
        END_ROUTE->end = i;
        END_ROUTE->length = best_dis_station + distances[0][i];
        END_ROUTE->capacity = demands[i];
        END_ROUTE->is_begin_with_station = true;
        is_throught_station[i] = false;
        END_ROUTE->next = NULL;
        continue;
      }
      else if (
          (max_length > distances[0][i]) && (max_length - distances[0][i] < distances[0][i]) && (max_length - distances[0][i] > distances[i][best_station]) && (max_length > distances[0][best_station]))
      {
        struct route *temp = (struct route *)malloc(sizeof(struct route));
        availble_energy[i] = max_length - distances[0][i];
        next_arr[i] = -1;
        pred_arr[i] = -1;
        is_interior[i] = false;
        visited[i] = true;
        // temp = temp->next;
        END_ROUTE->next = temp;
        END_ROUTE = temp;
        temp->start = i;
        temp->end = i;
        temp->length = dis_station[0][i] + distances[0][i];
        temp->capacity = demands[i];
        temp->is_begin_with_station = false;
        is_throught_station[i] = true;
        temp->next = NULL;
        continue;
      }
      else if (
          (max_length < distances[0][i]) && (max_length > distances[best_station][0]) && (max_length > distances[i][best_station]) && (max_length - distances[best_station][i] < distances[i][0]) && (max_length - distances[best_station][i] > distances[i][best_station]) && (max_length > distances[0][best_station]))
      {
        struct route *temp = (struct route *)malloc(sizeof(struct route));
        availble_energy[i] = max_length - distances[best_station][i];
        next_arr[i] = -1;
        pred_arr[i] = -1;
        is_interior[i] = false;
        visited[i] = true;
        // temp = temp->next;
        END_ROUTE->next = temp;
        END_ROUTE = temp;
        temp->start = i;
        temp->end = i;
        temp->length = 2 * dis_station[0][i];
        temp->capacity = demands[i];
        temp->is_begin_with_station = true;
        is_throught_station[i] = true;
        temp->next = NULL;
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
    if(x == i) is_found = true;
    if (is_found)
      break;

    tempi = tempi -> next;
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
    if(x == j) is_found = true;
    if (is_found)
      break;

    tempj = tempj -> next;
  }

  if (tempi != NULL && tempj != NULL)
  {
    return tempi->capacity + tempj->capacity < capacity;
  }
  else
    return false;
}

int temp_fix_index[1000];
int pred_available[1000];
int pred_throught_station[1000];

void refresh()
{
  for (int i = 0; i < num_customer; i++)
  {
    temp_fix_index[i] = -1;
    pred_available[i] = -1.0;
    pred_throught_station[i] = -2;
  }
}

void recovery_route()
{
  for (int i = 0; i < num_customer; i++)
  {
    if (pred_available[i] >= 0)
    {
      availble_energy[i] = pred_available[i];
    }
    if (pred_throught_station[i] == -1)
      is_throught_station[i] = false;
    if (pred_throught_station[i] == 1)
      is_throught_station[i] = true;
  }
}

bool fix_route(int i, struct route *r)
{
  refresh();
  int k = 0;
  int end = r->end;
  int x = next_arr[i];
  while (x != -1)
  {
    int pred_x = pred_arr[x];
    double avail_pred_x = availble_energy[pred_arr[x]];
    int best_station = station[pred_x][x];
    double best_dis_station = dis_station[pred_x][x];

    pred_available[x] = availble_energy[x];
    if (is_throught_station[x])
      pred_throught_station[x] = 1;
    else
      pred_throught_station[x] = -1;

    if (availble_energy[pred_x] <= distances[pred_x][x])
    {
      if (availble_energy[pred_x] > distances[pred_x][best_station] && max_length > distances[x][best_station])
      {
        pred_throught_station[pred_x] = is_throught_station[pred_x];
        is_throught_station[pred_x] = true;
        availble_energy[x] = max_length - distances[x][best_station];
      }
      else
      {
        recovery_route();
        return false;
      }
    }
    else
    {
      availble_energy[x] = availble_energy[pred_x] - distances[pred_x][x];
    }
    x = next_arr[x];
  }

  if (availble_energy[end] < distances[end][0])
  {
    if (availble_energy[end] > distances[station[end][0]][end] && max_length > distances[station[end][0]][0])
    {
      pred_throught_station[end] = is_throught_station[end];
      is_throught_station[end] = true;
    }
    else
    {
      recovery_route();
      return false;
    }
  }

  return true;
}

void add_new_node_to_route(int i, int j)
{
  int best_stationij = station[i][j];
  double best_dis_stationij = dis_station[i][j];
  int best_stationj0 = station[j][0];
  double best_dis_stationj0 = dis_station[j][0];

  if (next_arr[i] == -1)
  {
    if (availble_energy[i] > distances[i][j] + distances[j][0])
    {
      next_arr[i] = j;
      next_arr[j] = -1;
      pred_arr[j] = i;
      availble_energy[j] = availble_energy[i] - distances[i][j];
      is_throught_station[i] = false;
      is_throught_station[j] = false;
      update_bound(i, j, true, distances[i][j] + distances[j][0] - distances[i][0], false, 0);
    }
    else if ( // 1 0
        (availble_energy[i] < distances[i][j]) && (availble_energy[i] > distances[best_stationij][i]) && (max_length > distances[j][best_stationij] + distances[j][0]))
    {
      next_arr[i] = j;
      next_arr[j] = -1;
      pred_arr[j] = i;
      availble_energy[j] = max_length - distances[best_stationij][j];
      is_throught_station[i] = true;
      is_throught_station[j] = false;
      update_bound(i, j, true, best_dis_stationij + distances[j][0] - distances[i][0], false, 0);
    }
    else if ( // 0 1
        (availble_energy[i] > distances[i][j]) && (availble_energy[i] - distances[i][j] < distances[j][0]) && (availble_energy[i] - distances[i][j] > distances[best_stationj0][j]) && (max_length > distances[best_stationj0][0]))
    {
      next_arr[i] = j;
      next_arr[j] = -1;
      pred_arr[j] = i;
      availble_energy[j] = availble_energy[i] - distances[i][j];
      is_throught_station[i] = false;
      is_throught_station[j] = true;
      update_bound(i, j, true, distances[i][j] + best_dis_stationj0 - distances[i][0], false, 0);
    }
    else if ( // 1 1
        (availble_energy[i] < distances[i][j]) && (availble_energy[i] > distances[best_stationij][i]) && (max_length > distances[best_stationij][j]) && (max_length - distances[best_stationij][j] < distances[j][0]) && (max_length - distances[best_stationij][j] > distances[j][best_stationj0]) && (max_length > distances[best_stationj0][0]))
    {
      next_arr[i] = j;
      next_arr[j] = -1;
      pred_arr[j] = i;
      availble_energy[j] = max_length - distances[j][best_stationij];
      is_throught_station[i] = true;
      is_throught_station[j] = true;
      update_bound(i, j, true, best_dis_stationij + best_dis_stationj0 - distances[i][0], false, 0);
    }
  }
  else
  {
    if (max_length > distances[0][j] + distances[i][j])
    {
      next_arr[j] = i;
      pred_arr[j] = -1;
      pred_arr[i] = j;
      availble_energy[j] = max_length - distances[j][0];
      availble_energy[i] = availble_energy[j] - distances[i][j];
      is_throught_station[j] = false;
      struct route *r = update_bound(i, j, false, distances[i][j] + distances[j][0] - distances[i][0], false, 0);
      if (!fix_route(i, r))
      {
        pred_arr[i] = -1;
        availble_energy[i] = max_length - distances[i][0];
        update_bound(i, -1, false, distances[i][j] + distances[j][0] - distances[i][0], false, j);
      }
    }
  }
}

struct route *revert_route(int i)
{
  struct route *temp = ROUTE->next;
  bool is_found = false;
  while (temp != NULL && !is_found)
  {
    int start = temp->start;
    int x = start;
    int end = temp->end;
    while (x != end)
    {
      if (x == i)
        is_found = true;
      x = next_arr[x];
    }
    if (end == x)
      is_found = true;

    if (is_found)
    {
      x = end;
      while (x != start)
      {
        if (x == end)
        {
          next_arr[x] = pred_arr[x];
          pred_arr[x] = -1;
          availble_energy[x] = max_length - distances[end][0];
          // temp -> is_begin_with_station = is_throught_station[x];
        }
        else
        {
          int temp = pred_arr[x];
          pred_arr[x] = next_arr[x];
          next_arr[x] = temp;
          availble_energy[x] = availble_energy[pred_arr[x]] - distances[x][pred_arr[x]];
        }
      }

      temp->start = end;
      temp->end = start;
    }

    if (is_found)
      return temp;
    temp = temp->next;
  }
}

struct route *find_route(int i)
{
  struct route *temp = ROUTE;
  while (temp != NULL)
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
    if(x == i) return temp;
    temp = temp->next;
  }

  return NULL;
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
    // tao route khi 2 node chua dc visited
    // buoc init khong the tao route rieng cho moi node nen khong the tao route vs 2 node nay
    ////
    break;
  case 1:
    // merge j vao route neu thoa man
    if (!is_interior[i] && is_validate_new_node(i, j))
    {
      add_new_node_to_route(i, j);
    }
    break;
  case 2:
    add_new_node_to_route(j, i);
    break;
  case 3:
    if (is_valid_merge_route(i, j))
    {
      if (pred_arr[i] == -1 && pred_arr[j] == -1)
      {
        struct route *routej = revert_route(j);
        struct route *routei = find_route(i);
        double pred_avaiable = availble_energy[i];
        if (availble_energy[j] > distances[i][j])
        {
          availble_energy[i] = availble_energy[j] - distances[i][j];
          if (!fix_route(i, routei))
          {
            availble_energy[i] = pred_avaiable;
          }
          else
          {
            next_arr[j] = i;
            pred_arr[i] = j;
            is_throught_station[j] = false;
          }
        }
        else if (availble_energy[j] <= distances[i][j] && availble_energy[j] > distances[station[i][j]][j] && max_length > distances[station[i][j]][i])
        {
          availble_energy[i] = max_length - distances[station[i][j]][i];
          if (!fix_route(i, routei))
          {
            availble_energy[i] = pred_avaiable;
          }
          else
          {
            next_arr[j] = i;
            pred_arr[i] = j;
            is_throught_station[j] = true;
          }
        }
      }
      else if (next_arr[i] == -1 && next_arr[j] == -1)
      {
        struct route *routej = revert_route(j);
        double pred_avaiable = availble_energy[j];
        if (availble_energy[i] > distances[i][j])
        {
          availble_energy[j] = availble_energy[i] - distances[i][j];
          if (!fix_route(j, routej))
          {
            availble_energy[j] = pred_avaiable;
          }
          else
          {
            next_arr[i] = j;
            pred_arr[j] = i;
            is_throught_station[j] = false;
          }
        }
        else if (availble_energy[i] <= distances[i][j] && availble_energy[i] > distances[station[i][j]][i] && max_length > distances[station[i][j]][j])
        {
          availble_energy[j] = max_length - distances[station[i][j]][j];
          if (!fix_route(i, routej))
          {
            availble_energy[j] = pred_avaiable;
          }
          else
          {
            next_arr[i] = j;
            pred_arr[j] = i;
            is_throught_station[i] = true;
          }
        }
      }
      else if (next_arr[i] == -1 && pred_arr[j] == -1)
      {
        struct route *routej = find_route(j);
        double pred_avaiable = availble_energy[j];
        if (availble_energy[i] > distances[i][j])
        {
          availble_energy[j] = availble_energy[i] - distances[i][j];
          if (!fix_route(j, routej))
          {
            availble_energy[j] = pred_avaiable;
          }
          else
          {
            next_arr[i] = j;
            pred_arr[j] = i;
            is_throught_station[j] = false;
          }
        }
        else if (availble_energy[i] <= distances[i][j] && availble_energy[i] > distances[station[i][j]][i] && max_length > distances[station[i][j]][j])
        {
          availble_energy[j] = max_length - distances[station[i][j]][j];
          if (!fix_route(i, routej))
          {
            availble_energy[j] = pred_avaiable;
          }
          else
          {
            next_arr[i] = j;
            pred_arr[j] = i;
            is_throught_station[i] = true;
          }
        }
      }
      else if (next_arr[j] == -1 && pred_arr[i] == -1)
      {
        struct route *routei = find_route(i);
        double pred_avaiable = availble_energy[i];
        if (availble_energy[j] > distances[i][j])
        {
          availble_energy[i] = availble_energy[j] - distances[i][j];
          if (!fix_route(i, routei))
          {
            availble_energy[i] = pred_avaiable;
          }
          else
          {
            next_arr[j] = i;
            pred_arr[i] = j;
            is_throught_station[j] = false;
          }
        }
        else if (availble_energy[j] <= distances[i][j] && availble_energy[j] > distances[station[i][j]][j] && max_length > distances[station[i][j]][i])
        {
          availble_energy[i] = max_length - distances[station[i][j]][i];
          if (!fix_route(i, routei))
          {
            availble_energy[i] = pred_avaiable;
          }
          else
          {
            next_arr[j] = i;
            pred_arr[i] = j;
            is_throught_station[j] = true;
          }
        }
      }
    }
  }
}

void clarkewright()
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

void show_result()
{
  if (ROUTE == NULL)
  {
    printf("no data");
    return;
  }
  struct route *temp = ROUTE->next;
  while (temp != NULL)
  {
    int start = temp->start;
    int end = temp->end;
    int x = start;
    printf("start: %d - end: %d\n", start, end);
    while (x != end)
    {
      printf("%d ", x);
      x = next_arr[x];
    }
    printf("%d \n", end);
    temp = temp->next;
  }
}

// int main()
// {
//   read_file("./A-n36-k5.vrp");
//   prepare_data();
//   init();

//   clarkewright();

//   show_result();
//   return -1;
// }
