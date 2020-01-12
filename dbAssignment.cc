#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <sys/time.h>
#include <limits.h>
#include <vector>

#define DEFAULT_TABLE_SIZE 1000
#define DEFAULT_NUM_SPACE 10000
#define MAX_TABLE_NUM 16

struct timeval
  cur_time(void) {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t;
}

typedef struct _RECORD{
  std::vector<int> attr;
}  TUPLE;

typedef struct return_DP {
  int minimum_cost;
  std::vector<int> path;
} DP_Param;

std::vector<std::vector<TUPLE*> > rTables;
std::vector<std::vector<TUPLE*> > tables(MAX_TABLE_NUM);

int numSpace;
int tableSize;
int relationNumber;

void getCombination(int data[], int l, int r, int n, int *result, int * idx_combination){
    if(l==r){
        std::vector <int> tmp;
        for(int i=0; i < n; i++){
            *((result+ (*idx_combination) *n) + i) = data[i];
        }(*idx_combination) = (*idx_combination)  + 1;
    }else{
        for (int i = l; i <=r; i++){
            std::swap(data[l], data[i]);
            getCombination(data, l+1, r, n, result,idx_combination );
            std::swap(data[l], data[i]);
        }
    }
}

void alloc_join(int n, int *result){
  int data[n];
  for (int i = 0 ; i < n ; i++)
    data[i] = i;
  int idx_combination = 0;
  getCombination(data, 0, n-1, n, (int *)result, &idx_combination);
}

TUPLE* alloc_tuple(int j){
  TUPLE *tuple;
  if(!(tuple = (TUPLE *) calloc(1, sizeof(TUPLE))))
    printf("ERROR! (%s): calloc failed\n", __func__);
  //by value
  tuple->attr.push_back(j); //id_key
  tuple->attr.push_back(rand() % numSpace);
  return tuple;
}

void init_database(int n){
  for(int i = 0; i < n; i++){
    for(int j = 0; j < tableSize; j++)
      tables[i].push_back(alloc_tuple(j));
  }
}

//---get cost ----
int V(std::vector<TUPLE*> R){ //-> std::vector<TUPLE*>& R
 auto comp_ = [] ( const TUPLE *lts, const TUPLE *rts ) {return lts->attr[1] == rts->attr[1];};
 auto pred_ = []( const TUPLE *lts, const TUPLE *rts ) {return lts->attr[1] < rts->attr[1];};
 std::sort(R.begin(),R.end(),pred_);
 auto last_ = std::unique(R.begin(), R.end(),comp_);
 R.erase(last_, R.end());
  return R.size();
}

int row_R(std::vector<TUPLE*>& R){
  return R.size();
}

std::tuple<int, int> getCost(std::vector<TUPLE*>& table_observed){
  int row_r = row_R(table_observed);
  int v_r = V(table_observed);
  return {row_r, v_r};
}

//------ show the data
void print_database(int n){
  for(int i = 0 ; i < n; i++){
    std::cout << "Table " << i << '\n';
    for(long unsigned int j = 0; j < tables[i].size();j++){
      printf("[Row#: %ld] %d %d\n", j, tables[i][j]->attr[0], tables[i][j]->attr[1]);
    }
  }
}

void print_tuple(std::vector<TUPLE*>& data){
  for (int i=0; i < data.size(); i ++){
    for (int j = 0 ; j < data[1]->attr.size(); j++){
      std::cout << data[i]->attr[j]<< ";" ;
    }
    std::cout <<std::endl;
  }
}

void print_rel(std::vector<int> const &input){
  std::cout << "relation -> " ;
  for (int i = 0; i < input.size(); i++)
    std::cout << char('A' + input.at(i)) << " ";
  std::cout << '\n';
}
char nth_letter(int n)
{
  n = n+1;
  assert(n >= 1 && n <= 26);
    return "abcdefghijklmnopqrstuvwxyz"[n-1];
}

void print_first_memoization(std::vector<std::vector<int>> relations){
  for(int i=0; i<relations.size(); i++){
    for(int j=0; j<relations[0].size(); j++)
      std::cout << relations[i][j]<< " ";
  std::cout << "" << '\n';
  }
}

void print_performance(struct timeval begin, struct timeval end) {
  std::cout << "printing the performance";
  long diff = (end.tv_sec - begin.tv_sec) * 1000 * 1000
                + (end.tv_usec - begin.tv_usec);
  printf(" lat:%7ld usec\n", diff);
}

void print_join_probability(int relationNumber, int n_rel, int *pCombination){
  for (int i = 0 ; i < n_rel ; i++){
    for (int j = 0 ; j < relationNumber ; j++)
          std::cout << *((pCombination + i * relationNumber) + j) << ":";// std::cout << pCombination[i][j];
    std::cout << std::endl;
  }
}

std::tuple<int, int> brute_force(int relationNumber, int n_rel, int *pCombination){
  int cost, best;
  for (int i = 0 ; i < n_rel ; i++){
    int tmp = 0, l_tuples_1, l_V_1, max_V;
    for (int j = 1 ; j < relationNumber ; j++){
      if(j == 1){
        std::cout << *((pCombination + i * relationNumber) + 0) << ":";// std::cout << pCombination[i][j];
        int relA = *((pCombination + i * relationNumber) + 0);
        int relB = *((pCombination + i * relationNumber) + j);
        auto [l_tuples, l_V] = getCost(tables[relA]);
        auto [r_tuples, r_V] = getCost(tables[relB]);
        max_V = std::max(l_V, r_V);
        tmp = (l_tuples * r_tuples)/max_V;
        l_tuples_1 = tmp;
        l_V_1 = max_V;
        std::cout << "the left " << l_tuples << ";" <<  l_V << "the right " << r_tuples  << ";" <<  r_V << "the max " << max_V << " = " << tmp << "value float -> " << ((l_tuples * r_tuples)/max_V) <<  "or "<< ((20 * 20)/18) <<'\n';

      }else{
        int relB = *((pCombination + i * relationNumber) + j);
        auto [r_tuples, r_V] = getCost(tables[relB]);
        max_V = std::max(l_V_1, r_V);
        tmp += (l_tuples_1 * r_tuples)/max_V;
      }
      // std::cout << *((pCombination + i * relationNumber) + j) << ":";// std::cout << pCombination[i][j];
    }
      std::cout << "i " << i << "standard -> " << tmp << "cost -> " << cost << '\n';
    if(i == 0 || cost > tmp){
      std::cout << "i " << i << "best -> " << tmp << "cost -> " << cost << '\n';
      cost = tmp;
      best = i;
    }else{
      std::cout << "i " << i << "best -> " << tmp << "cost -> " << cost << '\n';
    }
    // std::cout << std::endl;
  }
  return {best, cost};
}


std::tuple<int, int> greedy(int relationNumber, int n_rel, int *pCombination){
  int cost, best, L, R, relA, relB;
  std::vector<int> R_taken(relationNumber, 0);
  for (int layer = 0 ;  layer < relationNumber-1 ; layer ++){
    int tmp = 0, l_tuples_1, l_V_1, max_V;
    if (layer == 0){
      for (int i = 0 ; i < n_rel ; i++){ //loop through all combinations
        std::cout << *((pCombination + i * relationNumber) + layer) << ":";// std::cout << pCombination[i][j];
        relA = *((pCombination + i * relationNumber) + layer);
        relB = *((pCombination + i * relationNumber) + layer+1);
        auto [l_tuples, l_V] = getCost(tables[relA]);
        auto [r_tuples, r_V] = getCost(tables[relB]);
        max_V = std::max(l_V, r_V);
        tmp = (l_tuples * r_tuples)/max_V;
        l_tuples_1 = tmp;
        l_V_1 = max_V;
        if(i == 0 || cost > tmp){
          cost = tmp;
          best = i ;
          R_taken[layer] = relA; //left
          R_taken[layer+1] = relB; //right
        }
      }
    }else{
      for (int i = 0 ; i < n_rel ; i++){ //loop through all combinations
        for (int j = 0 ;  j < R_taken.size() ; j ++){
          if (*((pCombination + i * relationNumber) + j) == R_taken[j]){
            int relB = *((pCombination + i * relationNumber) + layer+1);
            auto [r_tuples, r_V] = getCost(tables[relB]);
            max_V = std::max(l_V_1, r_V);
            tmp = (l_tuples_1 * r_tuples)/max_V;
          }
          if(i == 0 || cost > tmp){
            cost = tmp;
            best = i ;
            R_taken[layer] = relA; //left
            R_taken[layer+1] = relB; //right
          }
        }
      }
    }
  }
  return {best, cost};
}

std::tuple<int, std::vector<int>> selinger(const std::vector<std::vector<int>>& relations, int pos, int joined, std::vector<std::vector<int>>& state, std::vector<int> &solution)
{
  if(joined == ((1 << relations.size()) - 1)){
    solution.push_back(pos);
    return {relations[pos][0], solution}; //the last relation
  }

  if(state[pos][joined] != INT_MAX){
    return {state[pos][joined], solution};
  }
  std::vector<int> temp;
  solution.clear();
  // std::cout << "the length of relation pos -> " << pos << " -> " << state.size() << '\n';
  int best_cost = 1000000;
  int total_relations = relations.size();

  for(int i = 0; i < total_relations; ++i)
    {
      // std::cout << "joined " << (joined & (1 << i)) << " -> " << (1 << i) << '\n';
    if(i == pos || (joined & (1 << i)))
      continue;
    auto [join_cost, sol] = selinger(relations, i, joined | (1 << i), state, solution);
    int cost = (i == total_relations-1)? join_cost : relations[pos][i] + join_cost;
    solution.push_back(pos);

    // std::cout << "solution pos -> " << pos << "joined -> " << joined << '\n';


    if(cost < state[pos][joined]){
      state[pos][joined] = cost;
        if(pos == 0){
          if(cost < best_cost){
            best_cost = cost;
            temp.clear();
            for (int x=0; x<solution.size(); x++) {
              temp.push_back(solution[x]);
              // std::cout << "solution -> " << solution[x]<< '\n';
            }
            std::cout <<  '\n';
          }
        }
      }
    }
    return {state[pos][joined], temp};
}


void DP(int n) {
  std::vector<std::vector<int>> relations;
  int first_join_min = 100000000;
  int relation_joined[2];

  //making first layer memoization
  for(int i=0; i<n; i++){
    std::vector <int> rRow;
    int tmp = 0, max_V = 0;
    for(int j=0; j<n; j++){
      auto [l_tuples, l_V] = getCost(tables[i]);
      auto [r_tuples, r_V] = getCost(tables[j]);
      max_V = std::max(l_V, r_V);
      tmp = (i == j)?  0 : (l_tuples * r_tuples)/max_V;
      std::cout << "the left " << l_tuples <<  l_V << "the right " << r_tuples <<  r_V << "the max " << max_V << " = " << tmp << "value float -> " << ((l_tuples * r_tuples)/max_V) <<  "or "<< ((20 * 20)/18) <<'\n';
      rRow.push_back( tmp );
    }
    relations.push_back(rRow);
  }

  for(int i=0; i<n; i++){
    for (int j=0;j<n;j++){
      if (relations[i][j] < first_join_min && i != j){
        relation_joined[0] = i;
        relation_joined[1] = j;
        first_join_min = relations[i][j];
      }
    }
  }
  std::cout << "joined relation" << relation_joined[0] << relation_joined[1] << '\n';
  print_first_memoization(relations);

  //building the second memoization onward
  std::vector<std::vector<int>> state(relations.size());
  for(auto& neighbors : state)
    neighbors = std::vector<int>((1 << relations.size()) - 1, INT_MAX);
    std::vector<int> vect;
    // vect.push_back(relation_joined[1]);
    auto [join_cost, path] = selinger(relations, relation_joined[0], 1, state, vect);
    std::cout << "minimum tupple generated: " << std::endl << join_cost << std::endl;
    print_rel(path);
}


int main(int argc, char *argv[]){
  tableSize = DEFAULT_TABLE_SIZE;
  numSpace = DEFAULT_NUM_SPACE;
  relationNumber = MAX_TABLE_NUM;
  struct timeval begin, end;

  if (argc >= 2) tableSize = atoi(argv[1]); //total rows
  if (argc >= 3) numSpace = atoi(argv[2]); //the key
  if (argc == 4) relationNumber = atoi(argv[3]); //total relations

  if(relationNumber > MAX_TABLE_NUM || relationNumber < 2){
    printf("Invalid number of tables = %d, either too much or too little\n", relationNumber);
  }else{
    int n_sum = 1; //number of relation combination.
    for (int i = 1 ; i < relationNumber+1 ; i++)
        n_sum *= i;
    std::cout << "total summary" << n_sum << '\n';
    int pCombination[n_sum][relationNumber];


    init_database(relationNumber);
    print_database(relationNumber);
    alloc_join(relationNumber, (int *)pCombination);
    print_join_probability(relationNumber, n_sum, (int *)pCombination);

    begin = cur_time();
    auto [sol_native, tuples] = brute_force(relationNumber, n_sum, (int *)pCombination);
    end = cur_time();
    print_performance(begin, end);
    begin = cur_time();
    auto [sol_greedy, tuples_greedy] = greedy(relationNumber, n_sum, (int *)pCombination);
    end = cur_time();
    print_performance(begin, end);

    std::cout << "native" << sol_native <<  '\n';
    for (int i = 0 ; i < relationNumber ; i++) {std::cout << nth_letter(pCombination[sol_native][i]);}
    std::cout << "\ntotal value = " << tuples << '\n';

    std::cout << "greedy" << sol_native <<  '\n';
    for (int i = 0 ; i < relationNumber ; i++) {std::cout << nth_letter(pCombination[sol_greedy][i]);}
    std::cout << "\ntotal value = " << tuples << '\n';

    begin = cur_time();
    DP(relationNumber);
    end = cur_time();
  }
}
