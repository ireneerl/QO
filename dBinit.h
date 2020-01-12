#include iostream;
#include vector;

void getCombination(int data[], int l, int r, int n, int *result, int * idx_combination);
void alloc_db(int n);
void init_database(int n);
std::tuple<int, int> getCost(std::vector<TUPLE*>& table_observed);
void print_database(int n);
void print_tuple(std::vector<TUPLE*>& data);
void print_rel(std::vector<int> const &input);
void print_first_memoization(std::vector<std::vector<int>> relations);
void print_performance(struct timeval begin, struct timeval end);
