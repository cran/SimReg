#include <Rcpp.h>
#include <cmath>
#include "Utils.h"

#ifndef INCLUDED_HPO_H
#define INCLUDED_HPO_H

using namespace Rcpp;
using namespace std;

Rcpp::LogicalVector hpo_leaves(Rcpp::LogicalMatrix, Rcpp::IntegerVector);
Rcpp::IntegerVector minimal(Rcpp::LogicalMatrix, Rcpp::IntegerVector);
Rcpp::LogicalMatrix hpo_leaves_matrix(Rcpp::LogicalMatrix, Rcpp::IntegerMatrix);

int count_minimal(Rcpp::LogicalMatrix, int);

void first_combination(Rcpp::IntegerVector item, size_t n);

bool next_combination(Rcpp::IntegerVector item, size_t n, size_t N);

int num_orderings(IntegerVector items_of_each_kind);

int num_leaf_set_reps(int leaves, int num_exc_ancs, int vector_length);

bool next_capped_combo(IntegerVector item, int n, int cap);

#endif
