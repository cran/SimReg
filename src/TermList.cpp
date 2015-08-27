#include "TermList.h"

using namespace Rcpp;
using namespace std;

term_list::term_list(IntegerVector in_term_ids, IntegerVector in_case_ids, int in_num_cases) {
	term_ids = in_term_ids;
	case_ids = in_case_ids;
	num_cases = in_num_cases;
}


