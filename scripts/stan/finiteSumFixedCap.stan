array[] real finiteSumToCap(array[] real p, int Nmax, int n0) {
  vector[Nmax] storeVal;
  int k = 1;
  int xx = n0;
  storeVal[k] = logFunction(xx, p);
  while(k < Nmax){
    k = k + 1;
    xx = xx + 1;
    storeVal[k] = logFunction(xx, p);
  }
  return {log_sum_exp(sort_asc(storeVal)), 1. * Nmax};
}