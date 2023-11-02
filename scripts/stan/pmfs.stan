real size_lpmf(int x, real r, real omega) {
  real a = r/omega;
  real b = a/(1 + a)^(1 + omega);
  real lp = log(1+a)-log(a) + lchoose(omega*x + x - 2, omega*x -1) -log(x) + x*log(b);
  return lp;
}
real lbin(real x, int n, real p){
  if(n < x){
    return negative_infinity();
  }else{
    return lchoose(n, x) + x*log(p) + (n-x)*log1m(p);  
  }
}
real logFunction(int k, array[] real pars){
  real r = pars[1];
  real w = pars[2];
  real p = pars[3];
  real s = pars[4];
  real ans =  lbin(s, k, p) + size_lpmf(k | r, w);
  return(ans);
}
