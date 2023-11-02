 real compute_log_L(array[] real theta){
    real r = theta[1];
    real w = theta[2];
    real p = theta[3];
    real ans = (1 + w) * ( log(1 + w) - log(w + r) ) + log(r) + log1m(p);
    return(ans);
  }
  array[] real log_marg_prob_s(int s, real r, real w, real p,
                               real epsilon, int max_it, int trunc_cap,
                               real logL, int type){
    array[2] real lprob;                                 
    int n0 = max({s, 1});
    if(type == 1){// Adaptive
      lprob = infiniteAdaptive({r, w, p, s}, epsilon, max_it, logL, n0);
    }else{ // fixed cap
      lprob = finiteSumToCap({r, w, p, s}, trunc_cap, n0);
    }
    return(lprob);
  }
  array[] real compute_likelihood(int D, vector n, array[] int s,
                                  real r, real w, real p,
                                  real lp0, real logL, 
                                  real epsilon, int max_it, int trunc_cap,
                                  int type){
    real n_iters = 0;
    vector[D] log_Prs;
    for(i in 1:D){
      array[2] real temp = log_marg_prob_s(s[i],
                                          r, w, p,
                                          epsilon, max_it, trunc_cap,
                                          logL, type);
      n_iters = n_iters + temp[2];
      log_Prs[i] = temp[1] -  log1m_exp(lp0); 
    } 
    
    return({dot_product(n, log_Prs), n_iters});
  }
  real obj_fun(real x, real b, real w){
    real ans =  pow(b - x/pow(1 + x, 1 + w), 2);
    return(ans);
  }
  real obj_fun_grad(real x, real b, real w){
    real ans = 2*pow(x + 1, -2*w -3) * (w*x -1) *( pow(x+1, w)*(b*x + b) -x );
    return(ans);
  }
  real get_u_Newton(real r, real w, real p, real tol, int max_iter){
    // Inspired by https://rpubs.com/aaronsc32/newton-raphson-method
    real b = r/(w * pow(1 + r/w, w + 1)) * (1-p);
    real l = 0;
    real upr = 1/w; 
    real x0 = l;
    real x1 = 0;
    real ans = 0;
    int iter_counter = 0;
    while(iter_counter < max_iter){
      real deriv =  obj_fun_grad(x0, b, w);
      x1 = x0 - (obj_fun(x0, b, w) / deriv);
      // print("Doing iteration ", iter_counter + 1, " x1 = ", x1);
      if (abs(x1 - x0) < tol){
        ans = x1;
        break;
      }
      iter_counter += 1;
      x0 = x1;
      if(iter_counter >= max_iter){
        print("Exceeded maximum iterations; returning what I have");
        return(ans);
      }
      if(x1 > upr) print("R0=", r, " w=", w, " p=", p);
      if(x1 > upr) reject("Solution is out of bounds.");
    }
    return(ans);
  }
  real prob_s0(real r, real w, real p, real tol, int max_iter) {// Pr(S = 0 | R0, k, psi) via shitty implementation of Newton-Raphson. 
  real a;
  real u;
  real lp;
  a = r/w; // R0/k
  u =  get_u_Newton(r, w , p, tol, max_iter);
  lp = log1p(a) - log(a) + log(u)-log1p(u);// (a + 1)/a * u/(u + 1), see notes
  return(lp);
  }