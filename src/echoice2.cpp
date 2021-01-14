#include "RcppArmadillo.h"
#include <omp.h>

using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

//global timing objects
time_t itime;
time_t itime2;


void startTimer() {
  itime2 = time(NULL);
  Rcout << " Computation in progress \n";
}

void infoTimer(int rep, int R) {
  time_t ctime = time(NULL);    
  
  //time to end
  double timetoend = difftime(ctime, itime2) / 60.0 * (R - rep - 1) / (rep+1);
  
  //percent done, typecast
  double perc1=(double)rep/(double)R;
  
  //round percent done
  int perc = floor(perc1*100+0.5);
  
  //overwrite output with current status
  REprintf("\r");
  REprintf("Computing (%i percent), ETA: %.2f min.", perc, timetoend);
  
}


// utils
void startMcmcTimer() {
  itime = time(NULL);
  Rcout << " MCMC in progress \n";
}

void endMcmcTimer() {
  time_t ctime = time(NULL);
  REprintf("\n MCMC complete\n");
  Rprintf(" Total Time Elapsed: %.2f minutes\n", difftime(ctime, itime) / 60.0);     
  itime = 0;
}

void infoMcmcTimer(int rep, int R){
  time_t ctime = time(NULL);    
  
  //time to end
  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
  
  //percent done, typecast
  double perc1=(double)rep/(double)R;
  
  //round percent done
  int perc = floor(perc1*100+0.5);
  
  //overwrite output with current status
  REprintf("\r");
  REprintf("Iteration: %i of %i (%i percent), ETA: %.2f min.", rep, R, perc, timetoend);
}

void infoMcmcTimerRRLL(int rep, int R, double RejectionRate, double LogLL) {
  time_t ctime = time(NULL);    
  
  //time to end
  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
  
  //percent done, typecast
  double perc1=(double)rep/(double)R;
  
  //round percent done
  int perc = floor(perc1*100+0.5);
  
  //round RR
  int RR  = floor(RejectionRate*100+0.5);
  
  //overwrite output with current status
  REprintf("\r");
  REprintf("Iteration: %i of %i (%i percent), ETA: %.2f min., RR: %i, LogLL: %.1f", rep, R, perc, timetoend, RR, LogLL);
  
}

// constants
static double const log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
vec revd(int n=100, double loc=0, double scale=1){
  return(loc-scale*log(-log(runif(n))));
}


// [[Rcpp::export]]
vec revdx(vec locs, 
          vec scales){
  
  int n = locs.n_elem;
  vec out(n);
  // 
  // for(int i = 0; i < n; i++) {
  //   out(i) = loc-scale*log(-log(runif(1)));
  // }
  
  out = locs-scales%log(-log(randu(n) ));
  
  return(out);
}


//' @export
// [[Rcpp::export]]
vec revd0(int n, double scale){
  
  vec out(n);
  out = -scale*log(-log(randu(n) ));
  
  return(out);
}

//multinomial draw (~rmultinom(1,1,probs))
vec rmuno(vec const& probs){
  int n = probs.n_elem;
  vec out(n);
  out.fill(0);
  
  int k = sum(as_scalar(randu(1))>cumsum(probs));
  if(k<n){
    out(k)=1;
  }
  return(out);
}

//multinomial draw - index
int rmuno2(vec const& probs){
  int n = probs.n_elem;
  vec out(n);
  out.fill(0);
  
  int k = sum(as_scalar(randu(1))>cumsum(probs));
  return(k);
}

//upper lvl
void ULwishart(double nu, 
               mat const& V,
               mat& SIGMA,
               mat& CI){
  
  int m = V.n_rows;
  mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); 
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; 
    }}
  
  //output
  CI    = solve(trimatu(trans(T)*chol(V)),eye(m,m)); 
  SIGMA = CI * trans(CI);
}


void ULreg(mat const& Y, 
           mat const& X, 
           mat const& Bbar, 
           mat const& A, 
           double nu, 
           mat const& V,
           mat& MU,
           mat& SIGMA,
           mat& L){
  
  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  mat CI(m,m);
  
  //first draw Sigma
  mat RA = chol(A);
  mat W = join_cols(X, RA);
  mat Z = join_cols(Y, RA*Bbar);
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); 
  mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  mat E = Z-W*Btilde;
  mat S = trans(E)*E;
  
  // compute the inverse of V+S
  mat ucholinv  = solve(trimatu(chol(V+S)), eye(m,m));
  mat VSinv     = ucholinv*trans(ucholinv);
  
  ULwishart(nu+n, VSinv,
            SIGMA, CI);
  
  //output
  MU = Btilde + IR*randn(k,m)*trans(CI);
  L  = chol(SIGMA);
}


//Normal-Normal (univariate)
void ULnormnorm(double& prior_mean, double& prior_sd,
                vec x, 
                double mu_0, double nu, double alph, double bet){
  
  int n = x.n_elem;
  double xbar=mean(x);
  
  prior_mean=(nu*mu_0+n*xbar)/(nu+n)  +  randn(1)[0]/sqrt(nu+n);
  
  prior_sd= sqrt(1/randg(distr_param(alph+n/2,
                                1/(bet+0.5*(sum(pow(x-xbar,2)))+n*nu/(nu+n)*0.5*pow(xbar-mu_0,2) )  )));
}



//auto-tuner to achieve target rejection rate window
//this is a rather basic auto-tuner
void mh_tuner(vec tunes, vec rrs){
  int n = rrs.size();
  for(int ii=0; ii<n; ii++){
    //tuning
    if(rrs(ii)>0.8){
      tunes(ii)=tunes(ii)-tunes(ii)*.1;
    }else if(rrs(ii)<0.6){
      tunes(ii)=tunes(ii)+tunes(ii)*.1;
    }
    //guardrails
    if(tunes(ii)<0.0001){
      tunes(ii)=0.001;
    }
    
    if(tunes(ii)>4){
      tunes(ii)=4;
    }
  }
}


//crossprod functions for more readable code
mat crprod(arma::mat const& m1){
  return(arma::trans(m1)*m1);
}

mat tcrprod(arma::mat const& m1){
  return(m1*arma::trans(m1));
}


//normal density
// [[Rcpp::export]]
double lndMvnc(arma::vec const& x, 
               arma::vec const& mu, 
               arma::mat const& L){
  
  arma::mat rooti = arma::trans(arma::inv(trimatu(L)));
  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log2pi -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}












///////////////////////////////////////////////
// DD compensatory
///////////////////////////////////////////////


//log-likelihood
//' @export
// [[Rcpp::export]]
double ddl(arma::vec const& theta, 
           arma::uvec const& nalts,
           arma::vec const& X, 
           arma::vec const& P, 
           arma::mat const& A, 
           int ntask, 
           int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-2));
  double beta_p  = exp(theta(p-1));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    double denom=1;
    double aby=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      ab+= (-beta_p*p);
      denom+=exp(ab);
      
      if(x>0){
        aby+=ab;
      }
      
      xpicker+=1; //move up index for X,A,P
    }
    
    //add to LL
    ll+=(aby-log(denom));
    
  }
  return(ll);
}



vec ddLL(mat const&Theta,
         vec const& XX, 
         vec const& PP,
         mat const& AA,
         uvec const& nalts,
         ivec const& ntasks,  
         ivec const& xfr,  
         ivec const& xto,  
         ivec const& lfr,  
         ivec const& lto,
         int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= 
      ddl(Theta.col(n),
          nalts(span(lfr(n),lto(n))),
          XX(span(xfr(n),xto(n))), 
          PP(span(xfr(n),xto(n))), 
          AA(span(xfr(n),xto(n)),span::all), 
          ntasks(n), p);
  }
  
  return(ll_olds);
}


//' @export
// [[Rcpp::export]]
mat ddLLs(cube const&THETAS,
          vec const& XX, 
          vec const& PP,
          mat const& AA,
          uvec const& nalts,
          ivec const& ntasks,  
          ivec const& xfr,  
          ivec const& xto,  
          ivec const& lfr,  
          ivec const& lto,
          int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);
  
  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      ddLL(THETAS.slice(r),
           XX, 
           PP,
           AA,
           nalts,
           ntasks,  
           xfr,  
           xto,  
           lfr,  
           lto,
           p,
           N, cores);
  }
  
  return(ll_olds);
}




//generate one draw 

//i-level draws RWMH
void draw_dd_RWMH( arma::vec& ll_olds,    // vector of current log-likelihoods
                   arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                   arma::mat& theta_temp,    // container of current betas for all i
                   vec const& XX,             // data
                   vec const& PP,
                   mat const& AA,
                   uvec const& nalts,
                   ivec const& ntasks,  
                   ivec const& xfr,ivec const& xto,  
                   ivec const& lfr,ivec const& lto,
                   int p, int N, 
                   arma::vec const& mu,  // upper level mean
                   arma::mat const& L,   // upper level chol(sigma)
                   arma::vec& stay,      // rejection tracker, used for tuning
                   arma::vec& tunes,     // i-level tuning parameters
                   int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    //eval
    llnew = ddl(theta_cand,
                nalts(span(lfr(n),lto(n))),
                XX(span(xfr(n),xto(n))), 
                PP(span(xfr(n),xto(n))), 
                AA(span(xfr(n),xto(n)),span::all), 
                ntasks(n), p );
    
    lpnew = lndMvnc(theta_cand, mu, L);
    
    //A-R
    double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
    
    if(ldiff > log(randu(1)[0])){
      theta_temp.col(n)= theta_cand;
      ll_olds(n)       = llnew;
      lp_olds(n)       = lpnew;
    }else{
      stay(n)+=1;
    }
    
  }
}


// entire loop

//' @export
// [[Rcpp::export]]
List loop_dd_RWMH(  vec const& XX, 
                    vec const& PP,
                    mat const& AA,
                    uvec const& nalts,
                    ivec const& ntasks,  
                    ivec const& xfr,  
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    int p, int N,
                    int R, int keep, // MCMC parameters draws and keep interval
                    mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                    int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                    int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= 
      ddl(theta_temp.col(n),
          nalts(span(lfr(n),lto(n))),
          XX(span(xfr(n),xto(n))), 
          PP(span(xfr(n),xto(n))), 
          AA(span(xfr(n),xto(n)),span::all), 
          ntasks(n), p );
  }
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    draw_dd_RWMH(ll_olds,            // ll for current betas
                 lp_olds,
                 theta_temp,          // container of current betas for all i
                 XX,
                 PP,
                 AA,
                 nalts,
                 ntasks,
                 xfr,xto,lfr,lto,
                 p,N,
                 vectorise(MU),      // Mean Prior
                 Lprior,             // VarCov Prior (chol)
                 stay,               // tracking rejections
                 tunes,              // i-level tuning parameters
                 cores);       
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}




///////////////////////////////////////////////
// Discrete Demand with conjunctive screening rules
///////////////////////////////////////////////

//update screening prior
void drawdelta(vec& delta, 
               imat const& tauis, 
               int K, int N, 
               int cores, 
               double a0=1, double b0=1){
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int kk=0; kk<K; kk++){
    double tsum=sum(tauis.row(kk));
    delta(kk)=Rf_rbeta(tsum+a0, N-tsum+b0);
  }
}

//log-likelihood

//' @export
// [[Rcpp::export]]
double ddlsr(arma::vec const& theta, 
             arma::ivec const& taui,
             arma::uvec const& nalts,
             arma::vec const& X, 
             arma::vec const& P, 
             arma::mat const& A, 
             arma::mat const& Afull, 
             int ntask, 
             int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-2));
  double beta_p  = exp(theta(p-1));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    double denom=1;
    double aby=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      ab+=(-beta_p*p);
      
      //screening
      if(as_scalar(Afull.row(xpicker)*taui)>(0.01)){
        
      }else{
        denom+=exp(ab);
      }
      
      //chosen product
      if(x>0){
        aby+=ab;
      }
      
      xpicker+=1; //move up index for X,A,P
    }
    
    //add to LL
    ll+=(aby-log(denom));
  }
  return(ll);
}



//update 'betas' of DD model with conjunctive screening
void draw_ddsr_RWMH( arma::vec& ll_olds,       // vector of current log-likelihoods
                     arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                     arma::mat& theta_temp,    // container of current betas for all i
                       arma::imat const& tauis, 
                       arma::imat const& tauconsts,
                       vec const& XX,            // data
                       vec const& PP,
                       mat const& AA,
                       mat const& AAf,
                       uvec const& nalts,
                       ivec const& ntasks,  
                       ivec const& xfr,ivec const& xto,  
                       ivec const& lfr,ivec const& lto,
                       int p, int N, 
                       arma::vec const& mu,  // upper level mean
                       arma::mat const& L,   // upper level chol(sigma)
                         arma::vec& stay,      // rejection tracker, used for tuning
                         arma::vec& tunes,     // i-level tuning parameters
                       int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    //eval
    llnew = ddlsr(theta_cand,
                  tauis.col(n),
                  nalts(span(lfr(n),lto(n))),
                  XX(span(xfr(n),xto(n))), 
                  PP(span(xfr(n),xto(n))), 
                  AA(span(xfr(n),xto(n)),span::all), 
                  AAf(span(xfr(n),xto(n)),span::all), 
                  ntasks(n), p );
    
    lpnew = lndMvnc(theta_cand, mu, L);
    
    //A-R
    double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
    
    if(ldiff > log(randu(1)[0])){
      theta_temp.col(n)= theta_cand;
      ll_olds(n)       = llnew;
      lp_olds(n)       = lpnew;
    }else{
      stay(n)+=1;
    }
    
  }
}



//update screening for all individuals
void draw_dd_tau(  arma::vec& ll_olds,
                   imat& tauis,
                   arma::mat const& theta_temp,
                   imat const& tauconst,
                   mat const& delta,
                   vec const& XX, 
                   vec const& PP,
                   mat const& AA,
                   mat const& AAf,
                   uvec const& nalts,
                   ivec const& ntasks,  
                   ivec const& xfr,  
                   ivec const& xto,  
                   ivec const& lfr,  
                   ivec const& lto,
                   int p, int N, int cores){
  
  int K=tauconst.n_rows;
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    double ll0;
    double ll1;
    for(int kk=0; kk<K; kk++){
      
      if(tauconst(kk,n)==1){
      
        if(tauis(kk,n)==1){
          
          ll1=ll_olds(n);
            
          ivec tauk0=tauis.col(n);
          tauk0(kk)=0;
          ll0 = ddlsr( theta_temp.col(n),
                       tauk0,
                       nalts(span(lfr(n),lto(n))),
                       XX(span(xfr(n),xto(n))), 
                       PP(span(xfr(n),xto(n))), 
                       AA(span(xfr(n),xto(n)),span::all), 
                       AAf(span(xfr(n),xto(n)),span::all), 
                       ntasks(n), p );
            
        }else{
          ll0=ll_olds(n);
          
          ivec tauk1=tauis.col(n);
          tauk1(kk)=1;
          ll1 = ddlsr( theta_temp.col(n),
                       tauk1,
                       nalts(span(lfr(n),lto(n))),
                       XX(span(xfr(n),xto(n))), 
                       PP(span(xfr(n),xto(n))), 
                       AA(span(xfr(n),xto(n)),span::all), 
                       AAf(span(xfr(n),xto(n)),span::all), 
                       ntasks(n), p );
        }
          
        double probsc = (exp(ll1) * (delta(kk))) / 
               (exp(ll1) * (delta(kk)) + exp(ll0)*(1-delta(kk)  )); 
          
        tauis(kk,n)= Rf_rbinom( 1, probsc );
          
        if(tauis(kk,n)==1){
          ll_olds(n)=ll1;
        }else{
          ll_olds(n)=ll0;
        }      
          
      }
    }//k=loop
  }//nloop
}




//entire loop of DD conjunctive screening model
//' @export
// [[Rcpp::export]]
List loop_ddrs_RWMH(  vec const& XX, 
                      vec const& PP,
                      mat const& AA,
                      mat const& AAf,
                      imat const& tauconst,
                      uvec const& nalts,
                      ivec const& ntasks,  
                      ivec const& xfr,  
                      ivec const& xto,  
                      ivec const& lfr,  
                      ivec const& lto,
                      int p, int N,
                      int R, int keep, // MCMC parameters draws and keep interval
                      mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                      int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                      int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  int K=tauconst.n_rows;
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  imat tauis   = tauconst;
  vec delta(K); delta.fill(0.5);
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    ll_olds(n)= 
      ddlsr(theta_temp.col(n),
            tauis.col(n),
            nalts(span(lfr(n),lto(n))),
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            AAf(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  icube tauDraw(K,N,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  mat  deltaDraw(Rk,K);  
  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  omp_set_num_threads(cores);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    //theta
    draw_ddsr_RWMH(ll_olds,            // ll for current betas
                   lp_olds,
                   theta_temp,          // container of current betas for all i
                   tauis,
                   tauconst,
                   XX,
                   PP,
                   AA,
                   AAf,
                   nalts,
                   ntasks,
                   xfr,xto,lfr,lto,
                   p,N,
                   vectorise(MU),      // Mean Prior
                   Lprior,             // VarCov Prior (chol)
                   stay,               // tracking rejections
                   tunes,              // i-level tuning parameters
                   cores);       
    
    if(ir>1000){
      //tau
      draw_dd_tau(ll_olds,
                  tauis,
                  theta_temp,
                  tauconst,
                  delta,
                  XX, 
                  PP,
                  AA,
                  AAf,
                  nalts,
                  ntasks,  
                  xfr, xto, lfr,lto,
                  p, N, 
                  cores);
      
      // delta_tau
      drawdelta(delta, tauis, K, N, cores);
      
      //update LL
// #pragma omp parallel for schedule(static)
//       for(int n=0; n<N; n++){
//         ll_olds(n)= 
//           ddlsr(theta_temp.col(n),
//                 tauis.col(n),
//                 nalts(span(lfr(n),lto(n))),
//                 XX(span(xfr(n),xto(n))), 
//                 PP(span(xfr(n),xto(n))), 
//                 AA(span(xfr(n),xto(n)),span::all), 
//                 AAf(span(xfr(n),xto(n)),span::all), 
//                 ntasks(n), p );
//       }
      
    }
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      tauDraw.slice(mkeep)        = tauis;
      
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      deltaDraw.row(mkeep)        = trans(delta);
      
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds); //log-post not complete yet
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("tauDraw")        = tauDraw,
    Named("deltaDraw")      = deltaDraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}




//Log-Likelihood for several respondents and draws for DD conjunctive screening model
vec ddsrLL(mat const& Theta,
           imat const& tauis,
           vec const& XX, 
           vec const& PP,
           mat const& AA,
           mat const& AAf,
           uvec const& nalts,
           ivec const& ntasks,  
           ivec const& xfr,  
           ivec const& xto,  
           ivec const& lfr,  
           ivec const& lto,
           int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= 
      ddlsr(Theta.col(n),
            tauis.col(n),
            nalts(span(lfr(n),lto(n))),
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all),
            AAf(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p);
  }
  
  return(ll_olds);
}


//' @export
// [[Rcpp::export]]
mat ddsrLLs(cube const&THETAS,
            icube const&TAUIS,
            vec const& XX, 
            vec const& PP,
            mat const& AA,
            mat const& AAf,
            uvec const& nalts,
            ivec const& ntasks,  
            ivec const& xfr,  
            ivec const& xto,  
            ivec const& lfr,  
            ivec const& lto,
            int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);
  
  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      ddsrLL(THETAS.slice(r),
             TAUIS.slice(r),
             XX, 
             PP,
             AA,
             AAf,
             nalts,
             ntasks,  
             xfr,  
             xto,  
             lfr,  
             lto,
             p,
             N, cores);
  }
  
  return(ll_olds);
}





///////////////////////////////////////////////
// Discrete Demand - Screening w/price
///////////////////////////////////////////////


//' @export
// [[Rcpp::export]]
double ddlsrpr(arma::vec const& theta, 
               arma::ivec const& taui,
               double tau_pr,
               arma::uvec const& nalts,
               arma::vec const& X, 
               arma::vec const& P, 
               arma::mat const& A, 
               arma::mat const& Afull, 
               int ntask, 
               int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-2));
  double beta_p  = exp(theta(p-1));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    double denom=1;
    double aby=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      ab+=(-beta_p*p);
      
      //screening
      if(as_scalar(Afull.row(xpicker)*taui)>(0.01)){
        
      }else{
        if(p<=exp(tau_pr)){
          denom+=exp(ab);
        }
      }
      
      //chosen product
      if(x>0){
        aby+=ab;
      }
      
      xpicker+=1; //move up index for X,A,P
    }
    
    //add to LL
    ll+=(aby-log(denom));
  }
  return(ll);
}




void draw_ddsrpr_RWMH( arma::vec& ll_olds,       // vector of current log-likelihoods
                       arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                       arma::mat& theta_temp,    // container of current betas for all i
                       arma::imat const& tauis, 
                       vec const& tau_prs,
                       arma::mat const& tauconsts,
                       vec const& XX,            // data
                       vec const& PP,
                       mat const& AA,
                       mat const& AAf,
                       uvec const& nalts,
                       ivec const& ntasks,  
                       ivec const& xfr,ivec const& xto,  
                       ivec const& lfr,ivec const& lto,
                       int p, int N, 
                       arma::vec const& mu,  // upper level mean
                       arma::mat const& L,   // upper level chol(sigma)
                       arma::vec& stay,      // rejection tracker, used for tuning
                       arma::vec& tunes,     // i-level tuning parameters
                       int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    //eval
    llnew = ddlsrpr(theta_cand,
                    tauis.col(n),
                    tau_prs(n),
                    nalts(span(lfr(n),lto(n))),
                    XX(span(xfr(n),xto(n))), 
                    PP(span(xfr(n),xto(n))), 
                    AA(span(xfr(n),xto(n)),span::all), 
                    AAf(span(xfr(n),xto(n)),span::all), 
                    ntasks(n), p );
    
    lpnew = lndMvnc(theta_cand, mu, L);
    
    //A-R
    double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
    
    if(ldiff > log(randu(1)[0])){
      theta_temp.col(n)= theta_cand;
      ll_olds(n)       = llnew;
      lp_olds(n)       = lpnew;
    }else{
      stay(n)+=1;
    }
    
  }
}




void draw_dd_tauipr(arma::vec& ll_olds,
                    imat& tauis,
                    arma::mat const& theta_temp,
                    vec const& tau_prs,
                    imat const& tauconst,
                    mat const& delta,
                    vec const& XX, 
                    vec const& PP,
                    mat const& AA,
                    mat const& AAf,
                    uvec const& nalts,
                    ivec const& ntasks,  
                    ivec const& xfr,  
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    int p, int N, int cores){
  
  int K=tauconst.n_rows;
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    double ll0;
    double ll1;
    
    for(int kk=0; kk<K; kk++){
      
      if(tauconst(kk,n)==1){
        
        if(tauis(kk,n)==1){
          
          
          ll1=ll_olds(n);
          
          ivec tauk0=tauis.col(n);
          tauk0(kk)=0;
          ll0 = ddlsrpr( theta_temp.col(n),
                         tauk0,
                         tau_prs(n),
                         nalts(span(lfr(n),lto(n))),
                         XX(span(xfr(n),xto(n))), 
                         PP(span(xfr(n),xto(n))), 
                         AA(span(xfr(n),xto(n)),span::all), 
                         AAf(span(xfr(n),xto(n)),span::all), 
                         ntasks(n), p );
          
        }else{
          ll0=ll_olds(n);
          
          ivec tauk1=tauis.col(n);
          tauk1(kk)=1;
          ll1 = ddlsrpr( theta_temp.col(n),
                         tauk1,
                         tau_prs(n),
                         nalts(span(lfr(n),lto(n))),
                         XX(span(xfr(n),xto(n))), 
                         PP(span(xfr(n),xto(n))), 
                         AA(span(xfr(n),xto(n)),span::all), 
                         AAf(span(xfr(n),xto(n)),span::all), 
                         ntasks(n), p );
        }
        
        double probsc = (exp(ll1) * (delta(kk))) / 
          (exp(ll1) * (delta(kk)) + exp(ll0)*(1-delta(kk)  )); 
        
        tauis(kk,n)= Rf_rbinom( 1, probsc );
        
        if(tauis(kk,n)==1){
          ll_olds(n)=ll1;
        }else{
          ll_olds(n)=ll0;
        }
        
      }
    }//k=loop
  }//nloop
}


void draw_dd_taupr( vec& ll_olds,
                    imat const& tauis,
                    arma::mat const& theta_temp,
                    vec& tau_prs,
                    vec const& maxpaids,
                    double const& pr_mean,
                    double const& pr_sd,
                    vec& stay,
                    vec const& pricetune,
                    imat const& tauconst,
                    vec const& XX, 
                    vec const& PP,
                    mat const& AA,
                    mat const& AAf,
                    uvec const& nalts,
                    ivec const& ntasks,  
                    ivec const& xfr,  
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    int p, int N, int cores){
  
  //int K=tauconst.n_rows;
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    
    //candidate
    double tau_pr_cand = tau_prs(n) + pricetune(n)*randn(1)[0];
    
    //check
    if(tau_pr_cand>log(maxpaids(n))){
      
      
      double llnew = ddlsrpr(theta_temp.col(n),
                             tauis.col(n),
                             tau_pr_cand,
                             nalts(span(lfr(n),lto(n))),
                             XX(span(xfr(n),xto(n))), 
                             PP(span(xfr(n),xto(n))), 
                             AA(span(xfr(n),xto(n)),span::all), 
                             AAf(span(xfr(n),xto(n)),span::all),
                             ntasks(n), p);
      
      //A-R
      double ldiff = 
        llnew      + log_normpdf(tau_pr_cand, pr_mean, pr_sd ) - 
        ll_olds(n) - log_normpdf(tau_prs(n) , pr_mean, pr_sd );
      
      
      if(ldiff > log(randu(1)[0])){
        tau_prs(n)  = tau_pr_cand;
        ll_olds(n)  = llnew;
      }else{
        stay(n)+=1;
      }
    }else{
      stay(n)+=1;
    }
    
    
    
    
  }//nloop
}


//' @export
// [[Rcpp::export]]
List loop_ddrspr_RWMH(  vec const& XX, 
                        vec const& PP,
                        mat const& AA,
                        mat const& AAf,
                        imat const& tauconst,
                        uvec const& nalts,
                        ivec const& ntasks,  
                        ivec const& xfr,  
                        ivec const& xto,  
                        ivec const& lfr,  
                        ivec const& lto,
                        int p, int N,
                        int R, int keep, // MCMC parameters draws and keep interval
                        mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                        int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                        int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  
  vec maxpaids(N);
  for(int n=0; n<N; n++){
    maxpaids(n) = max(XX(span(xfr(n),xto(n)))%PP(span(xfr(n),xto(n))));
  }
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  int K=tauconst.n_rows;
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  imat tauis   = tauconst;
  vec delta(K); delta.fill(0.5);
  vec tau_prs = maxpaids*1.1;
  double pr_mean = mean(tau_prs);
  double pr_sd = 1;
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  arma::vec stay_prscr(N);
  stay_prscr.fill(0);
  arma::vec pricetunes(N);
  pricetunes.fill(.1);
  
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= 
      ddlsrpr(theta_temp.col(n),
              tauis.col(n),
              tau_prs(n),
              nalts(span(lfr(n),lto(n))),
              XX(span(xfr(n),xto(n))), 
              PP(span(xfr(n),xto(n))), 
              AA(span(xfr(n),xto(n)),span::all), 
              AAf(span(xfr(n),xto(n)),span::all), 
              ntasks(n), p );
  }
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  icube tauDraw(K,N,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  mat  deltaDraw(Rk,K);  
  mat  pricescreenPriorDraw(Rk,2);
  mat  tau_pr_draw(N,Rk);
  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  omp_set_num_threads(cores);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    //theta
    draw_ddsr_RWMH(ll_olds,            // ll for current betas
                   lp_olds,
                   theta_temp,          // container of current betas for all i
                   tauis,
                   tauconst,
                   XX,
                   PP,
                   AA,
                   AAf,
                   nalts,
                   ntasks,
                   xfr,xto,lfr,lto,
                   p,N,
                   vectorise(MU),      // Mean Prior
                   Lprior,             // VarCov Prior (chol)
                   stay,               // tracking rejections
                   tunes,              // i-level tuning parameters
                   cores);       
    
    if(ir>1000){
      //tau
      draw_dd_tauipr(ll_olds,
                     tauis,
                  theta_temp,
                  tau_prs, //
                  tauconst,
                  delta,
                  XX, 
                  PP,
                  AA,
                  AAf,
                  nalts,
                  ntasks,  
                  xfr, xto, lfr,lto,
                  p, N, 
                  cores);
      

      
      
      // delta_tau
      drawdelta(delta, tauis, K, N, cores);
      
      //update LL
// #pragma omp parallel for schedule(static)
//       for(int n=0; n<N; n++){
//         ll_olds(n)= 
//           ddlsrpr(theta_temp.col(n),
//                   tauis.col(n),
//                   tau_prs(n),
//                   nalts(span(lfr(n),lto(n))),
//                   XX(span(xfr(n),xto(n))), 
//                   PP(span(xfr(n),xto(n))), 
//                   AA(span(xfr(n),xto(n)),span::all), 
//                   AAf(span(xfr(n),xto(n)),span::all), 
//                   ntasks(n), p );
       // }
      
      
      
      draw_dd_taupr( ll_olds,
                     tauis,
                     theta_temp,
                     tau_prs,
                     maxpaids,
                     pr_mean,
                     pr_sd,
                     stay_prscr,
                     pricetunes,
                     tauconst,
                     XX,
                     PP,
                     AA,
                     AAf,
                     nalts,
                     ntasks,
                     xfr,xto,lfr,lto,
                     p,N,
                     cores);
      
      //price screening upper level mu_0, nu, alph, bet
      ULnormnorm(pr_mean, pr_sd,
                 tau_prs,
                 0.0, 0.01, 3.0, 3.0);  
      
    }
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      tauDraw.slice(mkeep)        = tauis;
      
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      deltaDraw.row(mkeep)        = trans(delta);
      
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
      tau_pr_draw.col(mkeep)        = tau_prs;
      pricescreenPriorDraw(mkeep,0) = pr_mean;
      pricescreenPriorDraw(mkeep,1) = pr_sd;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("tauDraw")        = tauDraw,
    Named("deltaDraw")      = deltaDraw,
    Named("tau_pr_draw")    = tau_pr_draw,
    Named("prscreenMuSigDraw")  = pricescreenPriorDraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}








///////////////////////////////////////////////
// DD Demand
///////////////////////////////////////////////


//' @export
//[[Rcpp::export]]
List dddem(vec const& PP,
           mat const& AA,
           uvec const& nalts,
           uvec const& tlens,
           ivec const& ntasks,  
           ivec const& xfr,
           ivec const& xto,  
           ivec const& lfr,  
           ivec const& lto,
           cube const& thetaDraw,
           int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);

  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs = PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-2));
        double beta_p = exp(theta(p-1));
        
        arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
        arma::vec pr = exp(ab)/(1+sum(exp(ab)));
        int ch = sum(as_scalar(randu(1))>cumsum(pr));

        //if not outside good, choose inside good
        if(ch<nalt){
          demcontainer(ch,ir)=1;
        }
        
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}




//' @export
//[[Rcpp::export]]
List ddsrdem(vec const& PP,
             mat const& AA,
             mat const& AAf,
             uvec const& nalts,
             uvec const& tlens,
             ivec const& ntasks,  
             ivec const& xfr,
             ivec const& xto,  
             ivec const& lfr,  
             ivec const& lto,
             cube const& thetaDraw,
             cube const& tauDraw, 
             int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);

  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs = PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-2));
        double beta_p = exp(theta(p-1));
        arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
        arma::vec pr = exp(ab)/(1+sum(exp(ab)));
        
        arma::vec taui  = tauDraw.slice(ir).col(n);
        pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
        
        //multinomial draw
        int pick_draw = rmuno2(pr);
        
        //if not outside good, choose inside good
        if(pick_draw!=nalt){
          demcontainer(pick_draw,ir)=1;
        }
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


//' @export
//[[Rcpp::export]]
List ddsrprdem(vec const& PP,
               mat const& AA,
               mat const& AAf,
               uvec const& nalts,
               uvec const& tlens,
               ivec const& ntasks,  
               ivec const& xfr,
               ivec const& xto,  
               ivec const& lfr,  
               ivec const& lto,
               cube const& thetaDraw,
               cube const& tauDraw, 
               mat const& tau_pr_Draw,
               int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      ivec nalt_space = linspace<ivec>(0, nalt-1); 
      
      arma::vec prcs = PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-2));
        double beta_p = exp(theta(p-1));
        arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
        arma::vec pr = exp(ab)/(1+sum(exp(ab)));
        
        arma::vec taui  = tauDraw.slice(ir).col(n);
        pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
        pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
        
        //multinomial draw
        int pick_draw = rmuno2(pr);
        
        //if not outside good, choose inside good
        if(pick_draw!=nalt){
          demcontainer(pick_draw,ir)=1;
        }
        
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      } 
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


///////////////////////////////////////////////
// DD Demand-Prob
///////////////////////////////////////////////


//' @export
//[[Rcpp::export]]
List ddprob(vec const& PP,
            mat const& AA,
            uvec const& nalts,
            uvec const& tlens,
            ivec const& ntasks,  
            ivec const& xfr,
            ivec const& xto,  
            ivec const& lfr,  
            ivec const& lto,
            cube const& thetaDraw,
            int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);

  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs = PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-2));
        double beta_p = exp(theta(p-1));
        
        arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
        arma::vec pr = exp(ab)/(1+sum(exp(ab)));

        demcontainer.col(ir)=pr;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}




//' @export
//[[Rcpp::export]]
List ddsrprob(vec const& PP,
              mat const& AA,
              mat const& AAf,
              uvec const& nalts,
              uvec const& tlens,
              ivec const& ntasks,  
              ivec const& xfr,
              ivec const& xto,  
              ivec const& lfr,  
              ivec const& lto,
              cube const& thetaDraw,
              cube const& tauDraw, 
              int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs = PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-2));
        double beta_p = exp(theta(p-1));
        arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
        arma::vec pr = exp(ab)/(1+sum(exp(ab)));
        
        arma::vec taui  = tauDraw.slice(ir).col(n);
        pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;

        demcontainer.col(ir)=pr;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }
      
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


//' @export
//[[Rcpp::export]]
List ddsrprprob(vec const& PP,
             mat const& AA,
             mat const& AAf,
             uvec const& nalts,
             uvec const& tlens,
             ivec const& ntasks,  
             ivec const& xfr,
             ivec const& xto,  
             ivec const& lfr,  
             ivec const& lto,
             cube const& thetaDraw,
             cube const& tauDraw, 
             mat const& tau_pr_Draw,
             int cores=1){
 
 //dimensions
 int R=thetaDraw.n_slices;
 int xdim = PP.n_rows;
 int N = tlens.n_elem;
 int p = thetaDraw.n_rows;
 
 //output init
 Rcpp::List XdL(xdim);
 
 //start timer
 startTimer();
 
 //resp-level
 for(int n=0; n<N; n++){
   infoTimer(n,N);
   
   int ntask = tlens(n);
   int xpick = xfr(n);
   
   //task-level
   for(int tt=0; tt<ntask; tt++){
     Rcpp::checkUserInterrupt();
     
     int nalt = nalts(tt);
     
     //temp storage
     mat demcontainer(nalt,R, fill::zeros);
     
     arma::vec prcs = PP(span(xpick,xpick+nalt-1));
     
     //draw-level
     omp_set_num_threads(cores);
     #pragma omp parallel for schedule(static)
     for (int ir = 0; ir <R; ++ir){
       
       //paras
       arma::vec theta = thetaDraw.slice(ir).col(n);
       arma::vec beta  = theta(arma::span(0,p-2));
       double beta_p = exp(theta(p-1));
       arma::vec ab = AA(span(xpick,xpick+nalt-1),span::all)*beta - prcs*beta_p;
       arma::vec pr = exp(ab)/(1+sum(exp(ab)));
       
       arma::vec taui  = tauDraw.slice(ir).col(n);
       pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
       pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
       
       demcontainer.col(ir)=pr;
     }
     
     for(int k=0; k<nalt; k++){
       XdL(k+xpick)=demcontainer.row(k);
     }
     
     xpick+=nalt;
   }
   
 }
 return(XdL);
}

               
               
               

/////////////////////////////////////// 
// VD Compensatory - EV error
///////////////////////////////////////

//' @export
// [[Rcpp::export]]
double vdl2(arma::vec const& theta, 
            arma::uvec const& nalts,
            arma::vec const& sumpxs, 
            arma::vec const& X, 
            arma::vec const& P, 
            arma::mat const& A, 
            int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-4));
  double bud = exp(theta(p-1));
  double gamma = exp(theta(p-2));
  double sigma = exp(theta(p-3));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        double gt = -ab+log(p)+lngx1-log(osg);
        
        ll+= -exp(-gt/sigma)-gt/sigma-log(sigma);        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        //corner
        double gt = -ab+log(p)-log(osg);
        ll+=-exp(-gt/sigma);        
        //screening
        //(afull*tau)>0
        
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}

//' @export
// [[Rcpp::export]]
vec vd2LL(mat const&Theta,
          vec const& XX, 
          vec const& PP,
          mat const& AA,
          uvec const& nalts,
          vec const& sumpxs,  
          ivec const& ntasks,  
          ivec const& xfr,  
          ivec const& xto,  
          ivec const& lfr,  
          ivec const& lto,
          int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= vdl2(Theta.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p);
  }
  
  return(ll_olds);
}

//' @export
// [[Rcpp::export]]
mat vd2LLs(cube const&THETAS,
           vec const& XX, 
           vec const& PP,
           mat const& AA,
           uvec const& nalts,
           vec const& sumpxs,  
           ivec const& ntasks,  
           ivec const& xfr,  
           ivec const& xto,  
           ivec const& lfr,  
           ivec const& lto,
           int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);
  
  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      vd2LL(THETAS.slice(r),
            XX, 
            PP,
            AA,
            nalts,
            sumpxs,  
            ntasks,  
            xfr,  
            xto,  
            lfr,  
            lto,
            p,
            N, cores);
  }
  
  return(ll_olds);
}



//i-level draws RWMH
void draw_vd2_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                   arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                   arma::mat& theta_temp,    // container of current betas for all i
                   vec const& XX,             // data
                   vec const& PP,
                   mat const& AA,
                   uvec const& nalts,
                   vec const& sumpxs,  
                   ivec const& ntasks,  
                   ivec const& xfr,ivec const& xto,  
                   ivec const& lfr,ivec const& lto,
                   vec const& maxspents,
                   int p, int N, 
                   arma::vec const& mu,  // upper level mean
                   arma::mat const& L,   // upper level chol(sigma)
                   arma::vec& stay,      // rejection tracker, used for tuning
                   arma::vec& tunes,     // i-level tuning parameters
                   int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdl2(theta_cand,
                   nalts(span(lfr(n),lto(n))),
                   sumpxs(span(lfr(n),lto(n))), 
                   XX(span(xfr(n),xto(n))), 
                   PP(span(xfr(n),xto(n))), 
                   AA(span(xfr(n),xto(n)),span::all), 
                   ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}


//' @export
// [[Rcpp::export]]
List loop_vd2_RWMH( vec const& XX, 
                    vec const& PP,
                    mat const& AA,
                    uvec const& nalts,
                    vec const& sumpxs,  
                    ivec const& ntasks,  
                    ivec const& xfr,  
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    int p, int N,
                    int R, int keep, // MCMC parameters draws and keep interval
                    mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                    int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                    int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= vdl2(theta_temp.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    draw_vd2_RWMH(ll_olds,            // ll for current betas
                  lp_olds,
                  theta_temp,          // container of current betas for all i
                  XX,
                  PP,
                  AA,
                  nalts,
                  sumpxs,
                  ntasks,
                  xfr,xto,lfr,lto,
                  maxspents,
                  p,N,
                  vectorise(MU),      // Mean Prior
                  Lprior,             // VarCov Prior (chol)
                  stay,               // tracking rejections
                  tunes,              // i-level tuning parameters
                  cores);       
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}






/////////////////////////////////////// 
// VD Compensatory - Normal Error
///////////////////////////////////////

//' @export
// [[Rcpp::export]]
double vdln(arma::vec const& theta, 
            arma::uvec const& nalts,
            arma::vec const& sumpxs, 
            arma::vec const& X, 
            arma::vec const& P, 
            arma::mat const& A, 
            int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-4));
  double bud = exp(theta(p-1));
  double gamma = exp(theta(p-2));
  double sigma = exp(theta(p-3));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg))/sigma)-log(sigma));
        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        //corner
        ll+=log(normcdf((-ab+log(p)-log(osg))/sigma));
        
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}

//' @export
// [[Rcpp::export]]
vec vdnLL(mat const&Theta,
          vec const& XX, 
          vec const& PP,
          mat const& AA,
          uvec const& nalts,
          vec const& sumpxs,  
          ivec const& ntasks,  
          ivec const& xfr,  
          ivec const& xto,  
          ivec const& lfr,  
          ivec const& lto,
          int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= vdln(Theta.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p);
  }
  
  return(ll_olds);
}

//' @export
// [[Rcpp::export]]
mat vdnLLs(cube const&THETAS,
           vec const& XX, 
           vec const& PP,
           mat const& AA,
           uvec const& nalts,
           vec const& sumpxs,  
           ivec const& ntasks,  
           ivec const& xfr,  
           ivec const& xto,  
           ivec const& lfr,  
           ivec const& lto,
           int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);
  
  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      vdnLL(THETAS.slice(r),
            XX, 
            PP,
            AA,
            nalts,
            sumpxs,  
            ntasks,  
            xfr,  
            xto,  
            lfr,  
            lto,
            p,
            N, cores);
  }
  
  return(ll_olds);
}



//i-level draws RWMH
void draw_vdn_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                   arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                   arma::mat& theta_temp,    // container of current betas for all i
                   vec const& XX,             // data
                   vec const& PP,
                   mat const& AA,
                   uvec const& nalts,
                   vec const& sumpxs,  
                   ivec const& ntasks,  
                   ivec const& xfr,ivec const& xto,  
                   ivec const& lfr,ivec const& lto,
                   vec const& maxspents,
                   int p, int N, 
                   arma::vec const& mu,  // upper level mean
                   arma::mat const& L,   // upper level chol(sigma)
                   arma::vec& stay,      // rejection tracker, used for tuning
                   arma::vec& tunes,     // i-level tuning parameters
                   int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdln(theta_cand,
                   nalts(span(lfr(n),lto(n))),
                   sumpxs(span(lfr(n),lto(n))), 
                   XX(span(xfr(n),xto(n))), 
                   PP(span(xfr(n),xto(n))), 
                   AA(span(xfr(n),xto(n)),span::all), 
                   ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}


//' @export
// [[Rcpp::export]]
List loop_vdn_RWMH( vec const& XX, 
                    vec const& PP,
                    mat const& AA,
                    uvec const& nalts,
                    vec const& sumpxs,  
                    ivec const& ntasks,  
                    ivec const& xfr,  
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    int p, int N,
                    int R, int keep, // MCMC parameters draws and keep interval
                    mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                    int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                    int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= vdln(theta_temp.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    draw_vdn_RWMH(ll_olds,            // ll for current betas
                  lp_olds,
                  theta_temp,          // container of current betas for all i
                  XX,
                  PP,
                  AA,
                  nalts,
                  sumpxs,
                  ntasks,
                  xfr,xto,lfr,lto,
                  maxspents,
                  p,N,
                  vectorise(MU),      // Mean Prior
                  Lprior,             // VarCov Prior (chol)
                  stay,               // tracking rejections
                  tunes,              // i-level tuning parameters
                  cores);       
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}








/////////////////////////////////////// 
// VD Screening model
///////////////////////////////////////

///tauconst [nattrf,N] - 0 if bought, 1 if not and thus screenable


//' @export
// [[Rcpp::export]]
double vdlsr2( arma::vec const& theta, 
               arma::ivec const& taui,
               arma::uvec const& nalts,
               arma::vec const& sumpxs, 
               arma::vec const& X, 
               arma::vec const& P, 
               arma::mat const& A, 
               arma::mat const& Afull, 
               int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-4));
  double bud = exp(theta(p-1));
  double gamma = exp(theta(p-2));
  double sigma = exp(theta(p-3));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg))/sigma)-log(sigma));
        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        
        //screening check
        if(as_scalar(Afull.row(xpicker)*taui)>(0.01)){
          
        }else{
          //corner
          ll+=log(normcdf((-ab+log(p)-log(osg))/sigma));
        }
        
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}

//' @export
// [[Rcpp::export]]
vec vdsr2LL( mat const&Theta,
             imat const&tauis,
             vec const& XX, 
             vec const& PP,
             mat const& AA,
             mat const& AAf,
             uvec const& nalts,
             vec const& sumpxs,  
             ivec const& ntasks,  
             ivec const& xfr,  
             ivec const& xto,  
             ivec const& lfr,  
             ivec const& lto,
             int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
  #pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= vdlsr2(
            Theta.col(n),
            tauis.col(n),
              nalts(span(lfr(n),lto(n))),
              sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            AAf(span(xfr(n),xto(n)),span::all), 
              ntasks(n), p);
  }
  
  return(ll_olds);
}


//' @export
// [[Rcpp::export]]
mat vdsr2LLs(cube const&THETAS,
             icube const&TAUIS,
            vec const& XX, 
            vec const& PP,
            mat const& AA,
            mat const& AAf,
            uvec const& nalts,
            vec const& sumpxs,  
            ivec const& ntasks,  
            ivec const& xfr,  
            ivec const& xto,  
            ivec const& lfr,  
            ivec const& lto,
            int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);

  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      vdsr2LL(THETAS.slice(r),
              TAUIS.slice(r),
                XX, 
                PP,
                AA,
                AAf,
                nalts,
                sumpxs,  
                ntasks,  
                xfr,  
                xto,  
                lfr,  
                lto,
                p,
                N, cores);
  }
  
  return(ll_olds);
}







void draw_vdsr2_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                     arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                     arma::mat& theta_temp,    // container of current betas for all i
                     arma::imat const& tauis, 
                     arma::imat const& tauconsts,
                     vec const& XX,             // data
                     vec const& PP,
                     mat const& AA,
                     mat const& AAf,
                     uvec const& nalts,
                     vec const& sumpxs,  
                     ivec const& ntasks,  
                     ivec const& xfr,ivec const& xto,  
                     ivec const& lfr,ivec const& lto,
                     vec const& maxspents,
                     int p, int N, 
                     arma::vec const& mu,  // upper level mean
                     arma::mat const& L,   // upper level chol(sigma)
                     arma::vec& stay,      // rejection tracker, used for tuning
                     arma::vec& tunes,     // i-level tuning parameters
                     int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdlsr2(theta_cand,
                     tauis.col(n),
                     nalts(span(lfr(n),lto(n))),
                     sumpxs(span(lfr(n),lto(n))), 
                     XX(span(xfr(n),xto(n))), 
                     PP(span(xfr(n),xto(n))), 
                     AA(span(xfr(n),xto(n)),span::all), 
                     AAf(span(xfr(n),xto(n)),span::all), 
                     ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}


void draw_tau(arma::vec& ll_olds,
              imat& tauis,
              arma::mat const& theta_temp,
              imat const& tauconst,
              mat const& delta,
              vec const& XX, 
              vec const& PP,
              mat const& AA,
              mat const& AAf,
              uvec const& nalts,
              vec const& sumpxs,  
              ivec const& ntasks,  
              ivec const& xfr,  
              ivec const& xto,  
              ivec const& lfr,  
              ivec const& lto,
              int p, int N, int cores){
  
  int K=tauconst.n_rows;
  
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    double ll0;
    double ll1;
    
    for(int kk=0; kk<K; kk++){
      
      if(tauconst(kk,n)==1){
        
        if(tauis(kk,n)==1){
          
          ll1=ll_olds(n);
          
          ivec tauk0=tauis.col(n);
          tauk0(kk)=0;
          ll0 = vdlsr2(theta_temp.col(n),
                       tauk0,
                       nalts(span(lfr(n),lto(n))),
                       sumpxs(span(lfr(n),lto(n))), 
                       XX(span(xfr(n),xto(n))), 
                       PP(span(xfr(n),xto(n))), 
                       AA(span(xfr(n),xto(n)),span::all), 
                       AAf(span(xfr(n),xto(n)),span::all), 
                       ntasks(n), p );
        }else{
          ll0=ll_olds(n);
          
          ivec tauk1=tauis.col(n);
          tauk1(kk)=1;
          ll1 = vdlsr2(theta_temp.col(n),
                       tauk1,
                       nalts(span(lfr(n),lto(n))),
                       sumpxs(span(lfr(n),lto(n))), 
                       XX(span(xfr(n),xto(n))), 
                       PP(span(xfr(n),xto(n))), 
                       AA(span(xfr(n),xto(n)),span::all), 
                       AAf(span(xfr(n),xto(n)),span::all), 
                       ntasks(n), p );  
        }
        
        double probsc = (exp(ll1) * (delta(kk))) / 
          (exp(ll1) * (delta(kk)) + exp(ll0)*(1-delta(kk)  )); 
        
        tauis(kk,n)= Rf_rbinom( 1, probsc );
        
        if(tauis(kk,n)==1){
          ll_olds(n)=ll1;
        }else{
          ll_olds(n)=ll0;
        }
        
      }
    }//k=loop
  }//nloop
}

//' @export
// [[Rcpp::export]]
List loop_vdrs2_RWMH( vec const& XX, 
                      vec const& PP,
                      mat const& AA,
                      mat const& AAf,
                      imat const& tauconst,
                      uvec const& nalts,
                      vec const& sumpxs,  
                      ivec const& ntasks,  
                      ivec const& xfr,  
                      ivec const& xto,  
                      ivec const& lfr,  
                      ivec const& lto,
                      int p, int N,
                      int R, int keep, // MCMC parameters draws and keep interval
                      mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                      int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                      int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  int K=tauconst.n_rows;
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  imat tauis   = tauconst;
  vec delta(K); delta.fill(0.5);
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= vdlsr2(theta_temp.col(n),
            tauis.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            AAf(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  icube tauDraw(K,N,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  mat  deltaDraw(Rk,K);  
  
  mat  loglikeM = zeros<mat>(N,Rk);
  
  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  omp_set_num_threads(cores);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    
    // n loop  *********
    //theta
    draw_vdsr2_RWMH(ll_olds,            // ll for current betas
                    lp_olds,
                    theta_temp,          // container of current betas for all i
                    tauis,
                    tauconst,
                    XX,
                    PP,
                    AA,
                    AAf,
                    nalts,
                    sumpxs,
                    ntasks,
                    xfr,xto,lfr,lto,
                    maxspents,
                    p,N,
                    vectorise(MU),      // Mean Prior
                    Lprior,             // VarCov Prior (chol)
                    stay,               // tracking rejections
                    tunes,              // i-level tuning parameters
                    cores);       
    
    if(ir>1000){
      //tau
      draw_tau(ll_olds,
               tauis,
               theta_temp,
               tauconst, 
               delta,
               XX, 
               PP,
               AA,
               AAf,
               nalts,
               sumpxs,  
               ntasks,  
               xfr, xto, lfr,lto,
               p, N, 
               cores);
      
      // delta_tau
      
      drawdelta(delta, tauis, K, N, cores);
      
      //update LL
      // #pragma omp parallel for schedule(static)
      // for(int n=0; n<N; n++){
      //   
      //   ll_olds(n)= 
      //     vdlsr2(theta_temp.col(n),
      //           tauis.col(n),
      //           nalts(span(lfr(n),lto(n))),
      //           sumpxs(span(lfr(n),lto(n))), 
      //           XX(span(xfr(n),xto(n))), 
      //           PP(span(xfr(n),xto(n))), 
      //           AA(span(xfr(n),xto(n)),span::all), 
      //           AAf(span(xfr(n),xto(n)),span::all), 
      //           ntasks(n), p );
      // }
      
    }
    // for(int kk=0; kk<K; kk++){
    //   double tsum=sum(tauis.row(kk));
    //   delta(kk)=rbeta(1,tsum, N-tsum+1)(0);
    // }
    // vec tsums = sum(tauis,1)
    
    
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      
      
      tauDraw.slice(mkeep)        = tauis;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      deltaDraw.row(mkeep)        = trans(delta);
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
      
      loglikeM.col(mkeep) =ll_olds;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("tauDraw")        = tauDraw,
    Named("deltaDraw")      = deltaDraw,
    Named("loglike")        = loglike,
    Named("loglikeM")       = loglikeM,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}






///////////////////////////////////////////////
// Volumetric - Conjunctive Screening, including price tag
///////////////////////////////////////////////

//' @export
// [[Rcpp::export]]
double vdlsrpr( arma::vec const& theta, 
                arma::ivec const& taui,
                double tau_pr,
                arma::uvec const& nalts,
                arma::vec const& sumpxs, 
                arma::vec const& X, 
                arma::vec const& P, 
                arma::mat const& A, 
                arma::mat const& Afull, 
                int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-4));
  double bud   = exp(theta(p-1));
  double gamma = exp(theta(p-2));
  double sigma = exp(theta(p-3));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg))/sigma)-log(sigma));
        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        //screening check; acceptable price and attributes
        if(as_scalar(Afull.row(xpicker)*taui)>(0.01)){
          
        }else{
          //check price tag screen
          if(p<=exp(tau_pr)){
            //corner
            ll+=log(normcdf((-ab+log(p)-log(osg))/sigma));
          }
        }
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}




void draw_taui_pr(arma::vec& ll_olds,
                  imat& tauis,
                  arma::mat const& theta_temp,
                  vec const& tau_prs,
                  imat const& tauconst,
                  mat const& delta,
                  vec const& XX, 
                  vec const& PP,
                  mat const& AA,
                  mat const& AAf,
                  uvec const& nalts,
                  vec const& sumpxs,  
                  ivec const& ntasks,  
                  ivec const& xfr,  
                  ivec const& xto,  
                  ivec const& lfr,  
                  ivec const& lto,
                  int p, int N, int cores){
  
  int K=tauconst.n_rows;
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    double ll0;
    double ll1;
    
    for(int kk=0; kk<K; kk++){
      
      
      if(tauconst(kk,n)==1){
        
        if(tauis(kk,n)==1){
          
          
          ll1=ll_olds(n);
          
          ivec tauk0=tauis.col(n);
          tauk0(kk)=0;
          ll0 = vdlsrpr(theta_temp.col(n),
                        tauk0,
                        tau_prs(n),
                        nalts(span(lfr(n),lto(n))),
                        sumpxs(span(lfr(n),lto(n))), 
                        XX(span(xfr(n),xto(n))), 
                        PP(span(xfr(n),xto(n))), 
                        AA(span(xfr(n),xto(n)),span::all), 
                        AAf(span(xfr(n),xto(n)),span::all), 
                        ntasks(n), p );
        }else{
          ll0=ll_olds(n);
          
          ivec tauk1=tauis.col(n);
          tauk1(kk)=1;
          ll1 = vdlsrpr(theta_temp.col(n),
                        tauk1,
                        tau_prs(n),
                        nalts(span(lfr(n),lto(n))),
                        sumpxs(span(lfr(n),lto(n))), 
                        XX(span(xfr(n),xto(n))), 
                        PP(span(xfr(n),xto(n))), 
                        AA(span(xfr(n),xto(n)),span::all), 
                        AAf(span(xfr(n),xto(n)),span::all), 
                        ntasks(n), p );
        }
        
        double probsc = (exp(ll1) * (delta(kk))) / 
          (exp(ll1) * (delta(kk)) + exp(ll0)*(1-delta(kk)  )); 
        
        tauis(kk,n)= Rf_rbinom( 1, probsc );
        
        if(tauis(kk,n)==1){
          ll_olds(n)=ll1;
        }else{
          ll_olds(n)=ll0;
        }
        
      }      
    }//k=loop
  }//nloop
}



void draw_taupr(arma::vec& ll_olds,    // vector of current log-likelihoods
                arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                arma::imat const& tauis,
                arma::mat const& theta_temp,
                vec& tau_prs,
                vec const& maxpaids,
                double const& pr_mean,
                double const& pr_sd,
                vec& stay,
                vec const& pricetune,
                vec const& XX, 
                vec const& PP,
                mat const& AA,
                mat const& AAf,
                uvec const& nalts,
                vec const& sumpxs,  
                ivec const& ntasks,  
                ivec const& xfr,  
                ivec const& xto,  
                ivec const& lfr,  
                ivec const& lto,
                int p, int N, int cores){
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //candidate
    double tau_pr_cand = tau_prs(n) + pricetune(n)*randn(1)[0];
    
    
    //check
    if(tau_pr_cand>log(maxpaids(n))){
      
      
      double llnew = vdlsrpr(theta_temp.col(n),
                             tauis.col(n),
                             tau_pr_cand,
                             nalts(span(lfr(n),lto(n))),
                             sumpxs(span(lfr(n),lto(n))), 
                             XX(span(xfr(n),xto(n))), 
                             PP(span(xfr(n),xto(n))), 
                             AA(span(xfr(n),xto(n)),span::all), 
                             AAf(span(xfr(n),xto(n)),span::all), 
                             ntasks(n), p );
      
      
      //A-R
      double ldiff = 
        llnew      + log_normpdf(tau_pr_cand, pr_mean, pr_sd ) - 
        ll_olds(n) - log_normpdf(tau_prs(n) , pr_mean, pr_sd );
      
      
      if(ldiff > log(randu(1)[0])){
        tau_prs(n)  = tau_pr_cand;
        ll_olds(n)  = llnew;
        // lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
    }else{
      stay(n)+=1;
    }
    
    
  }//nloop
}





void draw_vdspr_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                     arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                     arma::mat& theta_temp,    // container of current betas for all i
                     arma::imat& tauis, 
                     vec const& tau_prs,
                     arma::imat const& tauconsts,
                     vec const& XX,             // data
                     vec const& PP,
                     mat const& AA,
                     mat const& AAf,
                     uvec const& nalts,
                     vec const& sumpxs,  
                     ivec const& ntasks,  
                     ivec const& xfr,ivec const& xto,  
                     ivec const& lfr,ivec const& lto,
                     vec const& maxspents,
                     int p, int N, 
                     arma::vec const& mu,  // upper level mean
                     arma::mat const& L,   // upper level chol(sigma)
                     arma::vec& stay,      // rejection tracker, used for tuning
                     arma::vec& tunes,     // i-level tuning parameters
                     int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdlsrpr(theta_cand,
                      tauis.col(n),
                      tau_prs(n),
                      nalts(span(lfr(n),lto(n))),
                      sumpxs(span(lfr(n),lto(n))), 
                      XX(span(xfr(n),xto(n))), 
                      PP(span(xfr(n),xto(n))), 
                      AA(span(xfr(n),xto(n)),span::all), 
                      AAf(span(xfr(n),xto(n)),span::all), 
                      ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}







//' @export
// [[Rcpp::export]]
List loop_vdrspr_RWMH( vec const& XX, 
                       vec const& PP,
                       mat const& AA,
                       mat const& AAf,
                       imat const& tauconst,
                       uvec const& nalts,
                       vec const& sumpxs,  
                       ivec const& ntasks,  
                       ivec const& xfr,  
                       ivec const& xto,  
                       ivec const& lfr,  
                       ivec const& lto,
                       int p, int N,
                       int R, int keep, // MCMC parameters draws and keep interval
                       mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                       int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                       int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  vec maxpaids(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
    maxpaids(n) = max(sign(XX(span(xfr(n),xto(n))))%PP(span(xfr(n),xto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  int K=tauconst.n_rows;
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  imat tauis   = tauconst;
  vec delta(K); delta.fill(0.5);
  
  vec tau_prs = maxpaids*1.1;
  double pr_mean = mean(tau_prs);
  double pr_sd = 1;
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  arma::vec stay_prscr(N);
  stay_prscr.fill(0);
  arma::vec pricetunes(N);
  pricetunes.fill(.1);
  
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= vdlsrpr(theta_temp.col(n),
            tauis.col(n),
            tau_prs(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            AAf(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  icube tauDraw(K,N,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  mat  deltaDraw(Rk,K);  
  mat  pricescreenPriorDraw(Rk,2);
  mat  tau_pr_draw(N,Rk);
  
  mat  loglikeM = zeros<mat>(N,Rk);
  
  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  omp_set_num_threads(cores);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    
    // n loop  *********
    //theta
    draw_vdspr_RWMH(ll_olds,            // ll for current betas
                    lp_olds,
                    theta_temp,          // container of current betas for all i
                    tauis,
                    tau_prs,
                    tauconst,
                    XX,
                    PP,
                    AA,
                    AAf,
                    nalts,
                    sumpxs,
                    ntasks,
                    xfr,xto,lfr,lto,
                    maxspents,
                    p,N,
                    vectorise(MU),      // Mean Prior
                    Lprior,             // VarCov Prior (chol)
                    stay,               // tracking rejections
                    tunes,              // i-level tuning parameters
                    cores);       
    
    if(ir>1000){
      //tau
      draw_taui_pr(ll_olds,
                   tauis,
                   theta_temp,
                   tau_prs,
                   tauconst,
                   delta,
                   XX, 
                   PP,
                   AA,
                   AAf,
                   nalts,
                   sumpxs,  
                   ntasks,  
                   xfr, xto, lfr,lto,
                   p, N, 
                   cores);
      
      // delta_tau
      drawdelta(delta, tauis, K, N, cores);
      
      
      //update LL (integreate into draw_tau)
      // #pragma omp parallel for schedule(static)
      //       for(int n=0; n<N; n++){
      //         ll_olds(n)= 
      //           vdlsrpr(theta_temp.col(n),
      //                   tauis.col(n),
      //                   tau_prs(n),
      //                   nalts(span(lfr(n),lto(n))),
      //                   sumpxs(span(lfr(n),lto(n))), 
      //                   XX(span(xfr(n),xto(n))), 
      //                   PP(span(xfr(n),xto(n))), 
      //                   AA(span(xfr(n),xto(n)),span::all), 
      //                   AAf(span(xfr(n),xto(n)),span::all), 
      //                   ntasks(n), p );
      //       }
      
      //update price screening
      draw_taupr(ll_olds,    // vector of current log-likelihoods
                 lp_olds,       // vectors of lp's, just for tracking
                 tauis,
                 theta_temp,
                 tau_prs,
                 maxpaids,
                 pr_mean,
                 pr_sd,
                 stay_prscr,
                 pricetunes,
                 XX,
                 PP,
                 AA,
                 AAf,
                 nalts,
                 sumpxs,
                 ntasks,
                 xfr,
                 xto,
                 lfr,
                 lto,
                 p, N, cores);
      
      //price screening upper level mu_0, nu, alph, bet
      ULnormnorm(pr_mean, pr_sd,
                 tau_prs,
                 0.0, 0.01, 3.0, 3.0);
      
    }
    // for(int kk=0; kk<K; kk++){
    //   double tsum=sum(tauis.row(kk));
    //   delta(kk)=rbeta(1,tsum, N-tsum+1)(0);
    // }
    // vec tsums = sum(tauis,1)
    
    
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
        mh_tuner(pricetunes,stay_prscr/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      tauDraw.slice(mkeep)        = tauis;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      deltaDraw.row(mkeep)        = trans(delta);
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
      tau_pr_draw.col(mkeep)        = tau_prs;
      pricescreenPriorDraw(mkeep,0) = pr_mean;
      pricescreenPriorDraw(mkeep,1) = pr_sd;
      loglikeM.col(mkeep) =ll_olds;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("tauDraw")        = tauDraw,
    Named("deltaDraw")      = deltaDraw,
    Named("tau_pr_draw")    = tau_pr_draw,
    Named("prscreenMuSigDraw")  = pricescreenPriorDraw,
    Named("loglike")        = loglike,
    Named("loglikeM")       = loglikeM,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}



//' @export
// [[Rcpp::export]]
vec vdsrprLL(mat const&Theta,
             imat const&tauis,
             vec const&tau_prs,
             vec const& XX, 
             vec const& PP,
             mat const& AA,
             mat const& AAf,
             uvec const& nalts,
             vec const& sumpxs,  
             ivec const& ntasks,  
             ivec const& xfr,  
             ivec const& xto,  
             ivec const& lfr,  
             ivec const& lto,
             int p, int N, int cores=1){
  
  omp_set_num_threads(cores);
  
  vec ll_olds(N);
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    ll_olds(n)= vdlsrpr(
      Theta.col(n),
      tauis.col(n),
      tau_prs(n),
      nalts(span(lfr(n),lto(n))),
      sumpxs(span(lfr(n),lto(n))), 
      XX(span(xfr(n),xto(n))), 
      PP(span(xfr(n),xto(n))), 
      AA(span(xfr(n),xto(n)),span::all), 
      AAf(span(xfr(n),xto(n)),span::all), 
      ntasks(n), p);
  }
  
  return(ll_olds);
}

//' @export
// [[Rcpp::export]]
mat vdsrprLLs(cube const&THETAS,
              icube const&TAUIS,
              mat const&TAU_PR,
              vec const& XX, 
              vec const& PP,
              mat const& AA,
              mat const& AAf,
              uvec const& nalts,
              vec const& sumpxs,  
              ivec const& ntasks,  
              ivec const& xfr,  
              ivec const& xto,  
              ivec const& lfr,  
              ivec const& lto,
              int p, int N, int cores=1){
  
  int R = THETAS.n_slices;
  mat ll_olds(N,R+1);
  
  for(int r=0; r<R; r++){
    Rcpp::checkUserInterrupt();
    ll_olds.col(r)= 
      vdsrprLL(THETAS.slice(r),
               TAUIS.slice(r),
               TAU_PR.col(r),
               XX, 
               PP,
               AA,
               AAf,
               nalts,
               sumpxs,  
               ntasks,  
               xfr, xto, lfr, lto,
               p, N, cores);
  }
  
  return(ll_olds);
}



/////////////////////////////////////// 
// VD Set-Size Compensatory - EV error
///////////////////////////////////////


//' @export
// [[Rcpp::export]]
double vdl_ss(arma::vec const& theta, 
              arma::uvec const& nalts,
              arma::vec const& sumpxs, 
              arma::vec const& X, 
              arma::vec const& P, 
              arma::mat const& A, 
              int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-5));
  double bud      = exp(theta(p-1));
  double gamma    = exp(theta(p-2));
  double sigma    = exp(theta(p-3));
  double xi       = exp(theta(p-4));
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        //ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg)+log(xi*nalt+1))/sigma)-log(sigma));
        
        double gt = -ab+log(p)+lngx1-log(osg)+log(xi*nalt+1);
        // ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg)+log(xi*nalt+1))/sigma)-log(sigma));
        
        ll+= -exp(-gt/sigma)-gt/sigma-log(sigma);
        
        //sum(-gt.elem( find(xav>0) )/sigma-log(sigma));  
        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        //corner
        // ll+=log(normcdf((-ab+log(p)+log(gamma*x+1)-log(osg))/sigma));
        
        double gt = -ab+log(p)-log(osg)+log(xi*nalt+1);
        ll+=-exp(-gt/sigma);
        
        
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}



//i-level draws RWMH
void draw_vdl_ss_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                      arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                      arma::mat& theta_temp,    // container of current betas for all i
                      vec const& XX,             // data
                      vec const& PP,
                      mat const& AA,
                      uvec const& nalts,
                      vec const& sumpxs,  
                      ivec const& ntasks,  
                      ivec const& xfr,ivec const& xto,  
                      ivec const& lfr,ivec const& lto,
                      vec const& maxspents,
                      int p, int N, 
                      arma::vec const& mu,  // upper level mean
                      arma::mat const& L,   // upper level chol(sigma)
                      arma::vec& stay,      // rejection tracker, used for tuning
                      arma::vec& tunes,     // i-level tuning parameters
                      int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdl_ss(theta_cand,
                     nalts(span(lfr(n),lto(n))),
                     sumpxs(span(lfr(n),lto(n))), 
                     XX(span(xfr(n),xto(n))), 
                     PP(span(xfr(n),xto(n))), 
                     AA(span(xfr(n),xto(n)),span::all), 
                     ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}


//' @export
// [[Rcpp::export]]
List loop_vd_ss_RWMH( vec const& XX, 
                      vec const& PP,
                      mat const& AA,
                      uvec const& nalts,
                      vec const& sumpxs,  
                      ivec const& ntasks,  
                      ivec const& xfr,  
                      ivec const& xto,  
                      ivec const& lfr,  
                      ivec const& lto,
                      int p, int N,
                      int R, int keep, // MCMC parameters draws and keep interval
                      mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                      int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                      int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  theta_temp.row(p-4) += -6;
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= vdl_ss(theta_temp.col(n),
            nalts(span(lfr(n),lto(n))),
            sumpxs(span(lfr(n),lto(n))), 
            XX(span(xfr(n),xto(n))), 
            PP(span(xfr(n),xto(n))), 
            AA(span(xfr(n),xto(n)),span::all), 
            ntasks(n), p );
  }
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    draw_vdl_ss_RWMH(ll_olds,            // ll for current betas
                     lp_olds,
                     theta_temp,          // container of current betas for all i
                     XX,
                     PP,
                     AA,
                     nalts,
                     sumpxs,
                     ntasks,
                     xfr,xto,lfr,lto,
                     maxspents,
                     p,N,
                     vectorise(MU),      // Mean Prior
                     Lprior,             // VarCov Prior (chol)
                     stay,               // tracking rejections
                     tunes,              // i-level tuning parameters
                     cores);       
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}





/////////////////////////////////////// 
// VD Set-Size (with quadratic term) Compensatory - EV error
///////////////////////////////////////


//' @export
// [[Rcpp::export]]
double vdl_ssQ(arma::vec const& theta, 
               arma::uvec const& nalts,
               arma::vec const& sumpxs, 
               arma::vec const& X, 
               arma::vec const& P, 
               arma::mat const& A, 
               int ntask, int p ){
  
  //para
  arma::vec beta = theta(arma::span(0,p-6));
  double bud      = exp(theta(p-1));
  double gamma    = exp(theta(p-2));
  double sigma    = exp(theta(p-3));
  double xi       = exp(theta(p-4));
  double tau      = exp(theta(p-5));
  
  
  //init ll
  double ll=0;
  int xpicker = 0;
  
  //task level
  for(int tt=0; tt<ntask; tt++){
    int nalt = nalts(tt);
    
    double osg = bud-sumpxs(tt);
    double jactemp=0;
    
    
    //alternative level
    for(int kk=0; kk<nalt; kk++){
      
      double x = X(xpicker);
      double p = P(xpicker);
      double ab = as_scalar(A.row(xpicker)*beta);
      
      if(x>0){
        //inside
        double lngx1 = log(gamma*x+1);
        //ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg)+log(xi*nalt+1))/sigma)-log(sigma));
        
        double gt = -ab+log(p)+lngx1-log(osg)+log(xi*nalt+1+tau*(nalt^2));
        
        // double gt = -ab+log(p)+lngx1-log(osg)+
        //   log(xi*pow(nalt,tau)+1);
        
        
        // ll+=(log_normpdf((-ab+log(p)+lngx1-log(osg)+log(xi*nalt+1))/sigma)-log(sigma));
        
        ll+= -exp(-gt/sigma)-gt/sigma-log(sigma);
        
        //sum(-gt.elem( find(xav>0) )/sigma-log(sigma));  
        
        //jacobian
        ll+=log(gamma)-lngx1;
        jactemp+=(gamma*x+1)*p/(osg*gamma);
        
      }else{
        //corner
        // ll+=log(normcdf((-ab+log(p)+log(gamma*x+1)-log(osg))/sigma));
        
        double gt = -ab+log(p)-log(osg)+log(xi*nalt+1+tau*(nalt^2) );
        ll+=-exp(-gt/sigma);
        
        //screening
        //(afull*tau)>0
        
      }
      xpicker+=1; //move up index for X,A,P
    }
    //add 2nd part of jacobian at task level
    ll+=log(jactemp+1);
  }
  return(ll);
}


//i-level draws RWMH
void draw_vdl_ssQ_RWMH(arma::vec& ll_olds,    // vector of current log-likelihoods
                       arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                       arma::mat& theta_temp,    // container of current betas for all i
                       vec const& XX,             // data
                       vec const& PP,
                       mat const& AA,
                       uvec const& nalts,
                       vec const& sumpxs,  
                       ivec const& ntasks,  
                       ivec const& xfr,ivec const& xto,  
                       ivec const& lfr,ivec const& lto,
                       vec const& maxspents,
                       int p, int N, 
                       arma::vec const& mu,  // upper level mean
                       arma::mat const& L,   // upper level chol(sigma)
                       arma::vec& stay,      // rejection tracker, used for tuning
                       arma::vec& tunes,     // i-level tuning parameters
                       int cores=1){ 
  
  omp_set_num_threads(cores);
  
#pragma omp parallel for schedule(static)
  for(int n=0; n<N; n++){
    
    //local variables (thread-safe)
    double llnew;
    double lpnew;
    vec theta_cand = theta_temp.col(n);
    
    //lp old draw (with updated mu,L)
    lp_olds(n) = lndMvnc(theta_temp.col(n), mu, L);
    
    //candidate
    theta_cand+= tunes(n)*(trans(L) * arma::randn(p));
    
    if(exp(theta_cand(p-1))>maxspents(n)){
      
      //eval
      llnew = vdl_ssQ(theta_cand,
                      nalts(span(lfr(n),lto(n))),
                      sumpxs(span(lfr(n),lto(n))), 
                      XX(span(xfr(n),xto(n))), 
                      PP(span(xfr(n),xto(n))), 
                      AA(span(xfr(n),xto(n)),span::all), 
                      ntasks(n), p );
      
      lpnew = lndMvnc(theta_cand, mu, L);
      
      //A-R
      double ldiff = llnew + lpnew - ll_olds(n) - lp_olds(n);
      
      if(ldiff > log(randu(1)[0])){
        theta_temp.col(n)= theta_cand;
        ll_olds(n)       = llnew;
        lp_olds(n)       = lpnew;
      }else{
        stay(n)+=1;
      }
      
    }else{
      stay(n)+=1;
    }
    
  }
}


//' @export
// [[Rcpp::export]]
List loop_vd_ssQ_RWMH( vec const& XX, 
                       vec const& PP,
                       mat const& AA,
                       uvec const& nalts,
                       vec const& sumpxs,  
                       ivec const& ntasks,  
                       ivec const& xfr,  
                       ivec const& xto,  
                       ivec const& lfr,  
                       ivec const& lto,
                       int p, int N,
                       int R, int keep, // MCMC parameters draws and keep interval
                       mat const& Bbar, mat const& A, double nu, mat const& V, //Prior
                       int tuneinterval = 30, double steptunestart=.5, int tunelength=10000, int tunestart=500, //algo settings
                       int progressinterval=100, int cores=1){ //report interval
  
  //initialize i-parameters
  arma::mat theta_temp(p,N); // container of current betas for all i
  theta_temp.fill(0);
  
  vec maxspents(N);
  for(int n=0; n<N; n++){
    maxspents(n) = max(sumpxs(span(lfr(n),lto(n))));
  }
  
  theta_temp.row(p-1) = trans(log(maxspents+0.01));
  theta_temp.row(p-4) += -6;
  
  //no covariates (Z) for now
  mat Z(N,1);
  Z.fill(1);
  
  // dimensions ..................
  int Rk=R/keep;
  int mkeep;
  int m=1;//Z.n_cols; - no covariates for now
  
  // start values ..................
  mat MU      = Bbar;
  mat SIGMA   = eye(p,p);  
  mat Lprior  = trimatu(chol(SIGMA));
  
  // tuning ..................
  arma::vec stay(N);
  arma::vec stay_total(N);
  stay.fill(0);
  stay_total.fill(0);
  
  int tunecounter = 1;
  vec tunes = ones<vec>(N)*steptunestart;
  double currentRR=0;
  vec currentRRs(N);
  
  // initial log likelihood ..................
  vec ll_olds(N);
  for(int n=0; n<N; n++){
    
    ll_olds(n)= 
      vdl_ssQ(theta_temp.col(n),
              nalts(span(lfr(n),lto(n))),
              sumpxs(span(lfr(n),lto(n))), 
              XX(span(xfr(n),xto(n))), 
              PP(span(xfr(n),xto(n))), 
              AA(span(xfr(n),xto(n)),span::all), 
              ntasks(n), p );
  }
  // REprintf("Initial LL:");
  // Rcout << sum(ll_olds);
  // REprintf("\n");
  
  
  vec lp_olds(N);
  for(int n=0; n<N; n++){
    lp_olds(n)=lndMvnc(theta_temp.col(n),vectorise(MU),Lprior);
  }
  
  // draw storage ..................
  cube thetaDraw(p,N,Rk);
  cube SIGMADraw(p,p,Rk);
  
  vec  loglike = zeros<vec>(Rk);
  vec  logpost = zeros<vec>(Rk);
  mat  MUDraw(Rk,p*m);  
  arma::vec rrate(Rk);
  mat RRs(N,Rk);
  
  // loop ..................    
  startMcmcTimer();
  
  for(int ir=0; ir<R; ir++){
    Rcpp::checkUserInterrupt();
    
    // upper level *********
    ULreg(trans(theta_temp), 
          Z,                // upper level covariates if used
          Bbar,  A, nu, V,  // prior
          MU,               // outputs (MU, SIGMA, Lprior=chol(SIGMA))
          SIGMA,
          Lprior); 
    
    // n loop  *********
    draw_vdl_ssQ_RWMH(ll_olds,            // ll for current betas
                      lp_olds,
                      theta_temp,          // container of current betas for all i
                      XX,
                      PP,
                      AA,
                      nalts,
                      sumpxs,
                      ntasks,
                      xfr,xto,lfr,lto,
                      maxspents,
                      p,N,
                      vectorise(MU),      // Mean Prior
                      Lprior,             // VarCov Prior (chol)
                      stay,               // tracking rejections
                      tunes,              // i-level tuning parameters
                      cores);       
    
    // tuning stuff ..................
    tunecounter+=1; // update counter within tune interval
    
    //rejection rate in window
    if(tunecounter>=tuneinterval){
      
      //just for progress output
      currentRR=mean(stay/tunecounter); 
      //tune within tune range
      if( (ir>=tunestart) & (ir<(tunestart+tunelength))){
        mh_tuner(tunes,stay/tunecounter);
      }
      
      //reset
      tunecounter=1;
      stay_total+=stay;
      stay.fill(0);
    }
    // end of loop
    
    
    // save draws  ..................
    if((ir+1)%keep==0){
      mkeep = (ir+1)/keep-1;
      thetaDraw.slice(mkeep)      = theta_temp;
      MUDraw.row(mkeep)           = trans(vectorise(MU,0));
      SIGMADraw.slice(mkeep)      = SIGMA;
      loglike(mkeep)              = sum(ll_olds);
      logpost(mkeep)              = sum(ll_olds)+sum(lp_olds);
      rrate(mkeep)                = currentRR;
    }
    
    // display progress  ..................
    if((ir+1)%progressinterval==0){
      infoMcmcTimerRRLL(ir,R,currentRR,sum(ll_olds));
    }
    
  }
  endMcmcTimer();
  
  //return draws  ..................
  return List::create(
    Named("thetaDraw")      = thetaDraw,
    Named("MUDraw")         = MUDraw,
    Named("SIGMADraw")      = SIGMADraw,
    Named("loglike")        = loglike,
    Named("logpost")        = logpost,
    Named("reject")         = rrate,
    Named("RRs")            = stay_total/R);
}















//////////////////////////////////////////
// demand estimation stuff
//////////////////////////////////////////


//' @export
// [[Rcpp::export]]
vec vd_demand(arma::vec psi, double gamma, double E, vec prices) {
  arma::vec rho1    = prices / (psi);
  arma::uvec rho_ord= arma::stable_sort_index(rho1);
  rho1              = rho1.elem(rho_ord);
  
  int J = prices.size();
  int rhos = rho1.size();
  
  arma::vec rho(rhos+2);
  rho(0) = 0;
  rho(rhos+1) = R_PosInf;
  for (int i = 1; i <rhos+1; ++i){
    rho(i) = rho1(i-1);
  }
  
  arma::vec psi_ord     = psi.elem(rho_ord);
  arma::vec prices_ord  = prices.elem(rho_ord);
  double a = gamma*E;
  double b = gamma;
  double k = 0;
  double z = a/b;
  
  while( (z<=rho[k]) | (z>rho[k+1]) ){
    a=a+arma::as_scalar(prices_ord[k]);
    b=b+arma::as_scalar(psi_ord[k]);
    z=a/b;
    k=k+1;
  }
  arma::vec x = (psi_ord*z-prices_ord) /  (gamma*prices_ord);
  for (int i = k; i <J; ++i){
    x(i) = 0;
  }
  arma::vec X =x.elem(arma::stable_sort_index(arma::conv_to<arma::vec>::from(rho_ord)));
  return(X);
}






// ------------- - - - - - - - -
// standard volumetric - ev error
// ------------- - - - - - - - - 


//' @export
//[[Rcpp::export]]
List des_dem_vdm(vec const& PP,   //price (vectorised)
                 mat const& AA,   //attr-lvls (vectorised)
                 uvec const& nalts,
                 ivec const& ntasks,  
                 ivec const& xfr,  //indexing variables
                 ivec const& xto,  
                 ivec const& lfr,  
                 ivec const& lto,
                 ivec const& tlens,
                 cube const& thetaDraw, 
                 int cores=1){
  
  // vd2dem
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = xfr.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    

    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
          
          //paras
          arma::vec theta = thetaDraw.slice(ir).col(n);
          arma::vec beta  = theta(arma::span(0,p-4));
          double sigma    = exp(theta(p-3));
          double gamma    = exp(theta(p-2));
          double E        = exp(theta(p-1));
          
          //compute psi
          //arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +sigma*randn(nalt)); 
          arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
                              revd0(nalt, sigma)); 
          
          //find opt demand
          arma::vec dem = vd_demand(psi, gamma, E, prcs);
          
          demcontainer.col(ir)=dem;
          //Xd(span(xpick,xpick+nalt-1),ir)   = dem;
      }

      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


// ------------- - - - - - - - -
// standard volumetric - Normal error
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List des_dem_vdmn(vec const& PP,
                  mat const& AA,
                  uvec const& nalts,
                  ivec const& ntasks,  
                  ivec const& xfr,
                  ivec const& xto,  
                  ivec const& lfr,  
                  ivec const& lto,
                  ivec const& tlens,
                  cube const& thetaDraw, 
                  int cores=1){
  
  // vd2dem
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = xfr.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*randn(nalt)); 
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
        
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      } 
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


// ------------- - - - - - - - -
//standard volumetric - error argument
// ------------- - - - - - - - - 


//' @export
//[[Rcpp::export]]
List der_dem_vdm(vec const& PP,
                 mat const& AA,
                 uvec const& nalts,
                 ivec const& ntasks,  
                 ivec const& xfr,
                 ivec const& xto,  
                 ivec const& lfr,  
                 ivec const& lto,
                 ivec const& tlens,
                 cube const& thetaDraw, 
                 vec const& epsilon,
                 int cores=1){
  
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*epsilon(span(xpick,xpick+nalt-1)) ); 
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}



// ------------- - - - - - - - -
//conjunctive volumetric - normal error
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List des_dem_vdm_screen(vec const& PP,
                             mat const& AA,
                             mat const& AAf,
                             uvec const& nalts,
                             uvec const& tlens,
                             ivec const& ntasks,  
                             ivec const& xfr,
                             ivec const& xto,  
                             ivec const& lfr,  
                             ivec const& lto,
                             cube const& thetaDraw,
                             cube const& tauDraw, 
                             int cores=1){
  
  // vdsr2dem
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        arma::vec taui  = tauDraw.slice(ir).col(n);
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
                            sigma*randn(nalt)); 
        
        //screen
        psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;

        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }  
      xpick+=nalt;
    }
    
  }
  return(XdL);
}



// ------------- - - - - - - - -
// conjunctive volumetric - error argument
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List der_dem_vdm_screen(vec const& PP,
                             mat const& AA,
                             mat const& AAf,
                             uvec const& nalts,
                             uvec const& tlens,  //
                             ivec const& ntasks,  
                             ivec const& xfr,  //
                             ivec const& xto,  
                             ivec const& lfr,  
                             ivec const& lto,
                             cube const& thetaDraw,
                             cube const& tauDraw,
                             vec const& epsilon,
                             int cores=1){
  // vdsr2dem
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        arma::vec taui  = tauDraw.slice(ir).col(n);
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*epsilon(span(xpick,xpick+nalt-1)) ); 
        
        //screen
        psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


// ------------- - - - - - - - -
// conjunctive volumetric - normal error
// ------------- - - - - - - - - 
// screening includes screening on price tags

//' @export
//[[Rcpp::export]]
List des_dem_vdm_screenpr(vec const& PP,
                             mat const& AA,
                             mat const& AAf,
                             uvec const& nalts,
                             uvec const& tlens,
                             ivec const& ntasks,  
                             ivec const& xfr,
                             ivec const& xto,  
                             ivec const& lfr,  
                             ivec const& lto,
                             cube const& thetaDraw,
                             cube const& tauDraw, 
                             mat const& tau_pr_Draw,
                             int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        arma::vec taui  = tauDraw.slice(ir).col(n);
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*randn(nalt)); 
        
        //screen
        psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
        psi.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
                 
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}

// ------------- - - - - - - - -
// conjunctive volumetric - error argument
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List der_dem_vdm_screenpr(vec const& PP,
                             mat const& AA,
                             mat const& AAf,
                             uvec const& nalts,
                             uvec const& tlens,
                             ivec const& ntasks,  
                             ivec const& xfr,
                             ivec const& xto,  
                             ivec const& lfr,  
                             ivec const& lto,
                             cube const& thetaDraw,
                             cube const& tauDraw,
                             mat const& tau_pr_Draw,
                             vec const& epsilon,
                             int cores=1){
  
  // vdsr2dem
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-4));
        double gamma    = exp(theta(p-2));
        double E        = exp(theta(p-1));
        double sigma    = exp(theta(p-3));
        arma::vec taui  = tauDraw.slice(ir).col(n);
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*epsilon(span(xpick,xpick+nalt-1)) ); 
        
        //screen
        psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
        psi.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
          
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


// ******************************************* -------------------
// with set size variation
// ******************************************* -------------------


// ------------- - - - - - - - -
// standard volumetric - setsize - EV error
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List des_dem_vdm_ss(vec const& PP,
                             mat const& AA,
                             uvec const& nalts,
                             vec const& sumpxs,  
                             ivec const& ntasks,  
                             ivec const& xfr,
                             ivec const& xto,  
                             ivec const& lfr,  
                             ivec const& lto,
                             ivec const& tlens,
                             cube const& thetaDraw, 
                             int cores=1){
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = xfr.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
  
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-5));
        
        double E        = exp(theta(p-1));
        double gamma    = exp(theta(p-2));
        double sigma    = exp(theta(p-3));
        double xi       = exp(theta(p-4));
        
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          revd0(nalt, sigma) - 
          log(nalt*xi+1) ); 
        
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }  
      xpick+=nalt;
    }

  }
  return(XdL);
}







// ------------- - - - - - - - -
// standard volumetric - setsize - error argument
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List der_dem_vdm_ss(vec const& PP,
                    mat const& AA,
                    uvec const& nalts,
                    uvec const& xlens,  
                    uvec const& tlens,
                    ivec const& ntasks,  
                    ivec const& xfr,
                    ivec const& xto,  
                    ivec const& lfr,  
                    ivec const& lto,
                    cube const& thetaDraw, 
                    vec const& epsilon,
                    int cores=1){
  // vd2dem3
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-5));
        
        double E        = exp(theta(p-1));
        double gamma    = exp(theta(p-2));
        double sigma    = exp(theta(p-3));
        double xi       = exp(theta(p-4));
        
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +sigma*epsilon(span(xpick,xpick+nalt-1)) - log(nalt*xi+1) ); 
        
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}



// ------------- - - - - - - - -
// standard volumetric - setsize quadratic - EVerror 
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List des_dem_vdm_ssq(vec const& PP,
                     mat const& AA,
                     uvec const& nalts,
                     vec const& sumpxs,  
                     ivec const& ntasks,  
                     ivec const& xfr,
                     ivec const& xto,  
                     ivec const& lfr,  
                     ivec const& lto,
                     ivec const& tlens,
                     cube const& thetaDraw, 
                     int cores=1){
  // vd2Qdem3

  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-6));
        
        double E        = exp(theta(p-1));
        double gamma    = exp(theta(p-2));
        double sigma    = exp(theta(p-3));
        double xi       = exp(theta(p-4));
        double tau      = exp(theta(p-5));
        
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          revd0(nalt, sigma)- 
          log(nalt*xi+1+tau*(nalt^2)) ); 
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }   
      xpick+=nalt;
    }
    
  }
  return(XdL);
}


// ------------- - - - - - - - -
// standard volumetric - setsize quadratic - error argument
// ------------- - - - - - - - - 

//' @export
//[[Rcpp::export]]
List der_dem_vdm_ssq(vec const& PP,
                              mat const& AA,
                              uvec const& nalts,
                              uvec const& xlens,  
                              uvec const& tlens,
                              ivec const& ntasks,  
                              ivec const& xfr,
                              ivec const& xto,  
                              ivec const& lfr,  
                              ivec const& lto,
                              cube const& thetaDraw, 
                              vec const& epsilon,
                              int cores=1){
  // vd2Qdem3
  
  //dimensions
  int R=thetaDraw.n_slices;
  int xdim = PP.n_rows;
  int N = tlens.n_elem;
  int p = thetaDraw.n_rows;
  
  //output init
  Rcpp::List XdL(xdim);
  
  //start timer
  startTimer();
  
  //resp-level
  for(int n=0; n<N; n++){
    infoTimer(n,N);
    
    int ntask = tlens(n);
    int xpick = xfr(n);
    
    //task-level
    for(int tt=0; tt<ntask; tt++){
      Rcpp::checkUserInterrupt();
      
      int nalt = nalts(tt);
      
      //temp storage
      mat demcontainer(nalt,R, fill::zeros);
      
      arma::vec prcs= PP(span(xpick,xpick+nalt-1));
      
      //draw-level
      omp_set_num_threads(cores);
      #pragma omp parallel for schedule(static)
      for (int ir = 0; ir <R; ++ir){
        
        //paras
        arma::vec theta = thetaDraw.slice(ir).col(n);
        arma::vec beta  = theta(arma::span(0,p-6));
        
        double E      = exp(theta(p-1));
        double gamma    = exp(theta(p-2));
        double sigma    = exp(theta(p-3));
        double xi       = exp(theta(p-4));
        double tau       = exp(theta(p-5));
        
        
        //compute psi
        arma::vec psi = exp(AA(span(xpick,xpick+nalt-1),span::all)*beta +
          sigma*epsilon(span(xpick,xpick+nalt-1)) - 
          log(nalt*xi+1+tau*(nalt^2)) ); 
        //revd0(nalt, sigma)
        
        //find opt demand
        arma::vec dem = vd_demand(psi, gamma, E, prcs);
        demcontainer.col(ir)=dem;
      }
      
      for(int k=0; k<nalt; k++){
        XdL(k+xpick)=demcontainer.row(k);
      }  
      xpick+=nalt;
    }
    
  }
  return(XdL);
}

