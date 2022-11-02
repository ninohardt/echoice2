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


 // [[Rcpp::export]]
 vec revd(int n=100, double loc=0, double scale=1){
   return(loc-scale*log(-log(runif(n))));
 }


// [[Rcpp::export]]
vec revdx(vec locs, 
          vec scales){
  
  int n = locs.n_elem;
  vec out(n);
  out = locs-scales%log(-log(randu(n)));
  
  return(out);
}



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
  int k = sum(as_scalar(randu(1))>cumsum(probs));
  return(k);
}



 // [[Rcpp::export]]
 mat riwish(double v, mat S){
   return(iwishrnd( S, v ));
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



 // [[Rcpp::export]]
 mat rmvn(int n, vec mu, mat sig){
   return(trans(mvnrnd(mu, sig, n)));
 }









///////////////////////////////////////////////
// DD compensatory
///////////////////////////////////////////////


//log-likelihood

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
     maxpaids(n) = max(sign(XX(span(xfr(n),xto(n)))%PP(span(xfr(n),xto(n)))));
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



 //[[Rcpp::export]]
 arma::field<arma::vec> dddem(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }





 //[[Rcpp::export]]
 arma::field<arma::vec> ddsrdem(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }



 //[[Rcpp::export]]
 arma::field<arma::vec> ddsrprdem(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       } 
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


///////////////////////////////////////////////
// DD Demand-Prob
///////////////////////////////////////////////



 //[[Rcpp::export]]
 arma::field<arma::vec> ddprob(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }





 //[[Rcpp::export]]
 arma::field<arma::vec> ddsrprob(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }
       
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }



 //[[Rcpp::export]]
 arma::field<arma::vec> ddsrprprob(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }
       
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }











/////////////////////////////////////// 
// VD Compensatory - EV error
///////////////////////////////////////


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
   theta_temp.row(p-4) += -3;
   
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






 // [[Rcpp::export]]
 vec vdss_LL(mat const&Theta,
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
     ll_olds(n)= vdl_ss(Theta.col(n),
             nalts(span(lfr(n),lto(n))),
             sumpxs(span(lfr(n),lto(n))), 
             XX(span(xfr(n),xto(n))), 
             PP(span(xfr(n),xto(n))), 
             AA(span(xfr(n),xto(n)),span::all), 
             ntasks(n), p);
   }
   
   return(ll_olds);
 }


 // [[Rcpp::export]]
 mat vdss_LLs(cube const&THETAS,
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
       vdss_LL(THETAS.slice(r),
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










//////////////////////////////////////////
// demand estimation stuff
//////////////////////////////////////////



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
   
   while( (z<=rho[k]) || (z>rho[k+1]) ){
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



 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdm(vec const& PP,   //price (vectorised)
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
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


// ------------- - - - - - - - -
// standard volumetric - Normal error
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdmn(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       } 
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


// ------------- - - - - - - - -
//standard volumetric - error argument
// ------------- - - - - - - - - 



 //[[Rcpp::export]]
 arma::field<arma::vec> der_dem_vdm(vec const& PP,
                                    mat const& AA,
                                    uvec const& nalts,
                                    ivec const& ntasks,  
                                    ivec const& xfr,
                                    ivec const& xto,  
                                    ivec const& lfr,  
                                    ivec const& lto,
                                    ivec const& tlens,
                                    cube const& thetaDraw, 
                                    mat const& epsilon,
                                    int cores=1){
   
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
           sigma*epsilon(span(xpick,xpick+nalt-1), ir) ); 
         
         //find opt demand
         arma::vec dem = vd_demand(psi, gamma, E, prcs);
         demcontainer.col(ir)=dem;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }



// ------------- - - - - - - - -
//conjunctive volumetric - normal error
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdm_screen(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }  
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }



// ------------- - - - - - - - -
// conjunctive volumetric - error argument
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> der_dem_vdm_screen(vec const& PP,
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
                                           mat const& epsilon,
                                           int cores=1){
   // vdsr2dem
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
           sigma*epsilon(span(xpick,xpick+nalt-1),ir) ); 
         
         //screen
         psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
         
         //find opt demand
         arma::vec dem = vd_demand(psi, gamma, E, prcs);
         demcontainer.col(ir)=dem;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
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


 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdm_screenpr(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }

// ------------- - - - - - - - -
// conjunctive volumetric - error argument
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> der_dem_vdm_screenpr(vec const& PP,
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
                                             mat const& epsilon,
                                             int cores=1){
   
   // vdsr2dem
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
           sigma*epsilon(span(xpick,xpick+nalt-1), ir) ); 
         
         //screen
         psi.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))*=0;
         psi.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
         
         //find opt demand
         arma::vec dem = vd_demand(psi, gamma, E, prcs);
         demcontainer.col(ir)=dem;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
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


 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdm_ss(vec const& PP,
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
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = xfr.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }  
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }







// ------------- - - - - - - - -
// standard volumetric - setsize - error argument
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> der_dem_vdm_ss(vec const& PP,
                                       mat const& AA,
                                       uvec const& nalts,
                                       ivec const& ntasks,  
                                       ivec const& xfr,
                                       ivec const& xto,  
                                       ivec const& lfr,  
                                       ivec const& lto,
                                       uvec const& tlens,
                                       cube const& thetaDraw, 
                                       mat const& epsilon,
                                       int cores=1){
   // vd2dem3
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
           sigma*epsilon(span(xpick,xpick+nalt-1), ir) - log(nalt*xi+1) ); 
         
         
         //find opt demand
         arma::vec dem = vd_demand(psi, gamma, E, prcs);
         demcontainer.col(ir)=dem;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }



// ------------- - - - - - - - -
// standard volumetric - setsize quadratic - EVerror 
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> des_dem_vdm_ssq(vec const& PP,
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
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


// ------------- - - - - - - - -
// standard volumetric - setsize quadratic - error argument
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> der_dem_vdm_ssq(vec const& PP,
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
                                        mat const& epsilon,
                                        int cores=1){
   // vd2Qdem3
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
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
           sigma*epsilon(span(xpick,xpick+nalt-1), ir) - 
           log(nalt*xi+1+tau*(nalt^2)) ); 
         //revd0(nalt, sigma)
         
         //find opt demand
         arma::vec dem = vd_demand(psi, gamma, E, prcs);
         demcontainer.col(ir)=dem;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       }  
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }




// ------------- - - - - - - - -
// Screening probabilities (w/o price)
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> ec_screen_prob_cpp( vec const& PP,
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
   //int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
       //temp storage
       mat demcontainer(nalt,R, fill::zeros);
       
       //draw-level
       omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
       for (int ir = 0; ir <R; ++ir){
         arma::vec pr(nalt, fill::zeros);
         arma::vec taui  = tauDraw.slice(ir).col(n);
         
         pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))+=1;
         demcontainer.col(ir)=pr;
       }
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       }   
       
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


// ------------- - - - - - - - -
// Screening probabilities (w/ price)
// ------------- - - - - - - - - 


 //[[Rcpp::export]]
 arma::field<arma::vec> ec_screenpr_prob_cpp( vec const& PP,
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
   //int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = tlens(n);
     int xpick = xfr(n);
     uvec const& nalts_i=nalts.subvec(lfr(n),lto(n));
     
     //task-level
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       int nalt = nalts_i(tt);
       
       //temp storage
       mat demcontainer(nalt,R, fill::zeros);
       
       arma::vec prcs = PP(span(xpick,xpick+nalt-1));
       
       //draw-level
       omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
       for (int ir = 0; ir <R; ++ir){
         
         arma::vec pr(nalt);
         pr.fill(0);
         arma::vec taui  = tauDraw.slice(ir).col(n);
         
         pr.elem(find((AAf(span(xpick,xpick+nalt-1),span::all)*taui)>0.01))+=1;
         pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))+=1;
         pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))=sign(pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  )))));
         
         demcontainer.col(ir)=pr;
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       }
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }


















///////////////////////////////////////////////
// Discrete Demand - Screening w/price
///////////////////////////////////////////////



 // [[Rcpp::export]]
 double ddlpr(arma::vec const& theta, 
              double tau_pr,
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
       ab+=(-beta_p*p);
       
       //screening
       if(p<=exp(tau_pr)){
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




void draw_ddpr_RWMH( arma::vec& ll_olds,       // vector of current log-likelihoods
                     arma::vec& lp_olds,       // vectors of lp's, just for tracking 
                     arma::mat& theta_temp,    // container of current betas for all i
                     vec const& tau_prs,
                     vec const& XX,            // data
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
    llnew = ddlpr(theta_cand,
                  tau_prs(n),
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




void draw_dd_taupr1( vec& ll_olds,
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
      
      
      double llnew = ddlpr(theta_temp.col(n),
                           tau_pr_cand,
                           nalts(span(lfr(n),lto(n))),
                           XX(span(xfr(n),xto(n))), 
                           PP(span(xfr(n),xto(n))), 
                           AA(span(xfr(n),xto(n)),span::all), 
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



 // [[Rcpp::export]]
 List loop_ddpr_RWMH(  vec const& XX, 
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
   
   
   vec maxpaids(N);
   for(int n=0; n<N; n++){
     maxpaids(n) = max(sign(XX(span(xfr(n),xto(n)))%PP(span(xfr(n),xto(n)))));
   }
   
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
       ddlpr(theta_temp.col(n),
             tau_prs(n),
             nalts(span(lfr(n),lto(n))),
             XX(span(xfr(n),xto(n))), 
             PP(span(xfr(n),xto(n))), 
             AA(span(xfr(n),xto(n)),span::all), 
             ntasks(n), p );
   }
   
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
     draw_ddpr_RWMH(ll_olds,            // ll for current betas
                    lp_olds,
                    theta_temp,          // container of current betas for all i
                    tau_prs,
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
     
     if(ir>1000){
       
       draw_dd_taupr1( ll_olds,
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
       
       MUDraw.row(mkeep)           = trans(vectorise(MU,0));
       
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
     Named("tau_pr_draw")    = tau_pr_draw,
     Named("prscreenMuSigDraw")  = pricescreenPriorDraw,
     Named("loglike")        = loglike,
     Named("logpost")        = logpost,
     Named("reject")         = rrate,
     Named("RRs")            = stay_total/R);
 }
















///////////////////
///////////////////



//price screening only


 //[[Rcpp::export]]
 arma::field<arma::vec> ddprdem(vec const& PP,
                                mat const& AA,
                                uvec const& nalts,
                                uvec const& tlens,
                                ivec const& ntasks,  
                                ivec const& xfr,
                                ivec const& xto,  
                                ivec const& lfr,  
                                ivec const& lto,
                                cube const& thetaDraw,
                                mat const& tau_pr_Draw,
                                int cores=1){
   
   //dimensions
   int R=thetaDraw.n_slices;
   int xdim = PP.n_rows;
   int N = tlens.n_elem;
   int p = thetaDraw.n_rows;
   
   //output init
   arma::field<arma::vec>XdL(xdim);
   
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
         
         pr.elem(find((prcs>exp(  tau_pr_Draw(n,ir)  ))))*=0;
         
         //multinomial draw
         int pick_draw = rmuno2(pr);
         
         //if not outside good, choose inside good
         if(pick_draw!=nalt){
           demcontainer(pick_draw,ir)=1;
         }
         
       }
       
       for(int k=0; k<nalt; k++){
         XdL(k+xpick)=trans(demcontainer.row(k));
       } 
       xpick+=nalt;
     }
     
   }
   return(XdL);
 }




/////////////////// Experimental ///////////////////


// utility functions
int mod1(int a, int n){
  return(a - floor(a/n)*n);
}  


// from sequential id to permutation

// [[Rcpp::export]]
ivec index_id2alt(int id, ivec nalts) {
  int p = nalts.n_elem;
  ivec out(p);
  
  for(int k=0; k<p; k++){
    out(k) = floor(mod1(id/(prod(nalts.head(k))),nalts(k)));
  }
  return(out);
}



//implementation isn't super-efficient with memory, but should be fast


 //[[Rcpp::export]]
 arma::field<arma::vec> ddprdemseq1(arma::field<arma::vec> PPfield,
                                    arma::field<arma::mat> AAfield,
                                    arma::field<arma::uvec> naltfield,
                                    arma::field<arma::ivec> ntaskfield,
                                    arma::field<arma::ivec> xfrfield,
                                    arma::ivec pvecs,
                                    arma::imat pfrto,
                                    arma::field<arma::ivec> secpick,
                                    cube const& thetaDraw,
                                    mat const& tau_pr_Draw,
                                    int cores=1){
   
   int nstage = PPfield.n_elem;
   int xdim_all=0;
   for(int k=0; k<nstage; k++){
     xdim_all+=PPfield(k).n_rows;
   }
   
   //dimensions
   int R = thetaDraw.n_slices;
   int N = ntaskfield(0).n_elem;
   
   //output init
   arma::field<arma::vec>XdL(xdim_all);
   
   //start timer
   startTimer();
   
   //resp-level
   for(int n=0; n<N; n++){
     // Rcpp::Rcout << n << std::endl;    
     infoTimer(n,N);
     
     int ntask = ntaskfield(0)(n);
     
     mat prevstagespent(ntask, R, fill::zeros);
     
     //stage-level
     for(int ss=0; ss<nstage; ss++){
       
       int xpick = xfrfield(ss)(n);
       int p = pvecs(ss);
       
       //task-level - each task may have its own previous-stage-spent information
       //              which we have to pick from stagespents
       for(int tt=0; tt<ntask; tt++){
         Rcpp::checkUserInterrupt();
         
         int nalt = naltfield(ss)(tt);
         
         //temp storage
         mat demcontainer(nalt, R, fill::zeros);
         ivec nalt_space = linspace<ivec>(0, nalt-1); 
         arma::vec prcs = PPfield(ss).subvec(xpick,xpick+nalt-1);
         
         //draw-level
         omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
         for (int ir = 0; ir <R; ++ir){
           
           //paras - stage specific
           arma::vec theta = thetaDraw.slice(ir).col(n).subvec(pfrto(ss,0),pfrto(ss,1));		
           arma::vec beta  = theta(arma::span(0,p-2));
           double beta_p   = exp(theta(p-1));
           
           //precomputes
           arma::vec ab = (AAfield(ss)(span(xpick,xpick+nalt-1),span::all)) * beta - prcs*beta_p;
           arma::vec pr = exp(ab)/(1+sum(exp(ab)));
           
           pr.elem(find(prcs>(exp(tau_pr_Draw(n,ir)) - prevstagespent(tt, ir) )  ))*=0;
           
           //multinomial draw
           int pick_draw = rmuno2(pr);
           
           //if not outside good, choose inside good
           if(pick_draw!=nalt){
             demcontainer(pick_draw,ir)=1;
           }
           
           //update spent on previous stages info for each draw of theta
           prevstagespent(tt, ir)+=sum(demcontainer.col(ir) % prcs);
           
         }
         for(int k=0; k<nalt; k++){
           XdL( secpick(ss)(k+xpick) ) = trans( demcontainer.row(k) );
         } 
         xpick+=nalt;
         
       }//t loop
       
     }//stage loop
     
   }
   
   return(XdL);
 }





//simultaneous choice with budget


vec pick2demvec(int pick, 
                int nalt){
  vec out(nalt,fill::zeros);
  
  if(pick<nalt){
    out(pick)=1;
  }
  return(out);
}




 //[[Rcpp::export]]
 arma::field<arma::vec> ddprdemsimu1(arma::field<arma::vec>  PPfield,
                                     arma::field<arma::mat>  AAfield,
                                     arma::field<arma::uvec> naltfield,
                                     arma::field<arma::ivec> ntaskfield,
                                     arma::field<arma::ivec> xfrfield,
                                     arma::ivec pvecs,
                                     arma::imat pfrto,
                                     arma::field<arma::ivec> secpick,
                                     cube const& thetaDraw,
                                     mat const& tau_pr_Draw,
                                     int cores=1){
   
   int nstage = PPfield.n_elem;
   int xdim_all=0;
   for(int k=0; k<nstage; k++){
     xdim_all+=PPfield(k).n_rows;
   }
   
   //dimensions
   int R = thetaDraw.n_slices;
   int N = ntaskfield(0).n_elem;
   
   //output init
   arma::field<arma::vec> XdL(xdim_all);
   
   //start timer
   startTimer();
   
   //temp var tracking alternatives in each stage
   ivec nalt_temp(nstage);
   
   //resp-level ------------------------------------------------------
   for(int n=0; n<N; n++){
     infoTimer(n,N);
     
     int ntask = ntaskfield(0)(n);
     ivec xpicks(nstage);
     for(int ss=0; ss<nstage; ss++){
       xpicks(ss) = xfrfield(ss)(n);
     }
     
     //task-level -----------------------------------------------------
     for(int tt=0; tt<ntask; tt++){
       Rcpp::checkUserInterrupt();
       
       for(int ss=0; ss<nstage; ss++){
         nalt_temp(ss) = naltfield(ss)(tt);
       }
       
       arma::field<arma::vec> xbL(nstage);
       arma::field<arma::mat> demcontainers(nstage);
       
       //setup temp storage of draws
       for(int ss=0; ss<nstage; ss++){
         demcontainers(ss) = zeros(naltfield(ss)(tt),R);
       }
       
       //draw-level ------------------------------------------------------
       omp_set_num_threads(cores);
#pragma omp parallel for schedule(static)
       for (int ir = 0; ir <R; ++ir){
         
         //stage-level ------------------------------------------------------
         for(int ss=0; ss<nstage; ss++){
           
           //temp storage
           //ivec nalt_space = linspace<ivec>(0, nalt_temp(ss)-1); 
           arma::vec prcs  = PPfield(ss).subvec(xpicks(ss),xpicks(ss)+nalt_temp(ss)-1);
           
           int p = pvecs(ss);
           
           //paras - stage specific
           arma::vec theta = thetaDraw.slice(ir).col(n).subvec(pfrto(ss,0),pfrto(ss,1));		
           arma::vec beta  = theta(arma::span(0,p-2));
           double beta_p   = exp(theta(p-1));
           
           //precomputes
           vec nulli(1);
           nulli(0)=0;
           xbL(ss)= join_cols( vectorise((AAfield(ss)(span(xpicks(ss),xpicks(ss)+nalt_temp(ss)-1),span::all)) * beta - prcs*beta_p), nulli);
           
         } //stage loop =====================================================
         
         int ncombs = prod(nalt_temp+1); //include outside option
         vec logcombprops(ncombs, fill::zeros);
         vec combprops(ncombs, fill::zeros);
         
         //combinations -----
         for(int kk=0; kk<ncombs; kk++){
           ivec tempid = index_id2alt(kk, nalt_temp+1);
           double Etot=0;  //total cost for combination
           
           //combination probabilities
           for(int ss=0; ss<nstage; ss++){
             
             if(tempid(ss)<=nalt_temp(ss)){
               Etot+=PPfield(ss)(tempid(ss)); //0 cost osg
             }
             logcombprops(kk)+=xbL(ss)(tempid(ss));
           }
           
           //remove combinations beyond budget
           if(Etot<exp(tau_pr_Draw(n,ir))){
             combprops(kk)=exp(logcombprops(kk));
           }
           
         } // combinations  =====
         
         //probabilities of combinations          
         if(sum(combprops)>0){
           combprops = combprops / sum(combprops);
         }
         
         //draw
         int pick_draw = rmuno2(combprops);
         ///
         ivec tempid = index_id2alt(pick_draw, nalt_temp+1);
         
         for(int ss=0; ss<nstage; ss++){
           // demcontainers(ss).col(ir)+=tempid(ss);// = pick2demvec( tempid(ss), nalt_temp(ss) );
           Rcout << tempid(ss);
           
           demcontainers(ss).col(ir)+=pick2demvec( tempid(ss), nalt_temp(ss) );
           
         }
         
       } //draw-level =====================================================
       
       //each product has its own draw list
       for(int ss=0; ss<nstage; ss++){
         for(int k=0; k<nalt_temp(ss); k++){
           XdL( secpick(ss)(k+xpicks(ss)) ) = trans(demcontainers(ss).row(k));
         } 
         xpicks(ss)+=nalt_temp(ss);
       }
       
     } //t loop =====================================================
     
   } //n loop =====================================================
   
   return(XdL);
 }
