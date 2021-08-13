#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List inseq(arma::mat M, bool adjust=true)
{
  int i, m, trun;
  arma::mat mu=mean(M);
  //center the rows
  M.each_row() -= mu;
  int n=M.n_rows, p=M.n_cols;
  //Dtm is the vector of det(Sig)'s
  NumericVector Dtm;
  //gam_0 and gam_1 are for gam_2m and gam_2m+1, resp.
  arma::mat gam0(p,p), gam1(p,p), Gam(p,p), Sig(p,p), Sig1(p,p), Gamadj(p,p), eigvec(p,p);
  //for adjustment, set initial increment in Gam=0
  Gamadj.zeros();
  //store the eigenvalues and eigenvectors of each Gam
  arma::vec eigval(p),eigvalneg(p);
  double dtm;
  int sn= n/2; 
  for (m=0; m<n/2; m++)
  {
    gam0.zeros(); gam1.zeros();
    //calculate gam_2m (gam0) and gam_2m+1 (gam1)
    for(i=0; i<(n-2*m);i++) gam0+=trans(M.row(i))*M.row(i+2*m);
    gam0=gam0/n;
    for(i=0; i<(n-2*m-1);i++) gam1+=trans(M.row(i))*M.row(i+2*m+1);
    gam1=gam1/n;
    //Gam_m=gam_2m+gam_2m+1, then symmetrize
    Gam=gam0+gam1; Gam=(Gam+Gam.t())/2;
    
    if (m==0) Sig=-gam0+2*Gam;
    else Sig=Sig+2*Gam;
    
    if (eig_sym(Sig)(0)>0)
    {
      sn=m;
      break;
    }
  }
  if (sn>n/2-1) 
  {
    stop("Not enough samples.");
  }
  Dtm=det(Sig);
  for (m=sn+1; m<n/2; m++)
  {
    gam0.zeros(); gam1.zeros();
    //calculate gam_2m (gam0) and gam_2m+1 (gam1)
    for(i=0; i<(n-2*m);i++) gam0+=trans(M.row(i))*M.row(i+2*m);
    gam0=gam0/n;
    for(i=0; i<(n-2*m-1);i++) gam1+=trans(M.row(i))*M.row(i+2*m+1);
    gam1=gam1/n;
    //Gam_m=gam_2m+gam_2m+1, then symmetrize
    Gam=gam0+gam1; Gam=(Gam+Gam.t())/2;
    
    //Sig_m=Sig_m-1+2Gam_m
    Sig1=Sig+2*Gam;
    dtm=det(Sig1);
    //if dtm1>dtm, continue
    if (dtm<=Dtm(m-sn-1)) break;
    //update Sig
    Sig=Sig1;
    //record dtm
    Dtm.push_back(dtm);

    //to adjust the original Sig, subtract the negative part of Gam
    if (adjust) 
    {
      //calculate eigenvalues and eigenvectors of Gam
      eig_sym(eigval,eigvec,Gam);
      eigvalneg=eigval;
      eigvalneg.elem(find(eigvalneg>0)).zeros();
      Gamadj-=eigvec*diagmat(eigvalneg)*eigvec.t();
    }
  }
  trun = Dtm.size()-1+sn;
  List res; res["Sig"]=Sig; res["Dtm"]=Dtm; res["trunc"]= trun; res["sn"]=sn;
  if (adjust) res["Sigadj"]=Sig+2*Gamadj;
  return res;
}
