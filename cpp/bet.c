/*
######################################################################
#
# nli.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 6/25/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the computation of node level
# indices (NLIs).
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "bet.h"


SEXP betweenness_partial(SEXP mat, SEXP sn, SEXP sm, SEXP smeasure, SEXP sst, SEXP send)
/*
Compute betweenness (and some related measures) for the network in mat.  If sprecomp==TRUE, then sgd, ssigma, and spred are taken to hold geodesic distances, path counts, and predecessor lists (as returned by geodist_R or geodist_val_r); else, these are computed on the fly (which, BTW, saves memory, though it prohibits reuse).  If signoreevals==TRUE, then edge values are not used when computing paths (irrelevant if called with precomputed geodesics).
*/
{
  int n,i,j,k,wv,precomp,*npred,ignoreeval,measure,pc=0,st,end;
  double *gd, *sigma,*bet,*delta;
  element **pred,*w,*v;
  snaNet *g;
  SEXP sbet,lp,vp;

  /*Coerce inputs*/
  PROTECT(mat=coerceVector(mat,REALSXP)); pc++;
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(sm=coerceVector(sm,INTSXP)); pc++;
 // PROTECT(sprecomp=coerceVector(sprecomp,INTSXP)); pc++;
  PROTECT(smeasure=coerceVector(smeasure,INTSXP)); pc++;
  //PROTECT(signoreeval=coerceVector(signoreeval,INTSXP)); pc++;
  n=INTEGER(sn)[0];
  measure=INTEGER(smeasure)[0];
  ignoreeval = 1;
  st = INTEGER(sst)[0];
  end = INTEGER(send)[0];
  
  if (st < 0 || end > n)
    return NULL;


  /*Allocate memory*/
  PROTECT(sbet=allocVector(REALSXP,n)); pc++;
  npred=(int *)R_alloc(n,sizeof(int));
  pred=(element **)R_alloc(n,sizeof(element *));
  gd=(double *)R_alloc(n,sizeof(double));
  sigma=(double *)R_alloc(n,sizeof(double));
  delta=(double *)R_alloc(n,sizeof(double));
  bet=REAL(sbet);
  
  /*Set up stuff*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(mat),INTEGER(sn),INTEGER(sm));
  PutRNGstate();
  for(i=0;i<n;i++)
    bet[i]=0.0;

  /*Calculate betweenness*/
  for(i=st;i<end;i++){
    R_CheckUserInterrupt();
    /*Get geodesic information*/
    /*Compute on the fly*/
    if(ignoreeval)
      spsp(i,g,gd,sigma,pred,npred,1);
    else
      spsp_val(i,g,gd,sigma,pred,npred,1);
    
    /*Accumulate betweenness incremements*/
    switch(measure){
      case BETSTANDARD:    /*"Standard" form betweenness (a la Freeman)*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv];
        }
        break;
      case BETWENDPTS:   /*Betweenness including endpoints*/
        bet[i]+=npred[i]-1.0;
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv]+1.0;
        }
        break;
      case BETPROXIMALSRC:   /*Proximal source betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            if((int)(v->val)!=i)
              bet[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv];
          }
        }
        break;
      case BETPROXIMALTAR:   /*Proximal target betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
            if((int)(v->val)==i)
              bet[wv]+=delta[wv];
          }
        }
        break;
      case BETPROXIMALSUM:  /*Total proximal betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
            if((int)(v->val)!=i)
              bet[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv];
            else
              bet[wv]+=delta[wv];
          }
        }
        break;
      case BETLENSCALED:   /*Length-scaled betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0/gd[wv]+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv];
        }
        break;
      case BETLINSCALED:   /*Linearly-scaled betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0/gd[wv]+delta[wv]);
          if(i!=wv)
            bet[wv]+=gd[wv]*delta[wv];
        }
        break;
      case BETSTRESS:   /*Shimbel's stress centrality (not betweenness!)*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=1.0+delta[wv];
          if(i!=wv)
            bet[wv]+=sigma[wv]*delta[wv];
        }
        break;
      case BETLOAD:    /*Goh's load centrality (must be given transpose graph)*/
        for(j=0;j<n;j++)
          delta[j]=1.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=delta[wv]/(double)npred[wv];
          bet[wv]+=delta[wv];
        }
        break;
    }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
  return sbet;
}



