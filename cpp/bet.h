/*
######################################################################
#
# nli.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 3/29/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for nli.c.
#
######################################################################
*/
#ifndef NLI_H
#define NLI_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "geodist.h"

/*Definitions for various measures to be computed by betweenness_R (mostly
based on Brandes (2008)); note that some are not forms of betweenness, but
can be calculated using that routine.*/
#define BETSTANDARD     0         /*"Standard" form betweenness (a la Freeman)*/
#define BETWENDPTS      1                    /*Betweenness including endpoints*/
#define BETPROXIMALSRC  2                        /*Proximal source betweenness*/
#define BETPROXIMALTAR  3                        /*Proximal target betweenness*/
#define BETPROXIMALSUM  4                         /*Total proximal betweenness*/
#define BETLENSCALED    5                          /*Length-scaled betweenness*/
#define BETLINSCALED    6                        /*Linearly-scaled betweenness*/
#define BETSTRESS       7     /*Shimbel's stress centrality (not betweenness!)*/
#define BETLOAD       8/*Goh's load centrality (must be given transpose graph)*/
      
      
/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

SEXP betweenness_partial(SEXP mat, SEXP sn, SEXP sm, SEXP smeasure, SEXP st, SEXP end);



#endif
