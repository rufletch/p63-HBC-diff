#!/bin/bash

ncores=23

R --vanilla < oeHBCdiff_filtering.R --args --expt oeHBCdiff  
R --vanilla < oeHBCdiff_norm.R --args --expt oeHBCdiff --ncores $ncores



