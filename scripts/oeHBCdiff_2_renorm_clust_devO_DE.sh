#!/bin/bash

ncores=23
nrm="none_fq_qc1_nobio_nobatch"

R --vanilla < oeHBCdiff_makeSE.R --args --expt oeHBCdiff --nrm none,fq,qc_k=1,no_bio,no_batch --ncores $ncores
R --vanilla < oeHBCdiff_exclude.R --args --expt oeHBCdiff --ncores $ncores --nrm $nrm
R --vanilla < oeHBCdiff_filtering.R --args --expt oeHBCdiff --norm $nrm
R --vanilla < oeHBCdiff_norm.R --args --expt oeHBCdiff --ncores $ncores

R --vanilla < oeHBCdiff_makeSE.R --args --expt oeHBCdiff --nrm none,fq,qc_k=1,no_bio,no_batch --ncores $ncores

R --vanilla < oeHBCdiff_clust.R --args --expt oeHBCdiff --nrm $nrm --ncores $ncores



R --vanilla <<<"rmarkdown::render('oeHBCdiff_slingshot.Rmd')"