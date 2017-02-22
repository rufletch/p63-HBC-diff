#!/bin/bash

R --vanilla <<<"rmarkdown::render('oeHBCdiff_clusterPlots.Rmd')"
R --vanilla <<<"rmarkdown::render('oeHBCdiff_devorderplots.Rmd')"
R --vanilla <<<"rmarkdown::render('oeHBCdiff_cellCycle.Rmd')"
R --vanilla < oeHBCdiff_tf_hm.R 
R --vanilla <<<"rmarkdown::render('oeHBCdiff_tf.Rmd')"
R --vanilla <<<"rmarkdown::render('oeHBCdiff_geneClustHeatmaps.Rmd')"
R --vanilla <<<"rmarkdown::render('oeHBCdiff_genePlots.Rmd')"
R --vanilla <<<"rmarkdown::render('oeHBCdiff_GSEAplots.Rmd')"
R --vanilla < oeHBCdiff_volcano.R 
R --vanilla < oeHBCdiff_OR.R 