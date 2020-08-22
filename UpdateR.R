if (!require('devtools')) install.packages('devtools'); require('devtools')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr')
#install_Rtools()
devtools::install_github('talgalili/installr')

require(installr)
updateR()


old.packages()
update.packages()
install.packages("packagename")

sessionInfo()
browseVignettes("ggplot2")