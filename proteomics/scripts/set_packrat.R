################################################################################
# CureMILS project:                                                            #
# PROTEOMICS: Set packrat environment for the proteomics analysis              #
################################################################################
packrat::off()

if(!require(packrat, quietly = T)){
        install.packages("packrat",
                         repos = 'http://cran.us.r-project.org')
}
if(!require(ggplot2, quietly = T)){
        install.packages("ggplot2",
                         repos = 'http://cran.us.r-project.org')
}

if(!require(argparser, quietly = T)){
        install.packages("argparser",
                         repos = 'http://cran.us.r-project.org')
}
library(argparser)


# Terminal argument parser
################################################################################
parser <- arg_parser("Set packrat environment given the root directory.")

parser <- add_argument(parser = parser,
                       arg = "input",
                       help = "Root directory of the project",
                       flag = F)

parsed <- parse_args(parser)

# Directory stuff
################################################################################

rootDir <- parsed$input

# Create packrat environment
################################################################################

packrat::init(rootDir)

packrat::on(project = rootDir)

packrat::set_opts(external.packages = c("BiocManager", "GO.db", "org.Hs.eg.db", "readxl", "argparser"))


install.packages("data.table", repos = "https://pbil.univ-lyon1.fr/CRAN/")
install.packages("ndjson")
install.packages("neat",
                 repos = "https://pbil.univ-lyon1.fr/CRAN/")
#BiocManager::install("ropls", update = F)
install.packages("caret",
                 repos= "https://pbil.univ-lyon1.fr/CRAN/")

packrat::snapshot(project = rootDir)