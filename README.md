# general information 

## ISA-tab
<img src="https://isa-tools.org/wp-content/uploads/2015/12/ISAmodel-structure.png" width="400"/>


linked from <https://isa-tools.org/format/specification.html>

## pISA-tree structure
pISA-tree: Standard project directory tree (ISA-tab compliant) 

<https://github.com/NIB-SI/pISA>

### seekr
pISA-tree: Standard project directory tree (ISA-tab compliant)

<https://github.com/NIB-SI/seekr>

### pisaR
R support for pISA-tree

<https://github.com/NIB-SI/pisar>

## FAIRDOMHub
FAIRDOMHub: a repository and collaboration environment for sharing systems biology research

<https://fairdomhub.org/>

# recommended OS

LTS versions of:

  <https://linuxmint.com/>

  <https://ubuntu.com/>
  
  # R on linux
  ## install
  ```
  apt-get update
  apt-get install r-base
  apt-get install pandoc
  ```
  ## invoke R from the command line
  ```
  R
  ```
  ## quit R
  ```
  q()
  ```
  ## help
  ```
  R --help
  ```
  ## install packages/libraries
  ```
  R -e 'install.packages("rmarkdown", repos="https://cran.rstudio.com/")'
  ```
  ## run app, script, ...
  ```
  R -e "shiny::runApp('./pathToShinyApp/name.R')"
  Rscript -e "rmarkdown::render('./pathToScript/scriptName.Rmd')"
  ```
  
  # Biopython and conda
  ```
  pip install biopython
  ```
  [biopython.org](https://biopython.org/)
  
  [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
  
  # other packages/libraries
  ```
  
  apt-get install dos2unix
  
  apt-get install moreutils
  
  R -e 'install.packages(c("Rcpp", "httpuv"))'

  R -e 'install.packages("shiny")'

  R -e 'install.packages("devtools", repos="https://cran.rstudio.com/")'

  apt-get install libssl-dev libxml2-dev libcurl4-openssl-dev libcurl4-gnutls-dev curl
  ```
  
  # cheatsheet
  <img src="https://www.git-tower.com/blog/media/pages/posts/command-line-cheat-sheet/1073300074-1586415841/command-line-cheat-sheet-large01.png" width="600"/>
  
  command line cheatsheet linked from [www.git-tower.com](www.git-tower.com)
  
  always use relative path
  
. represents the current directory
.. represents the parent directory
ls .. lists one level up
cd ../.. moves two levels up
