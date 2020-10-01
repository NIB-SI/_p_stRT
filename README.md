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
  ## install R and pandoc
  ```
  apt-key adv --keyserver [from this location or server] --recv-keys [retrieve key(s)]
  add-apt-repository ‘deb https://cloud.r-project.org/bin/linux/ubuntu [type appropriate selection from https://cloud.r-project.org/bin/linux/ubuntu/]’

  apt-get update
  apt-get install r-base
  apt-get install r-base-dev
  apt install build-essential
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
  ## command line cheatsheet
  
  <img src="https://www.git-tower.com/blog/media/pages/posts/command-line-cheat-sheet/1073300074-1586415841/command-line-cheat-sheet-large01.png" width="600"/>
  
  linked from [www.git-tower.com](www.git-tower.com)
  
  ## relative path
  
  always use relative path
  
      . represents the current directory

      .. represents the parent directory

   ```ls .. ``` list information about the files one level up, 

   ```ls ../.. ``` or ```ls ../../``` list information about the files two levels up, etc

<br/>

## 📝 Useful links and code lines

## :clipboard: Obtaining materials from the GitHub

- use wget to pull down the ```raw``` file
```
$ wget https://raw.githubusercontent.com/username/reponame/path/to/file
```
- use git clone to pull the complete repository (prerequisites: user with ```sudo``` privileges)
```
$ sudo apt update
$ sudo apt install git
$ git clone https://github.com/username/reponame.git
```

## :clipboard: Install Perl on Ubuntu-like Linux OS
- Install Perl on Ubuntu-like Linux OS (prerequisites: user with ```sudo``` privileges)
```
$ sudo apt update
$ sudo apt-get install perl
# check version 
$ perl -v
```
More details for installing Perl can be found at [perl.org](www.perl.org/get.html#unix_like)

- Use [CPAN](www.cpan.org) (‘Comprehensive Perl Archive Network’) to install Perl modules

- Install different version of Perl
```
$ sudo cpan App::perlbrew
$ perlbrew init

# see which versions are available:
$ perlbrew available

# install version 5.X.Y
$ perlbrew install perl-5.X.Y

# list all installed versions
$ perlbrew list

# change Perl for the current shell # (or per your sessions)
$ perlbrew use perl-5.X.Y # (or perlbrew switch perl-5.X.Y)
$ which perl

# revert version to default for the current shell # (or per your sessions)
$ perlbrew off # (or perlbrew switch-off)
```

## :clipboard: pip 
Installing pip for Python 3 and Python2 on Ubuntu-like OS (prerequisites: user with sudo privileges)
```
# Installing pip for Python 3
$ sudo apt update
$ sudo apt install python3-pip

# Installing pip for Python 2
$ sudo apt update 
$ sudo apt install python2
$ curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py
$ sudo python2 get-pip.py
```


<br/>
