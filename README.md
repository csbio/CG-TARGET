#**CG-TARGET**

CG-TARGET is a collection of R scripts to be used for the purpose of predicting gene-set targets from chemical-genetic interaction profiles.

##Features

##Installation

###Requirements

This software is written in R, and thus requires a working R installation. Microsoft Open R (or any other R with an optimized BLAS library, such as OpenBLAS) is recommended, as speed at which some steps finish depends highly on the speed of the matrix multiplications involved.

__**The following libraries are required:**__

	data.table
	digest
	ggplot2
	grid
	gridExtra
	optparse
	reshape2
	tools
	yaml

__**The following libraries are optional:**__

(they are only required for the `visualize_gene_targets.r` script, which is optional)

	ctc
	fastcluster

###Downloading CG-TARGET

####Basic

Head on over to https://github.com/csbio/CG-TARGET/releases/ and download the latest release. Extract the compressed folder to a good location from which to run the software (i.e. get it out of your downloads folder!)

####Advanced

If you know what you are doing and want to keep up-to-date with the latest version, clone the repository (git clone https://github.com/csbio/CG-TARGET.git or windows equivalent).


###Setting up environment variables

Required: **TARGET_PATH**

Set the value of this environment variable to the path of the CG-TARGET folder you downloaded and extracted. The scripts from BEAN-counter will look for this variable in your environment, so it must be set!

Optional, but strongly recommended: adding `$TARGET_PATH/scripts/` to your PATH

Adding the scripts folder inside of the CG-TARGET folder to your **PATH** environment variable allows you to execute the commands in the pipeline by calling them only by their names (all demos will assume this has been done).

####How do I set my environment variables?

#####Linux/Mac
The best way to do this is by adding code to the scripts that run every time you open a new shell. If you use the bash shell, then add the following line to either your ~/.bashrc or ~/.bash_profile files:

```
export TARGET_PATH=/your/path/to/CG-TARGET/
```

If you use the c shell (csh), then add the following line to your ~/.cshrc file:

```
setenv TARGET_PATH /your/path/to/CG-TARGET/
```

(but please, **please** do us all a favor and ditch the c shell already)

To append a directory to your PATH variable, add this line to your ~/.bashrc or ~/.bash_profile (or equivalent for ~/.cshrc):

```
export PATH=$PATH:/your/path/to/CG-TARGET/scripts
```

#####Windows

[Tutorial for changing path variable](http://www.computerhope.com/issues/ch000549.htm)

#####Environments

I am an advocate of using virtual environments to manage environment variables in ways that keep one's default environment clean from the many different variables that may be required for different software packages. I currently use anaconda, a combined virtual environment and package manager, to manage my Microsoft R Open and GNU R installations. Venv is another good virtual environment manager.

##How to use

After downloading the software and setting up your environment variables, it is quite simple to run the commands. First, you will need to set up your configuration file (here called `config_file.yaml`). The `example_config_file.yaml` is provided as a template in the main folder.

####View available genetic interaction datasets

`gi_datasets.r`

####View available gene sets and their associated genetic interaction datasets

To view only the gene sets associated with a particular genetic interaction dataset, give the name of that genetic interaction dataset as an argument.

```
gene_sets.r
gene_sets.r SGA_DMA_array_no-supps_sig-only
```

####Generate resampled profiles

`gen_randomized_profiles.r config_file.yaml`

####Predict gene targets for both the real and resampled profiles

```
predict_gene_targets.r config_file.yaml
predict_gene_targets.r --rand config_file.yaml
```

####**Optional:** Create a clustered heat map of gene target prediction scores

This will be a *.CDT file, viewable with Java TreeView. You can choose to view the predictions for the resampled profiles too.

```
visualize_gene_targets.r config_file.yaml
visualize_gene_targets.r --rand config_file.yaml
```

####Predict the gene-set targets

```predict_gene-set_targets.r config_file.yaml```

####Estimate false discovery rate and export final gene-set target prediction tables.

```gene_set_FDR_analysis.r config_file.yaml```


##License

Use of this code is free for academic use and for a 60-day commercial trial period. Sustained commercial use requires the purchase of a license.

In any case, please head on over to z.umn.edu/cgtarget to register for a license. This helps me keep track of the user base :)
