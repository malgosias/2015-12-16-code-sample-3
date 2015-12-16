## About

Create XSPEC/Sherpa style table model based on a grid of precomputed
spectra read from an input file.

```
$ import tabmod
$ tabmod.table_model(modelname, userparfile, specfile, outfile, clobber)
```

## Parameters

* ```modelname``` - name of the output table model.

* ```specfile``` - file containing grid of precomputed models in ph/cm2/s.

Format of the ```specfile```: 1st column is energy grid, consecutive columns contain grid models.

* ```userparfile``` - file containing user definitions of model parameters.

Format of the ```userparfile``` (see ```userparams_test.txt``` for an example):

```
addmodel True
method 0. 0. ...
name parameter1 parameter2 ...
bottom 0.1 2.0 ...
minimum 0.1 2.0 ...
top 10. 500. ...
maximum 10. 500. ... 
initial 5. 250. ...
delta 0.02 5. ...
numbvals 100 20 ...
value <values of parameter1 separated with spaces> <values of parameter2 separated
with spaces> ...
```

where

```addmodel:``` True for additive model, False for multiplicative model.

```method:``` Interpolation method: 0 for linear; 1 for logarithmic.

```name:``` A sequence of M parameter names, where M is the total number of model parameters.

```top/maximum/bottom/minimum:``` Corresponding boundaries of the model parameters.

```initial:```  Initial Sherpa/XSPEC values of the model parameters.

```delta:``` Initial step in model fitting.

```numbvals:``` A sequence of M numbers, each number is the total no. of grid values of the given model parameter.

```value:``` M sequencies, each sequence contains all the grid values of the given model parameter.

**Note:** The keywords order may be random, but the ```userparfile``` file **must** contain the above 11 entries.

* ```outfile``` - name of the output FITS file.

* ```clobber``` - defaults to False, if True then ```outfile``` will be overwritten.

## Use in Sherpa

```
$ load_table_model('mytabmod', 'mymodel.fits')
$ set_model(mymodel)
```

where ```mymodel.fits``` is the ```outfile``` parameter of ```tabmod.table_model```.

## Use in XSPEC

* Additive model:

```
$ model atable{mymodel.fits}
```

* Multiplicative model:

```
$ model mtable{mymodel.fits}
```
