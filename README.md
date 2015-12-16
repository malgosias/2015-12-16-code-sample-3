## About

Create XSPEC/Sherpa style table model based on a grid of precomputed
spectra read from an input file.

```
$ import tabmod
$ tabmod.table_model(modelname, userparfile, specfile, outfile, clobber)
```

## Parameters

* ```modelname``` - name of the user table model.

* ```specfile``` - file containing grid of precomputed models in ph/cm2/s.

Format of the ```specfile```: 1st column is energy grid; consecutive columns contain grid models.

* ```userparfile``` - file containing user definitions of model parameters.

Format of the ```userparfile``` (the order may be random, but must containt these 11 entries):

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

```
addmodel - True for additive model, False for multiplicative model
method - interpolation method: 0 for linear; 1 for logarithmic
top/maximum/bottom/minimum - corresponding boundaries of the model parameters
initial - initial fit values of the model parameters
delta - initial step in model fitting
numbvals - <no. of values for par1> <no. of values for par2> ...
```

* ```outfile``` - name of the output FITS file.

* ```clobber``` - defaults to False, if True then ```outfile``` will be overwritten.

## Use in Sherpa:

```
$ load_table_model('mytabmod', 'mymodel.fits')
$ set_model(mymodel)
```

where ```mymodel.fits``` is the ```outfile``` parameter of ```tabmod.table_model```.

## Use in XSPEC:

* Additive model:

```
$ model atable{mymodel.fits}
```

* Multiplicative model:

```
$ model mtable{mymodel.fits}
```
