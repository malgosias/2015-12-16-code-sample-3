# Create XSPEC/Sherpa style table models
#
# Used in:
#
# You et al. (2015), ApJ submitted, arXiv:1506.03959
# Adhikari et al. (2016), in preparation
#
# table_model
# Create Xspec/Sherpa style table model
#
# read_userparams
# Read model parameters from a file provided by the user
#
# get_paramval
# Find combinations of model paramaters corresponding to the input
# spectra
#
# read_input_spectra
# Read grid of model spectra from a file prvided by the user


from astropy.io import fits
import numpy as np
import os

def table_model(modelname, userparfile, specfile, outfile,
                clobber=False):
    
# Create Xspec/Sherpa style table model
#
# :rtype: None, <outfile> fits file created.
#
# :param modelname: name of the table model displayed in Xspec/Sherpa
# :param userparfile: file with user keywords for the table model
# :specfile: file with grid of energy spectra, 1st column: energy grid
#            in keV; consecutive columns: model spectra
# :param outfile: name of the output fits file
# :param clobber: T/F, if T outfile will be overwritten
#
# Tmp file 'tmp_tabmod.fits' created and removed.

# ----------- READ IN USER INPUT --------------------

    if (os.path.isfile(outfile) and not clobber):
        raise NameError('Output file ' + outfile +
                        ' exists and clobber set to false\n')

    userdict = read_userparams(userparfile)
  
    # Convert userdict['value'] into a list of tuples; one tuple
    # for each model parameter

    # value_not_padded: used in get_paramval to calculate array with
    # combinations of model parameters
    idx = 0
    value_not_padded = []
    for item in userdict['numbvals']:
        value_not_padded.append(tuple(userdict['value'][idx:idx+item]))
        idx += item

    # value_padded: format required for col10 of the fits file
    value_padded = []
    maxnum = max(userdict['numbvals'])
    for val in value_not_padded:
      if (len(val) != maxnum):
        n = maxnum - len(val)
        value_padded.append(val + (0.,)*n)
      else:
        value_padded.append(val)

    # model energy grid (bin edges!)
    energy = np.loadtxt(specfile)[:,0] # energy in keV
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    paramval = get_paramval(value_not_padded)

    # list of tuples with spectra
    input_spectra = read_input_spectra(specfile, userdict['numbvals'])

# ----------- END ---------------------------------

    # initialize fits file by creating user parameters extension

    col1 = fits.Column(name='NAME', format='12A',
                       array=userdict['name'])
    col2 = fits.Column(name='METHOD', format='J',
                       array=userdict['method'])
    col3 = fits.Column(name='INITIAL', format='E',
                       array=userdict['initial'])
    col4 = fits.Column(name='DELTA', format='E',
                       array=userdict['delta'])
    col5 = fits.Column(name='MINIMUM', format='E',
                       array=userdict['minimum'])
    col6 = fits.Column(name='BOTTOM', format='E',
                       array=userdict['bottom'])
    col7 = fits.Column(name='TOP', format='E',
                       array=userdict['top'])
    col8 = fits.Column(name='MAXIMUM', format='E',
                       array=userdict['maximum'])
    col9 = fits.Column(name='NUMBVALS', format='J',
                       array=userdict['numbvals'])
    col10 = fits.Column(name='VALUE',
                        format=np.str(np.max(userdict['numbvals']))+'E',
                        array=value_padded)

    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8,
                         col9, col10])

    tbhdu = fits.BinTableHDU.from_columns(cols)

    # update header of user parameters extension

    nintparm = len(userdict['numbvals'])
    
    tbhdr = tbhdu.header

    tbhdr.set('EXTNAME', 'PARAMETERS',
              'name of this binary table extension')
    tbhdr.set('HDUCLASS', 'OGIP',
              'format conforms to OGIP standard')
    tbhdr.set('HDUCLAS1', 'XSPEC TABLE MODEL',
              'model spectra for XSPEC')
    tbhdr.set('HDUCLAS2', 'PARAMETERS',
              'extension containing parameter info')
    tbhdr.set('HDUVERS1', '1.0.0',
              'version of format')
    tbhdr.set('NINTPARM', nintparm,
              'Number of interpolation parameters')
    tbhdr.set('NADDPARM', 0,
              'Number of additional parameters')

    if (os.path.isfile('tmp_tabmod.fits')):
        os.remove('tmp_tabmod.fits')

    tbhdu.writeto('tmp_tabmod.fits')

    # update primary header

    hdulist = fits.open('tmp_tabmod.fits')

    prihdr = hdulist[0].header
    prihdr['bitpix'] = 16
    prihdr.set('modlname', modelname, 'model name')
    prihdr.set('modlunit', 'photons/cm^2/s', 'model units')
    prihdr.set('redshift', True,
               'If true then redshift will be included as a par')
    prihdr.set('addmodel', userdict['addmodel'],
               'If true then this is an additive table model')
    prihdr.set('hduclass', 'OGIP',
               'format conforms to OGIP standard')
    prihdr.set('hduclas1', 'XSPEC TABLE MODEL',
               'model spectra for XSPEC')
    prihdr.set('hduvers1', '1.0.0', 'version of format')

    if (os.path.isfile(outfile)):
        os.remove(outfile)

    hdulist.writeto(outfile)
 
    hdulist.close()
    os.remove('tmp_tabmod.fits')
  
    # append extension energies and update its header

    col1 = fits.Column(name='ENERG_LO', format='E', array=energ_lo,
                       unit='keV')
    col2 = fits.Column(name='ENERG_HI', format='E', array=energ_hi,
                       unit='keV')

    cols = fits.ColDefs([col1, col2])

    tbhdu_energies = fits.BinTableHDU.from_columns(cols)

    hdr = tbhdu_energies.header
    hdr.set('EXTNAME', 'ENERGIES',
            'name of this binary table extension')
    hdr.set('HDUCLASS', 'OGIP',
            'format conforms to OGIP standard')
    hdr.set('HDUCLAS1', 'XSPEC TABLE MODEL',
            'model spectra for XSPEC')
    hdr.set('HDUCLAS2', 'ENERGIES',
            'extension containing energy bins info')
    hdr.set('HDUVERS1', '1.0.0', 'version of format')

    fits.append(outfile, tbhdu_energies.data, hdr)

    # append extension spectra and update its header

    col1 = fits.Column(name='PARAMVAL',
                       format=np.str(nintparm)+'E', array=paramval)
    col2 = fits.Column(name='INTPSPEC',
                       format=np.str(len(energ_lo))+'E',
                       array=input_spectra, unit='photons/cm^2/s')

    cols = fits.ColDefs([col1, col2])

    tbhdu_spectra = fits.BinTableHDU.from_columns(cols)

    hdr = tbhdu_spectra.header
    hdr.set('EXTNAME', 'SPECTRA',
            'name of this binary table extension')
    hdr.set('HDUCLASS', 'OGIP',
            'format conforms to OGIP standard')
    hdr.set('HDUCLAS1', 'XSPEC TABLE MODEL',
            'model spectra for XSPEC')
    hdr.set('HDUCLAS2', 'MODEL SPECTRA',
            'extension containing model spectra')
    hdr.set('HDUVERS1', '1.0.0', 'version of format')

    fits.append(outfile, tbhdu_spectra.data, hdr)
  
    return


def read_userparams(fname):
# Read model parameters from the <fname> file
#
# :rtype: dictionary, required keyword : value
# :param fname: name of the file containing required keywords and values

    # Build a dictonary with user input
    
    with open(fname, 'r') as f:
        # ignore empty lines
        lines = [ l for l in f.readlines() if l.strip() ]

    expected_no_lines = 11

    if (len(lines)!= expected_no_lines):
        raise NameError(fname + ': expected ' +
                        np.str(expected_no_lines) + ' lines, but ' +
                        np.str(len(lines)) + ' lines found\n')

    lines_split = []
    for l in lines:
        separate = l.rstrip().split()
        if (len(separate)<2):
            raise NameError(fname + ': cannot split line, ' + l + '\n')
        lines_split.append(separate)

    # Recover items that will be the keys and check if all provided

    keys_required = ['addmodel', 'name', 'method', 'initial',
                     'delta', 'minimum', 'bottom', 'top',
                     'maximum', 'numbvals', 'value']
    keys = []
    for item in lines_split:
        keys.append(item[0].lower())

    if not (set(keys) == set(keys_required)):
        raise NameError(fname + ': missing entry, required keys: ',
                        keys_required)

    userdict = {}
    for item in lines_split:
        if ('numbval' in item[0].lower()):
            userdict[item[0].lower()] = [ int(float(i))
                                         for i in item[1:] ]
        elif ('addmodel' in item[0].lower()):
            userdict[item[0].lower()] = bool(item[1])
        elif ('name' in item[0].lower()) or ('unit' in item[0].lower()):
            userdict[item[0].lower()] = item[1:]
        else:
            userdict[item[0].lower()] = [ float(i) for i in item[1:] ]

    if ( len(userdict['value']) != np.sum(userdict['numbvals']) ):
        raise NameError(fname +
                        ': numbvals and value entries do not match\n')

    return userdict


def get_paramval(arrays):
# Find combinations of model paramaters corresponding to the input
# spectra, needed for extension spectra
#
# :rtype: list of tuples, each tuple contains a combination of model
#         parameters corresponding to one of the input grid spectrum.
# :param arrays: list of tuples, each tuple contains grid values of
#                model parameters, not padded with 0
#
    arrays = [np.asarray(a) for a in arrays]
    shape = (len(x) for x in arrays)

    tmp = np.indices(shape, dtype=int)
    ix = tmp.reshape(len(arrays), -1).T
    out = tmp.reshape(len(arrays), -1).T.astype('float')
    
    for n, arr in enumerate(arrays):
        out[:, n] = arrays[n][ix[:, n]]

    param_combinations = []
    for item in out:
        param_combinations.append(tuple(item))
  
    return param_combinations


def read_input_spectra(fname, numbvals):
# Read grid of model spectra from the <fname> file
#
# :rtype: list of tuples containing model spectra.
#
# :param fname: name of the file containing the grid of model spectra.
# :param numbvals: list with total grid number of model parameters.
    
    d = np.loadtxt(fname)
    parameter_combinations = np.prod(numbvals)
    numspectra = len(d.T)-1
    # spectra in columns, 1st column is energy grid, skip it
    if ( numspectra != parameter_combinations ):
        raise NameError(fname + ': No. of spectra ' +
                        np.str(numspectra) +
                        ', different from declared param combinations '
                        + np.str(parameter_combinations) + '\n')
    input_spectra = []
    for i in range( len(d.T)-1 ):
        # d[:-1] - match no. of energy bins;
        # T[1:] - skip the 1st column with energy vector
        input_spectra.append(tuple(d[:-1].T[1:][i]))

    return input_spectra



# -------------  miscellaneous -----------

'''
    
def generate_spectra():
# Tool to generate a grid of test spectra and test parameter file

    energy = np.linspace(0.5, 10., 30)   # adjust to generate
    par1 = np.linspace(1., 3., 5)       # parameter 1
    par2 = np.linspace(2., 50., 9)      # parameter 2

    # Uncomment the following 7 lines, enter your model expression
    #spec = [energy]
    #for p1 in par1:
    #    for p2 in par2:
    #        tmp = < model expression >
    #        spec.append(tmp)
    #
    #np.savetxt('specfile_test.txt', np.array(spec).T, fmt='%10.4f' )

    f = open('userparams_test.txt', 'w+')

    f.write('addmodel True\n')
    f.write('name par1 par2\n')
    f.write('method 0. 0.\n')
    f.write('top ' + np.str(par1[-1]) + ' ' + np.str(par2[-1]) + '\n')
    f.write('maximum ' + np.str(par1[-1]) + ' ' + np.str(par2[-1]) +
                                                                '\n')
    f.write('bottom ' + np.str(par1[0]) + ' ' + np.str(par2[0]) + '\n')
    f.write('minimum ' + np.str(par1[0]) + ' ' + np.str(par2[0]) + '\n')
    f.write('initial ' + np.str(par1.mean()) + ' ' +
                                            np.str(par2.mean()) + '\n')
    # update the definition of delta is required
    f.write('delta ' + np.str( (par1[1:] - par1[:-1]).min()/5. ) + ' ' +
                       np.str( (par2[1:] - par2[:-1]).min()/5. ) + '\n')
    f.write('numbvals ' + np.str(len(par1)) + ' ' + np.str(len(par2)) +
                                                                '\n')

    par1vals = ''
    for p1 in par1:
        par1vals += np.str(p1) + ' '

    par2vals = ''
    for p2 in par2:
        par2vals += np.str(p2) + ' '

    f.write('value ' + par1vals + par2vals + '\n')

    f.close()

    return

'''


