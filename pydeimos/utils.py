"""Some utility functions I often use for development and tests

"""

from __future__ import division, print_function
from future_builtins import *

import os
import cPickle as pickle
import gzip
import astropy.io.fits

import logging
logger = logging.getLogger(__name__)


def read_fits(filepath):
    """Reads simple 1-hdu FITS file into a numpy arrays

    A transposition is done, so that the indexes [x,y] of the numpy array follow the orientation of x and y in DS9
    and SExtractor.

    Parameters
    ----------
    filepath : str
        Filepath to read the array from

    """
    a = astropy.io.fits.getdata(filepath).transpose()
    logger.info("Read FITS images %s from file %s" % (a.shape, filepath))
    return a


def write_fits(a, filepath, overwrite=True):
    """Writes a simple 2D numpy array into a FITS file

    As for read_fits, a transposition is applied to conserve the orientation.

    Parameters
    ----------
    a         : array
    filepath  : str
    overwrite : bool
        Set this to False if an existing file should not be overwritten

    """
    if os.path.exists(filepath) and overwrite:
        logger.info("File %s exists, I will overwrite it!" % (filepath))

    astropy.io.fits.writeto(filepath, a.transpose(), overwrite=overwrite)
    logger.info("Wrote %s array into %s" % (a.shape, filepath))



def write_pickle(obj, filepath, protocol = -1):
    """Write a python object into a pickle file

    Parameters
    ----------
    obj : object
        A python object to be written into the pickle file
    filepath : str
        Filepath to write to.
        If filepath ends with .gz, uses gzip to compress the pickle.
    protocol : int
        Protocol to use. Leave this to -1 : it will use the latest binary protocol of pickle.
    """
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath, 'wb')
    else:
        pkl_file = open(filepath, 'wb')

    pickle.dump(obj, pkl_file, protocol)
    pkl_file.close()
    logger.info("Wrote '{}'".format(filepath))


def read_pickle(filepath):
    """Reads a pickle file and returns the object it contains.

    Parameters
    ----------
    filepath : str
        Filepath to read from.
        If filepath ends with .gz, uses gzip to uncompress the pickle.

    """
    if os.path.splitext(filepath)[1] == ".gz":
        pkl_file = gzip.open(filepath,'rb')
    else:
        pkl_file = open(filepath, 'rb')
    obj = pickle.load(pkl_file)
    pkl_file.close()
    logger.info("Read '{}'".format(filepath))
    return obj
