#import re,sys,os, glob
#from string import *
#import math
import numpy
#import random
#from numpy import array
import scipy
#from scipy import stats
#import csv
#import pvalue_combine
#from sets import Set
#import math

#from scipy.stats import *

def combine_pvalues(pvalues, method='fisher', weights=None):
    """
    Methods for combining the p-values of independent tests bearing upon the
    same hypothesis.
    Parameters
    ----------
    pvalues : array_like, 1-D
        Array of p-values assumed to come from independent tests.
    method : {'fisher', 'stouffer'}, optional
        Name of method to use to combine p-values. The following methods are
        available:
        - "fisher": Fisher's method (Fisher's combined probability test),
          the default.
        - "stouffer": Stouffer's Z-score method.
    weights : array_like, 1-D, optional
        Optional array of weights used only for Stouffer's Z-score method.
    Returns
    -------
    statistic: float
        The statistic calculated by the specified method:
        - "fisher": The chi-squared statistic
        - "stouffer": The Z-score
    pval: float
        The combined p-value.
    Notes
    -----
    Fisher's method (also known as Fisher's combined probability test) [1]_ uses
    a chi-squared statistic to compute a combined p-value. The closely related
    Stouffer's Z-score method [2]_ uses Z-scores rather than p-values. The
    advantage of Stouffer's method is that it is straightforward to introduce
    weights, which can make Stouffer's method more powerful than Fisher's
    method when the p-values are from studies of different size [3]_ [4]_.
    Fisher's method may be extended to combine p-values from dependent tests
    [5]_. Extensions such as Brown's method and Kost's method are not currently
    implemented.
    .. versionadded:: 0.15.0
    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Fisher%27s_method
    .. [2] http://en.wikipedia.org/wiki/Fisher's_method#Relation_to_Stouffer.27s_Z-score_method
    .. [3] Whitlock, M. C. "Combining probability from independent tests: the
           weighted Z-method is superior to Fisher's approach." Journal of
           Evolutionary Biology 18, no. 5 (2005): 1368-1373.
    .. [4] Zaykin, Dmitri V. "Optimally weighted Z-test is a powerful method
           for combining probabilities in meta-analysis." Journal of
           Evolutionary Biology 24, no. 8 (2011): 1836-1841.
    .. [5] https://en.wikipedia.org/wiki/Extensions_of_Fisher%27s_method
    """
    pvalues = numpy.asarray(pvalues)
    if pvalues.ndim != 1:
        raise ValueError("pvalues is not 1-D")

    if method == 'fisher':
        Xsq = -2 * numpy.sum(numpy.log(pvalues))
        pval = scipy.stats.distributions.chi2.sf(Xsq, 2 * len(pvalues))
        return (Xsq, pval)
    elif method == 'stouffer':
        if weights is None:
            weights = numpy.ones_like(pvalues)
        elif len(weights) != len(pvalues):
            raise ValueError("pvalues and weights must be of the same size.")

        weights = numpy.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("weights is not 1-D")

        Zi = scipy.stats.distributions.norm.isf(pvalues)
        Z = numpy.dot(weights, Zi) / numpy.linalg.norm(weights)
        pval = scipy.stats.distributions.norm.sf(Z)

        return (Z, pval)
    else:
        raise ValueError(
            "Invalid method '%s'. Options are 'fisher' or 'stouffer'", method)


if __name__ == "__main__":
    pvalues = [0.1,0.01,0.02]
    print(combine_pvalues(pvalues, method='fisher'))
    #print combine_pvalues( pvalues, method='stouffer' )
