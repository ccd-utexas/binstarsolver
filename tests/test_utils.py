"""Tests for binstarsolver/utils.py.

"""

from __future__ import absolute_import, print_function, division
import sys
import pytest
import binstarsolver

# TODO: make test with pytest and redo functions below.

def is_calc_close_to_ref(calc, ref, tol=0.001, name=None, verbose=False):
    """Test that a calulated value is within a tolerance of the reference value.
    
    Parameters
    ----------
    calc : float
        Calculated value to be tested.
    ref : float
        Reference value to be tested against.
    tol : {0.001}, float, optional
        Tolerance value to determine if `calc` is close to `ref`.
    name : {None}, string, optional
        Name of quantity being tested for printed output.
    verbose : {False, True}, bool, optional
        Print outcome of test if `True`.
    
    Returns
    -------
    is_close : bool
        Boolean for whether or not calculated value was
        within the tolerance of the reference value.
    
    """
    is_close = None
    abs_diff = abs(ref - calc)
    if (abs_diff < tol) or (abs_diff == tol):
        is_close = True
        if verbose:
            print(("INFO: Reference and calculated values match within tolerance.\n" +
                   "name          = {name}\n" +
                   "ref           = {ref}\n" +
                   "calc          = {calc}\n" +
                   "abs(ref-calc) = {arc}\n" +
                   "tol           = {tol}").format(name=name, ref=ref, calc=calc, arc=abs_diff, tol=tol))        
    else:
        is_close = False
        if verbose:
            print(("ERROR: Reference and calculated values do not match within tolerance.\n" +
                   "name          = {name}\n" +
                   "ref           = {ref}\n" +
                   "calc          = {calc}\n" +
                   "abs(ref-calc) = {arc}\n" +
                   "tol           = {tol}").format(name=name, ref=ref, calc=calc, arc=abs_diff, tol=tol),
                  file=sys.stderr)
    if is_close is None:
        raise AssertionError("Program error. `is_close` was not set to `True` or `False`.")
    return is_close


def summarize_failures(tests_failed):
    """Print summary of tests that have failed.
    
    Parameters
    ----------
    tests_failed : array
        Array of string names of tests that have failed.
        
    Returns
    -------
    None
    
    Notes
    -----
    The names from `tests_failed` will be formatted with a message and printed to stdout.
    If `tests_failed` == [], then the message will state that all tests passed.
    
    """
    if tests_failed == []:
        print("INFO: Tests complete. All tests passed.")
    else:
        print(("ERROR: Tests complete. Some tests failed:\n" +
               "tests_failed = {tf}").format(tf=tests_failed),
              file=sys.stderr)
    return None
