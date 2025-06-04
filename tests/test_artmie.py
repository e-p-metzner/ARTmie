import pytest
import os
import numpy as np
import json
import ARTmie


def load_test_data(file_name,inner_path):
    '''
      loads test data directly from json-file on the fly
    '''
    folder = os.path.abspath(os.path.dirname(__file__))
    jsonfile = os.path.join(folder, file_name)
    with open(jsonfile) as file:
        content = json.load(file)
        for token in inner_path.split(','):
            content = content[token]
        if 'bessel' in inner_path.split(','):
            print('uses bessel test data')
            data = [ (item['nu'],item['z'][0]+item['z'][1]*1j,item['res'][0]+item['res'][1]*1j) for item in content ]
            return data
        return content


@pytest.mark.parametrize("nu,z,expectation", load_test_data('testdata_bessel.json','bessel,J'))
def test_besselj(nu,z,expectation):
    if np.abs(expectation) < 0.0001:
        assert np.abs(ARTmie.besselj(nu,z+0j)-expectation)<1.0e-8
    else:
        assert np.abs(ARTmie.besselj(nu,z+0j)/expectation-1.0)<1.0e-8


@pytest.mark.parametrize("nu,z,expectation", load_test_data('testdata_bessel.json','bessel,Y'))
def test_bessely(nu,z,expectation):
    if np.abs(expectation) < 0.0001:
        assert np.abs(ARTmie.bessely(nu,z+0j)-expectation)<1.0e-8
    else:
        assert np.abs(ARTmie.bessely(nu,z+0j)/expectation-1.0)<1.0e-8


@pytest.mark.parametrize("nu,z,expectation", load_test_data('testdata_bessel.json','bessel,H'))
def test_hankel(nu,z,expectation):
    if np.abs(expectation) < 0.0001:
        assert np.abs(ARTmie.hankel(nu,z+0j,1)-expectation)<1.0e-8
    else:
        assert np.abs(ARTmie.hankel(nu,z+0j,1)/expectation-1.0)<1.0e-8

