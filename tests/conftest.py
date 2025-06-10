import pytest
import os
import json

def load_test_data(file_name,inner_path):
    '''
      loads test data directly from json-file on the fly
    '''
    folder = os.path.abspath(os.path.dirname(__file__))
    jsonfile = os.path.join(folder, file_name)
    with open(jsonfile) as file:
        content = json.load(file)
        tokens = inner_path.split(',')
        for token in tokens:
            content = content[token]
        if 'bessel' in tokens:
            print('uses bessel test data')
            data = [ (item['nu'],item['z'][0]+item['z'][1]*1j,item['res'][0]+item['res'][1]*1j) for item in content ]
            return data
        data = [ item for item in content ]
        return data