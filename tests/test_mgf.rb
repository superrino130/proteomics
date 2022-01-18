require 'minitest/autorun'

# import os
# import numpy as np
require 'numpy'
np = Numpy
# import pyteomics

# pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]

# import tempfile
# import unittest
# import pickle
# import shutil
# import json
# from collections import OrderedDict
# from pyteomics import mgf, auxiliary as aux
# import data

class MGFTest < Minitest::Test
  def setup
    @maxDiff = nil
    @_encoding = 'utf-8'
    @path = 'test.mgf'
    @header = read_header(@path)
    
  end
end