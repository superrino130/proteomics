require 'minitest/autorun'
require_relative 'data'

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
require_relative '../rbteomics/mgf'
require_relative '../rbteomics/auxiliary/file_helpers'
require_relative '../rbteomics/auxiliary/structures'
require_relative '../rbteomics/auxiliary/utils'
require 'set'


class MGFTest < Minitest::Test
  def setup
    @maxDiff = nil
    @_encoding = 'utf-8'
    @path = 'tests/test.mgf'
    @header = Mgf.read_header(@path)
    f = Mgf::Read.call(@path)
    # @spectra = f.to_a
  end

  def test_read
    # assert_equal Mgf_spectra_long, Mgf::Read.call(@path)
    # assert_equal Mgf_spectra_long, Mgf::MGF.new(@path)
    assert_equal Mgf_spectra_long, Mgf::IndexedMGF.new(@path)


  end
end