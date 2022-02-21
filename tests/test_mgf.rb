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
    # for func in [mgf.read, mgf.MGF, mgf.IndexedMGF]:
    #   # http://stackoverflow.com/q/14246983/1258041
    #   self.assertEqual(data.mgf_spectra_long, list(func(self.path)))
    #   self.assertEqual(data.mgf_spectra_short, list(func(self.path, False)))
    #   with func(self.path) as reader:
    #       self.assertEqual(data.mgf_spectra_long, list(reader))
    #   with func(self.path, False) as reader:
    #       self.assertEqual(data.mgf_spectra_short, list(reader))

    assert_equal Mgf_spectra_long, Mgf::Read.call(@path)
      # self.assertEqual(data.mgf_spectra_short, list(func(self.path, False)))
      # with func(self.path) as reader:
      #     self.assertEqual(data.mgf_spectra_long, list(reader))
      # with func(self.path, False) as reader:
      #     self.assertEqual(data.mgf_spectra_short, list(reader))


  end
end