# import os
# import pyteomics
# pyteomics.__path__ = [os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, 'pyteomics'))]
# from itertools import product
# import unittest
# from pyteomics.mzxml import MzXML, read, chain
# from pyteomics import xml
# from data import mzxml_spectra
# import tempfile
# import shutil
# import numpy as np

require 'minitest/autorun'
require_relative '../rbteomics/mzxml'
require_relative 'data'
require 'set'
require 'tempfile'
require 'numpy'

class TestMzxml < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'test.mzXML'
  end

  def testReadSpectrum

  end

  def test_decoding
    MzXML.new(@path, decode_binary: true) do |reader|
      
    end
  end
end
