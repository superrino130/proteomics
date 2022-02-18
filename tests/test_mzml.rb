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
require_relative '../rbteomics/mzml'
require_relative '../rbteomics/xml'
require_relative 'data'
require 'set'
require 'tempfile'
require 'numpy'

class TestMzml < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'tests/test.mzML'
  end

  # def test_read
  #   [true, false].repeated_permutation(3) do |rs, it, ui|
  #     next if rs
      # for func in [MzML, read, chain,
      #   lambda x, **kw: chain.from_iterable([x], **kw), PreIndexedMzML]:

      # r = Mzml::MzML.new(@path, 'read_schema' => rs, 'iterative' => it, 'use_index' => ui)
      # assert_equal Mzml_spectra, r.to_s
    #   r = Mzml::MzMLtest.new(@path)
    #   assert_equal Mzml_spectra, r
    # end

  # end

  def test_mp_read
    key = Itemgetter.new('index')
    f = Mzml::MzMLtest.new(@path)
    assert_equal Mzml_spectra.to_a.sort_by{ _1[key] }, f.map().to_a.sort_by{ _1[key] }
  end

  def test_mp_requires_index
    r = Mzml::MzMLtest.new(@path, 'use_index' => false)
    p [49, r.map]
    assert r.map
  end


end
