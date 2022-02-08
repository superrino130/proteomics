# from os import path
# import pyteomics
# pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
# from pyteomics.mass import unimod
# import unittest

require 'minitest/autorun'
require_relative '../rbteomics/mass/unimod'
require_relative '../rbteomics/auxiliary/constants'
require 'zlib'

class TestUnimod < Minitest::Test
  def setup
    @handle = Unimod::Unimod.new
    
  end

  def test_modifications_have_composition

  end
end