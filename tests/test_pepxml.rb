# # # from os import path
# # # import pyteomics
# # # pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
# # from itertools import product
# import unittest
# from pyteomics.pepxml import PepXML, read, chain, filter
# from data import pepxml_results

require 'minitest/autorun'
require_relative '../rbteomics/pepxml'
require_relative 'data'
require 'set'

class TestPepxml < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'tests/test.pep.xml'

    @_kw = {'full_output' => false, 'fdr' => 1,
      'key' =>  lambda { |x| x['search_hit'].map{ |sh| [sh['search_score']['expect'], 1].min } }
    }
  end

  def test_ReadPSM
    [true, false].repeated_permutation(2).each do |rs, it|
      r = Pepxml::PepXML.new(@path, 'read_schema' => rs, 'iterative' => it)
      assert_equal r, Pepxml_results
    end
    # [true, false].repeated_permutation(2).each do |rs, it|
    #   r = Pepxml.read(@path, 'read_schema' => rs, 'iterative' => it)
    #   assert_equal r, Pepxml_results
    # end
    # [true, false].repeated_permutation(2).each do |rs, it|
    #   r = Pepxml.chain(@path, 'read_schema' => rs, 'iterative' => it)
    #   assert_equal r, Pepxml_results
    # end
    # [true, false].repeated_permutation(2).each do |rs, it|
      # [PepXML].each do |func|
      # [PepXML.new, read, chain,
      #   lambda { |x, **kw| chain },
      #   lambda { |x, **kw| filter.call() }
      # ].each do |func|
        # r = Pepxml.read(@path, 'read_schema' => rs, 'iterative' => it)
        # assert_equal r, Pepxml_results
      # end
  # end
  end
end
