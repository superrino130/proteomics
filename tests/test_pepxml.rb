# # # from os import path
# # # import pyteomics
# # # pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
# # from itertools import product
# import unittest
# from pyteomics.pepxml import PepXML, read, chain, filter
# from data import pepxml_results

require 'minitest/autorun'
require_relative '../rbteomics/pepxml'
require 'set'

class TestPepxml < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'test.pep.xml'

    @_kw = {'full_output' => false, 'fdr' => 1,
      'key' =>  lambda { |x| x['search_hit'].map{ |sh| [sh['search_score']['expect'], 1].min } }
    }
  end

  def testReadPSM
    [true, false].repeated_permutation(2).each do |rs, it|
      [PepXML.new, read, chain,
        lambda { |x, **kw| chain }
      ]
    end
  end
end
