require 'minitest/autorun'
require_relative '../rbteomics/mzid'
require_relative '../rbteomics/xml'
require_relative '../rbteomics/auxiliary/structures'
require_relative '../rbteomics/auxiliary/file_helpers'
require_relative 'data'

class TestMzid < Minitest::Test
  include Xml
  def setup
    @maxDiff = nil
    @path = 'tests/test.mzid'
  end
  
  def test_ReadPSM

  end

  def test_unit_info

  end

  def test_structure_normalization

  end

  def test_map
    mz = Mzid::MzIdentML.new(@path)
    p Mzid_spectra[[true, true]].size

    p mz.map
    cnt = 0
    loop do
      mz.next
      cnt += 1
    end
    assert_equal Mzid_spectra[[true, true]].size, MZID::MzIdentML.new(@path).size
  end
end
