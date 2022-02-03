require 'minitest/autorun'
require_relative '../rbteomics/mzid'
require_relative '../rbteomics/auxiliary/structures'
require_relative 'data'

class TestMzid < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'test.mzid'
  end
  
  def test_ReadPSM

  end

  def test_unit_info

  end

  def test_structure_normalization

  end

  def test_map
    mz = MZID::MzIdentML.new(@path).to_enum

    cnt = 0
    loop do
      mz.next
      cnt += 1
    end
    assert_equal Mzid_spectra[[true, true]].size, MZID::MzIdentML.new(@path).size
  end

  def test_iterfind_map

  end
end
