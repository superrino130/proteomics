require 'minitest/autorun'
require_relative '../rbteomics/parser'
require_relative '../rbteomics/usi'
require_relative '../rbteomics/auxiliary/structures'
require_relative 'data'
require 'set'

class TestUsi < Minitest::Test
  def test_parse
    usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
    inst = Usi::USI.new(usi_str).parse
    assert_equal inst.USI, usi_str
    assert_equal inst.protocol, 'mzspec'
    assert_equal inst.dataset, "MSV000085202"
    assert_equal inst.datafile, "210320_SARS_CoV_2_T"
    assert_equal inst.scan_identifier_type, "scan"
    assert_equal inst.scan_identifier, "131256"
    assert_equal inst.interpretation, nil
  end
end