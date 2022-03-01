require 'minitest/autorun'
require_relative '../rbteomics/parser'
require_relative '../rbteomics/usi'
require_relative '../rbteomics/auxiliary/structures'
require_relative 'data'
require 'set'

class TestProxi < Minitest::Test
  def test_request
    usi_str = "mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256"
    response = Usi.proxi(usi_str, 'backend' =>'peptide_atlas')

    assert usi_proxi_data.keys <= response.keys

    # for a, b in zip(response['m/z array'], usi_proxi_data['m/z array']):
    #     self.assertAlmostEqual(a, b, 3)

    # for a, b in zip(response['intensity array'], usi_proxi_data['intensity array']):
    #     self.assertAlmostEqual(a, b, 3)
    # end
  end

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