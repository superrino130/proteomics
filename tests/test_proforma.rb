require 'minitest/autorun'
require_relative '../rbteomics/proforma'
require 'set'

class TestPepxml < Minitest::Test
  def setup
    @maxDiff = nil
  end

  # def test_complicated_short
  #   complicated_short = "<[Carbamidomethyl]@C><13C>[Hydroxylation]?{HexNAc}[Hex]-ST[UNIMOD:Oxidation](EPP)[+18.15]ING"
  #   tokens, properties = Proforma::ProForma.parse(complicated_short)
    # p [13, tokens,properties]
    # assert_equal tokens.size, 8
    # assert_equal properties['n_term'].size == 1
    # assert_equal properties['n_term'][0], 'Hex'
    # assert_equal properties['intervals'].size, 1
    # assert_equal properties['intervals'][0], TaggedInterval.new(2, 5, [MassModification(18.15)])
    # assert_equal properties['isotopes'].size, 1
    # assert_equal properties['isotopes'][0], StableIsotope.new("13C")
    # assert_equal properties['fixed_modifications'][0], ModificationRule.new(
    #     GenericModification.new('Carbamidomethyl', nil, nil), ['C'])
    # assert_equal to_proforma(tokens, **properties), complicated_short
    # assert_equal ProForma.new(tokens, properties).mass.round(3), 1210.5088.round(3)
  # end

  def test_range
    seq = "PRQT(EQC[Carbamidomethyl]FQRMS)[+19.0523]ISK"
    parsed = Proforma::ProForma.parse(seq)
    assert_equal parsed.to_s, seq
    chunk = parsed[0...6]
    assert chunk.intervals
  end
end
