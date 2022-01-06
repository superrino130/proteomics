require 'minitest/autorun'
require_relative '../rbteomics/electrochem'
require_relative '../rbteomics/auxiliary/structures'

class TestEctrochem < Minitest::Test
  def setup; end
  
  def test_charge_calculations_str
    assert charge(
      'AAA',
      5.0,
      **{
        'pK' => {
          'H-' => [[9.0, 1]],
          '-OH' => [[8.0, -1]]
        },
        'pK_nterm' => {
          'H-' => {'A' => [[3.0, 1]]}
        }
      }
    ).abs < 0.01
    assert (charge('H-AAA-OH', 0.0) - 1.0).abs < 0.01
    assert (charge('H-AAA-OH', 14.0) + 1.0).abs < 0.01
    assert (charge('H-AAA-OH', (2.34 + 9.69) / 2.0)).abs < 0.01
  end

  def test_charge_calculations_list
    e = assert_raises PyteomicsError do
      charge(
        ['A','A','A'],
        5.0,
        **{
          'pK' => {
            'H-' => [[9.0, 1]],
            '-OH' => [[8.0, -1]]
          },
          'pK_nterm' => {
            'H-' => {'A' => [[3.0, 1]]}
          }
        }
      )
    end
    assert e.message.include?("Parsed sequences must contain terminal groups at 0-th and last positions.")

    assert (charge(['H-','A','A','A','-OH'], 0.0) - 1.0).abs < 0.01
    assert (charge(['H-','A','A','A','-OH'], 14.0) + 1.0).abs < 0.01
    assert (charge(['H-','A','A','A','-OH'], (2.34 + 9.69) / 2.0)).abs < 0.01
  end

  def test_charge_calculations_dict
    e = assert_raises PyteomicsError do
      charge(
        {'H-' => 1, '-OH' => 1, 'E' => 1},
        7,
        **{'pK_nterm' => {
          'H-' => {'A' => [[9.0, 1]]}}}
      )
    end
    assert e.message.include?("Two terminal residues must be present in peptide (designated as 'ntermX' and 'ctermX', where 'X' is the one-letter residue label). Use 'term_aa=true' when calling 'parser.amino_acid_composition'.")

    assert (charge({'A' => 3, 'H-' => 1, '-OH' => 1}, 14.0) + 1.0).abs < 0.01
    assert (charge({'A' => 1, 'H-' => 1, '-OH' => 1, 'ntermB' => 1, 'ctermA' => 1},
      14.0,
      **{
        'pK' => {'H-' => [[9.0, 1]], '-OH': [[8.0, -1]]},
        'pK_nterm' => {'H-' => {'A' => [[3.0, 1]], 'B' => [[3.0, 1]]}}})).abs < 0.01
    e = assert_raises PyteomicsError do
      charge({
        'A' => 1, 'H-' => 1, '-OH' => 1, 'ctermA' => 1},
        14.0,
        **{
          'pK' => {'H-' => [[9.0, 1]], '-OH' => [[8.0, -1]]},
          'pK_nterm' => {'H-' => {'A' => [[3.0, 1]]}}}
      )
    end
    assert e.message.include?("Two terminal residues must be present in peptide (designated as 'ntermX' and 'ctermX', where 'X' is the one-letter residue label). Use 'term_aa=true' when calling 'parser.amino_acid_composition'.")
    e = assert_raises PyteomicsError do
      charge({
        'A' => 1, 'H-' => 1, '-OH' => 1, 'ntermA' => 1},
        14.0,
        **{
          'pK' => {'H-' => [[9.0, 1]], '-OH' => [[8.0, -1]]},
          'pK_nterm' => {'H-' => {'A' => [[3.0, 1]]}}}
      )
    end
    assert e.message.include?("Two terminal residues must be present in peptide (designated as 'ntermX' and 'ctermX', where 'X' is the one-letter residue label). Use 'term_aa=true' when calling 'parser.amino_acid_composition'.")
    e = assert_raises PyteomicsError do
      charge({
        'A' => 1, 'H-' => 1, '-OH' => 1, 'ntermA' => 2, 'ctermA' => 1},
        14.0,
        **{
          'pK' => {'H-' => [[9.0, 1]], '-OH' => [[8.0, -1]]},
          'pK_nterm' => {'H-' => {'A' => [[3.0, 1]]}}}
      )
    end
    assert e.message.include?("More that one N-terminal residue in {")
    e = assert_raises PyteomicsError do
      charge({
        'A' => 1, 'H-' => 1, 'ntermA' => 1, 'ctermA' => 1},
        14.0,
        **{
          'pK' => {'H-' => [[9.0, 1]], '-OH' => [[8.0, -1]]},
          'pK_nterm' => {'H-' => {'A' => [[3.0, 1]]}}}
      )
    end
    assert e.message.include?("Peptide must have two explicit terminal groups")
  end

  def test_pI_calculations
    p [110, pI('H-AAA-OH'), 2.34+9.69]
    assert (pI('H-AAA-OH') - (2.34 + 9.69) / 2.0).abs < 0.01
  end

  def test_pI_precision
    pI_best = pI('PEPTIDE', precision_pI: 1e-15)
    16.times do |i|
      precision = 10 ** -i
      assert (pI('PEPTIDE', precision_pI: precision) - pI_best).abs < precision
    end
  end

  def test_charge_input
    (0...14).each do |i|
      a = charge('H-ACDEFGH-OH', i)
      b = charge(['H-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', '-OH'], i)
      assert (a - b).round(7).zero?
    end

    (0...14).each do |i|
      a = charge('H-ACDEFGH-OH', i)
      b = charge({'H-' => 1, 'A' => 1, 'C' => 1, 'D' => 1, 'E' => 1, 'F' => 1, 'G' => 1, 'H' => 1, '-OH' => 1}, i)
      assert (a - b).round(7).zero?
    end
  end
end