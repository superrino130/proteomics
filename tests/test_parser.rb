require 'minitest/autorun'
require_relative '../rbteomics/parser'

class TestMass < Minitest::Unit::TestCase
  def setup
    @simple_sequences = 10.times.map{ rand(1..20).times.map{ ('A'..'Z').to_a.sample }.join('') }
    @labels = ['A', 'B', 'C', 'N', 'X']
    @extlabels = @labels[0..-1]
    @potential = {'pot' => ['X', 'A', 'B'], 'otherpot' => ['A', 'C'], 'N-' => ['N'], '-C' => ['C']}
    @constant = {'const' => ['B']}
    @extlabels.concat(['pot', 'otherpot', 'const', '-C', 'N-'])
  end

  def test_parse_simple
    @simple_sequences.each do |seq|
      assert seq == parse(seq, 'labels' => 'upcase').join('')
    end
  end

  def test_parse
    assert [['P'], ['E'], ['P'], ['T'], ['I'], ['D'], ['E']] == parse('PEPTIDE', splitflg: true)
    assert ['P', 'E', 'P', 'T', 'I', 'D', 'E'] == parse('H-PEPTIDE')
    ['PEPTIDE', 'H-PEPTIDE', 'PEPTIDE-OH', 'H-PEPTIDE-OH'].each do |seq|
      assert ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH'] == parse(seq, show_unmodified_termini: true)
    end
    assert ['T', 'E', 'pS', 'T', 'oxM'] == parse('TEpSToxM', 'labels' => STD_labels + ['pS', 'oxM'])
    assert [['H-', 'z', 'P'], ['E'], ['P'], ['z', 'T'], ['I'], ['D'], ['z', 'E', '-OH']] == parse('zPEPzTIDzE', show_unmodified_termini: true, splitflg: true, 'labels' => STD_labels + ['z'])
  end

  def test_tostring
    @simple_sequences.each do |seq|
      assert seq == tostring(parse(seq, 'labels' => 'uppercase'))
      p [33, seq, tostring(parse(seq, show_unmodified_termini: true, splitflg: true, 'labels' => 'uppercase'), show_unmodified_termini: false)]
      assert seq == tostring(parse(seq, show_unmodified_termini: true, splitflg: true, 'labels' => 'uppercase'), show_unmodified_termini: false)
    end
    
  end

end