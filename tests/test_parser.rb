require 'minitest/autorun'
require_relative '../rbteomics/parser'
require 'set'

class TestMass < Minitest::Test
  def setup
    @simple_sequences = 10.times.map{ rand(1..20).times.map{ ('A'..'Z').to_a.sample }.join('') }
    @labels = ['A', 'B', 'C', 'N', 'X']
    @extlabels = @labels[0..-1].dup
    @potential = {'pot' => ['X', 'A', 'B'], 'otherpot' => ['A', 'C'], 'N-' => ['N'], '-C' => ['C']}
    @constant = {'const' => ['B']}
    @extlabels.concat(['pot', 'otherpot', 'const', '-C', 'N-'])
  end

  def test_parse_simple
    @simple_sequences.each do |seq|
      assert_equal seq, parse(seq, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ').join('')
    end
  end

  def test_parse
    assert_equal [['P'], ['E'], ['P'], ['T'], ['I'], ['D'], ['E']], parse('PEPTIDE', split: true)
    assert_equal ['P', 'E', 'P', 'T', 'I', 'D', 'E'], parse('H-PEPTIDE')
    ['PEPTIDE', 'H-PEPTIDE', 'PEPTIDE-OH', 'H-PEPTIDE-OH'].each do |seq|
      assert_equal ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH'], parse(seq, show_unmodified_termini: true)
    end
    assert_equal ['T', 'E', 'pS', 'T', 'oxM'], parse('TEpSToxM', 'labels' => STD_labels + ['pS', 'oxM'])
    assert_equal [['H-', 'z', 'P'], ['E'], ['P'], ['z', 'T'], ['I'], ['D'], ['z', 'E', '-OH']], parse('zPEPzTIDzE', show_unmodified_termini: true, split: true, 'labels' => STD_labels + ['z'])
  end

  def test_tostring
    @simple_sequences.each do |seq|
      assert_equal seq, tostring(parse(seq, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
      assert_equal seq, tostring(parse(seq, show_unmodified_termini: true, split: true, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'), show_unmodified_termini: false)
    end 
  end

  def test_amino_acid_composition_simple
    @simple_sequences.each do |seq|
      comp = amino_acid_composition(seq, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
      seq.split('').to_set.each do |aa|
        assert_equal seq.count(aa), comp[aa]
      end
    end
  end

  def test_amino_acid_composition
    @simple_sequences.each do |seq|
      comp = amino_acid_composition(seq, term_aa: true, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
      comp_default = amino_acid_composition(seq, 'labels' => 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
      assert_equal 1, comp['nterm' + seq[0]]
      if seq.size > 1
        assert_equal 1, comp['cterm' + seq[-1]]
      end
      assert_equal comp_default.values.sum, comp.values.sum
    end
  end

  # def test_cleave
  #   assert_equal xcleave('PEPTIDEKS', EXPASY_rules['trypsin']), [[0, 'PEPTIDEK'], [8, 'S']]
  # end

  # def test_cleave_semi

  # end

  # def test_cleave_min_length

  # end

  # def test_num_sites
  #   assert_equal num_sites('RKCDE', 'K'), 1
  #   assert_equal num_sites('RKCDE', 'E'), 0
  #   assert_equal num_sites('RKCDE', 'R'), 1
  #   assert_equal num_sites('RKCDE', 'Z'), 0
  # end

  def test_isoforms_simple
    assert_equal(
      isoforms('PEPTIDE', 'variable_mods' => {'xx' => ['A', 'B', 'P', 'E']}),
      ['PEPTIDE', 'PEPTIDxxE', 'PExxPTIDE', 'PExxPTIDxxE', 'PxxEPTIDE', 'PxxEPTIDxxE', 'PxxExxPTIDE',
      'PxxExxPTIDxxE', 'xxPEPTIDE', 'xxPEPTIDxxE', 'xxPExxPTIDE', 'xxPExxPTIDxxE', 'xxPxxEPTIDE',
      'xxPxxEPTIDxxE', 'xxPxxExxPTIDE', 'xxPxxExxPTIDxxE'])
  end

  def test_isoforms_fixed_simple
    assert_equal(
      isoforms('PEPTIDE', 'fixed_mods' => {'n-' => true, '-c' => true, 'x' => ['P', 'T']}),
      ['n-xPExPxTIDE-c'])
  end

  def test_isoforms_simple_2
    assert_equal(isoforms('PEPTIDE', 'variable_mods' => {'x' => 'T', 'y' => 'T'}),
            ['PEPTIDE', 'PEPxTIDE', 'PEPyTIDE'])
  end

  def test_isoforms_universal
    assert_equal(Set.new(isoforms('PEPTIDE', 'variable_mods' => {'xx-' => true})), Set.new(['PEPTIDE', 'xx-PEPTIDE']))
    assert_equal(Set.new(isoforms('PEPTIDE', 'variable_mods' => {'-xx' => true})), Set.new(['PEPTIDE', 'PEPTIDE-xx']))
    @simple_sequences.each do |seq|
      assert_equal(isoforms(seq, 'variable_mods' => {'x' => true}).size, 2**seq.size)
    end
  end

  def test_isoforms_terminal
    assert_equal(Set.new(isoforms('PEPTIDE', 'variable_mods' => {'xx' => ['ntermP'], 'yy-' => 'P'})),
            Set.new(['PEPTIDE', 'xxPEPTIDE', 'yy-PEPTIDE', 'yy-xxPEPTIDE']))
  end

  def test_isoforms_len
    50.times do |j|
      l = rand(1..10)
      peptide = l.times.map{ @labels.sample }.join('')
      modseqs = isoforms(peptide, 'variable_mods' => @potential, 'fixed_mods' => @constant, 'labels' => @labels)
      pp = parse(peptide, 'labels' => @extlabels)
      n = (pp[0] == 'N' ? 1 : 0) + (pp[-1] == 'C' ? 1 : 0)
      modseqs.each do |pm|
        assert_equal(pp.size, length(pm, 'labels' => @extlabels))
      end
      assert_equal(modseqs.size, (3 ** pp.count('A')) * (2 ** (pp.count('X') + pp.count('C') + n)))
    end
  end

  # テストがおかしい
  # pyteomics の def test_isoforms_maxmods では
  # modseqs が empty な為、for ms in modseqs が走らない
  #
  # def test_isoforms_maxmods
  #   50.times do
  #     l = rand(1..10)
  #     m = rand(1..10)
  #     peptide = l.times.map{ @labels.sample }.join('')
  #     modseqs = isoforms(peptide, 'variable_mods' => @potential, 'labels' => @labels, 'max_mods' => m, 'format' => 'split')
  #     pp = parse(peptide, 'labels' => @extlabels, split: true)
  #     modseqs.each do |ms|
  #       assert_equal(pp.size, ms.size)
  #       assert(pp.zip(ms).select{ _1 != _2 }.size <= m)
  #     end
  #   end
  # end

  def test_fast_valid
    50.times do |j|
      l = rand(1..10)
      peptide = l.times.map{ @labels.sample }.join('')
      assert fast_valid(peptide, labels: Set.new(@labels))
      assert valid(peptide, 'labels' => @labels)
      assert valid(peptide)
      peptide.split('').to_set.each do |aa|
        bad = peptide.gsub(aa, 'Z')
        assert_equal fast_valid(bad, labels: Set.new(@labels)), false
        assert_equal valid(bad, 'labels' => @labels), false
      end
    end
  end

  def test_valid
    50.times do |j|
      l = rand(1..10)
      peptide = l.times.map{ @labels.sample }.join('')
      modseqs = isoforms(peptide, 'variable_mods' => @potential, 'fixed_mods' => @constant, 'labels' => @labels)
      assert_equal valid('H-' + peptide, 'labels' => @labels), false
      modseqs.each do |s|
        assert valid(s, 'labels' => @extlabels)
        Set.new(peptide.split('')).each do |aa|
          bad = s.gsub(aa, 'Z')
          assert_equal fast_valid(bad, labels: Set.new(@labels)), false
          assert_equal valid(bad, 'labels' => @labels), false
        end
      end
    end
  end
end