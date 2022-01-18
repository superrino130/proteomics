require 'minitest/autorun'
require_relative '../rbteomics/mass/mass'
require_relative '../rbteomics/auxiliary/constants'
require 'zlib'

class TestMass < Minitest::Test
  def setup
    @mass_data = {
      'A' => {0 => [1.0, 1.0],
             1 => [1.0, 0.5],
             2 => [2.0, 0.5]},
      'B' => {0 => [2.0, 1.0],
             2 => [2.0, 0.5],
             3 => [3.0, 0.5]},
      'C' => {0 => [3.0, 1.0],
             3 => [3.0, 0.5],
             4 => [4.0, 0.5]},
      'D' => {0 => [4.0, 1.0],
             4 => [4.0, 0.5],
             5 => [5.0, 0.5]},
      'E' => {0 => [5.0, 1.0],
             5 => [5.0, 0.5],
             6 => [6.0, 0.5]},
      'F' => {0 => [6.0, 1.0],
             6 => [6.0, 0.7],
             7 => [7.0, 0.3]},
      'H+'=> {0 => [5.0, 1.0],
             5 => [5.0, 1.0]},
    }

    @mass_H = $_nist_mass['H'][0][0]
    @mass_O = $_nist_mass['O'][0][0]
    @test_aa_mass = {'X' => 1.0, 'Y' => 2.0, 'Z' => 3.0}
    @random_peptides = (0...10).to_a.map{ (0...20).to_a.map{ 'XYZ'.split('').sample }.join('') }

    @aa_comp = {
            'X' =>   Composition.new({'A' => 1}, 'mass_data' => @mass_data),
            'Y' =>   Composition.new({'B' => 1}, 'mass_data' => @mass_data),
            'Z' =>   Composition.new({'C' => 1}, 'mass_data' => @mass_data),
            'F' =>   Composition.new({'F' => 1}, 'mass_data' => @mass_data),
            'H-' =>  Composition.new({'D' => 1}, 'mass_data' => @mass_data),
            '-OH' => Composition.new({'E' => 1}, 'mass_data' => @mass_data),
        }

    @ion_comp = {
      # 'M' => Composition.new({}, 'mass_data' => @mass_data),
      'a' => Composition.new({'A' => -1}, 'mass_data' => @mass_data)
    }

    @mods = {'a' => Composition.new({'A' => 1}), 'b' => Composition.new({'B' => 1})}
    @d = {'A' => 1, 'B' => 1, 'C' => 1, 'D' => 1, 'E' => 1}
  end

  def test_fast_mass
    @random_peptides.each do |pep|
      assert (fast_mass(pep, **{'aa_mass' => @test_aa_mass}) -
        (@test_aa_mass.sum{ |aa, m| pep.count(aa) * m } + @mass_H * 2.0 + @mass_O)).round(7).zero?
    end
  end

  def test_fast_mass2
    @random_peptides.each do |pep|
      assert_equal (fast_mass2(pep, 'aa_mass' => @test_aa_mass) -
        (@test_aa_mass.sum{ |aa, m| pep.count(aa) * m } + @mass_H * 2.0 + @mass_O)).round(7), 0.0
    end
  end

  def test_Composition_dict
    assert_equal Composition.new(@d, 'mass_data' => @mass_data), @d
  end

  def test_Composition_formula
    assert_equal @d, Composition.new('formula' => 'ABCDE', 'mass_data' => {"A"=>{0=>[1.0, 1.0]}, "B"=>{0=>[1.0, 1.0]}, "C"=>{0=>[1.0, 1.0]}, "D"=>{0=>[1.0, 1.0]}, "E"=>{0=>[1.0, 1.0]}})
  end

  def test_Composition_seq
    assert_equal @d, Composition.new('sequence' => 'XYZ', 'aa_comp' => @aa_comp)
  end

  def test_Composition_pseq
    assert_equal Composition.new('split_sequence' => [['X'], ['Y'], ['Z']], 'aa_comp' => @aa_comp),
      {'A' => 1, 'B' => 1, 'C' => 1}
  end

  def test_Composition_sum
    assert_equal Composition.new('sequence' => 'XXY', 'aa_comp' => @aa_comp) + Composition.new('sequence' => 'YZZ', 'aa_comp' => @aa_comp),
      {'A' => 2, 'B' => 2, 'C' => 2, 'D' => 2, 'E' => 2}
  end

  def test_Composition_sub
    assert_equal Composition.new() - Composition.new('sequence' => 'XYZ', 'aa_comp' => @aa_comp),
      {'A' => -1, 'B' => -1, 'C' => -1, 'D' => -1, 'E' => -1}
  end

  def test_Composition_mul
    # Integer * Hash は未定義
    assert_equal Composition.new('sequence' => 'XYZ', 'aa_comp' => @aa_comp) * 2,
      {'A' => 2, 'B' => 2, 'C' => 2, 'D' => 2, 'E' => 2}
  end

  def test_Composition_positional
    ac = @aa_comp.dup
    ac.merge!(@mods)
    assert_equal Composition.new('aXbYZ', 'aa_comp' => ac), {'A' => 2, 'B' => 2, 'C' => 1, 'D' => 1, 'E' => 1}
    assert_equal Composition.new('AB2C3', 'mass_data' => @mass_data), {'A' => 1, 'B' => 2, 'C' => 3}
  end

  def test_calculate_mass
    assert_equal calculate_mass('formula' => 'ABCDE', 'mass_data' => @mass_data),
      'ABCDE'.split('').map{ @mass_data[_1][0][0] }.sum
    assert_equal calculate_mass('sequence' => 'XYZ',
                            'aa_comp' => @aa_comp,
                            'mass_data' => @mass_data),
      'ABCDE'.split('').map{ @mass_data[_1][0][0] }.sum
    assert_equal calculate_mass('parsed_sequence' => ['H-', 'X', 'Y', 'Z', '-OH'],
       'aa_comp' => @aa_comp,
       'mass_data' => @mass_data),
       'ABCDE'.split('').map{ @mass_data[_1][0][0] }.sum

    a = []
    'ABCDE'.split('').each do |x|
      @mass_data[x].keys.each do |y|
          a << @mass_data[x][y][0] * @mass_data[x][y][1] if y != 0
      end
    end
    assert_equal calculate_mass('formula' => 'ABCDE',
      'average' => true,
      'mass_data' => @mass_data),
      a.sum
    
      [1, 2, 3].each do |charge|
        assert_equal calculate_mass('formula' => 'ABCDE', 'ion_type' => 'M', 'charge' => charge, 'mass_data' => @mass_data),
          calculate_mass('formula' => 'ABCDE' + "H+#{charge}", 'mass_data' => @mass_data)
      
        assert_equal calculate_mass('formula' => 'ABCDE', 'ion_type' => 'M', 'charge' => charge, 'mass_data' => @mass_data),
          (calculate_mass('formula' => 'ABCDE', 'mass_data' => @mass_data) + @mass_data['H+'][0][0] * charge) / charge

        e = assert_raises PyteomicsError do
          calculate_mass('formula' => "ABCDEH+#{charge}",
            'ion_type' => 'M', 'charge' => charge, 'mass_data' => @mass_data)
        assert e.message.include?("Charge is specified both by the number of protons and 'charge' in kwargs")
      end
    end

    @random_peptides.each do |pep|
      assert_equal calculate_mass('sequence' => pep, 'aa_comp' => @aa_comp, 'mass_data' => @mass_data, 'ion_comp' => @ion_comp),
        calculate_mass('parsed_sequence' => parse(pep, 'labels' => ['X', 'Y', 'Z'], show_unmodified_termini: true),
          'aa_comp' => @aa_comp, 'mass_data' => @mass_data, 'ion_comp' => @ion_comp)
    end
  end

  def test_most_probable_isotopic_composition
    assert_equal most_probable_isotopic_composition('formula' => 'F', 'mass_data' => @mass_data),
      [Composition.new(['F[6]' => 1], ['F[7]'=> 0], **{'mass_data' => @mass_data}), 0.7]

    assert_equal most_probable_isotopic_composition('formula' => 'F10', 'mass_data' => @mass_data),
      [Composition.new(*['F[6]' => 7, 'F[7]' => 3], **{'mass_data' => @mass_data}), (0.3)**3 * (0.7)**7 * 120]
    
    
    assert_equal most_probable_isotopic_composition('formula' => 'A20F10', 'elements_with_isotopes' => ['F'], 'mass_data' => @mass_data),
        [Composition.new(*['A' => 20, 'F[6]' => 7, 'F[7]' => 3], 'mass_data' => @mass_data), (0.3)**3 * (0.7)**7 * 120]
  end

  def test_isotopic_composition_abundance
    (1...10).to_a.each do |peplen|
      assert_equal isotopic_composition_abundance('formula' => 'F[6]' * peplen, 'mass_data' => @mass_data).round(7),
        @mass_data['F'][6][1] ** peplen.round(7)

      assert_equal isotopic_composition_abundance('formula' => 'AF[6]' * peplen, 'mass_data' => @mass_data).round(7),
        @mass_data['F'][6][1] ** peplen.round(7)

      assert_equal isotopic_composition_abundance('formula' => 'A[1]F[6]' * peplen, 'mass_data' => @mass_data).round(7),
          ((@mass_data['A'][1][1] * @mass_data['F'][6][1]) ** peplen).round(7)
    end
  end

  # def test_Unimod_mass
  #   db = Unimod.new(source: Zlib::GzipReader.open('unimod.xml.gz'))
  #   db.mods.each do |x|
  #     assert 0.00001 < (x['mono_mass'] - calculate_mass(x['composition'], 'mass_data' = > db.mass_data)).abs
  #   end
  # end

  def test_Unimod_methods

  end

  def test_nist_mass
    assert NIST_mass.values.map{ |g| g[0][1] - 1 < 1e-6 }.all?

    NIST_mass.values.each do |g|
      s = g.sum{ |num, p| num != 0 ? p[1] : 0 }
      assert (s - 1).abs < 1e-6 || s.abs < 1e-6
    end
  end

  def test_composition_objects_are_pickleable
    dict_ = Composition.new(@d, 'mass_data' => @mass_data)
    formula = Composition.new('formula' => 'ABCDE',
      'mass_data' => 'ABCDE'.split('').map{ [_1, {0 => {0 => [1.0, 1.0]}}] }.to_h)
    sequence = Composition.new('sequence' => 'XYZ', 'aa_comp' => @aa_comp)
    parsed_sequence = Composition.new('parsed_sequence' => ['X', 'Y', 'Z'],
      'aa_comp' => @aa_comp)
    split_sequence = Composition.new('split_sequence' => [['X'], ['Y'], ['Z']],
      'aa_comp' => @aa_comp)
    
    assert_equal dict_, Marshal.load(Marshal.dump(dict_))
    assert_equal formula, Marshal.load(Marshal.dump(formula))
    assert_equal sequence, Marshal.load(Marshal.dump(sequence))
    assert_equal parsed_sequence, Marshal.load(Marshal.dump(parsed_sequence))
    assert_equal split_sequence, Marshal.load(Marshal.dump(split_sequence))
  end

  def test_aa_mass
    h2o = calculate_mass('formula' => 'H2O')
    STD_aa_mass.each do |aa, m|
      assert_equal m + h2o, fast_mass(aa)
    end
  end

  def test_isotopologues
    peptide = 'XYF'
    states = [{'F[6]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}, {'F[7]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}]
    abundances = [0.7, 0.3]
    kw_common = {'elements_with_isotopes' => 'F', 'aa_comp' => @aa_comp, 'mass_data' => @mass_data}
    kwlist = [
      {},
      {'sequence' => 'XYF'},
      {'parsed_sequence' => parse('XYF', show_unmodified_termini: true)},
      {'split_sequence' => parse('XYF', show_unmodified_termini: true, split: true)},
      {'formula' => 'ABDEF'},
      {'composition' => Composition.new('sequence' => 'XYF', 'aa_comp' => @aa_comp)}]
    arglist = [[peptide], [], [], [], [], []]
    arglist.zip(kwlist).each do |args, kw|
      kwargs = kw_common.dup
      kwargs.merge(kw)
      isotopologues = isotopologues(*args, **kwargs)
      isotopologues.each do |state|
        i = states.index(state)
        assert i != -1
        assert_equal abundances[i].round(7), isotopic_composition_abundance(state,
          'aa_comp' => @aa_comp, 'mass_data' => @mass_data).round(7)
      end
    end
  end

  def test_isotopologues_with_abundances
    peptide = 'XYF'
    states = [{'F[6]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}, {'F[7]': 1, 'A': 1, 'B': 1, 'D': 1, 'E': 1}]
    abundances = [0.7, 0.3]
    isotopologues(peptide, 'elements_with_isotopes' => 'F',
      'aa_comp' => @aa_comp, 'mass_data' => @mass_data, 'report_abundance' => true).each do |state, abundance|
      i = states.index(state)
      assert i != -1
      assert_equal abundances[i].round(7), abundance.round(7)
    end
  end

  def test_std_aa_mass
    STD_aa_mass.each do |key, value|
      assert_equal value.round(4), calculate_mass('parsed_sequence' => [key]).round(4)
    end
  end
end
