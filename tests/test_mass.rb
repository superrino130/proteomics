require 'minitest/autorun'
require_relative '../rbteomics/mass/mass'
require_relative '../rbteomics/auxiliary/constants'

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
      assert (fast_mass2(pep, **{'aa_mass' => @test_aa_mass}) -
        (@test_aa_mass.sum{ |aa, m| pep.count(aa) * m } + @mass_H * 2.0 + @mass_O)).round(7).zero?
    end
  end

  def test_Composition_dict
    p [Composition.new(@d, 'mass_data' => @mass_data), @d]
    assert Composition.new(@d, 'mass_data' => @mass_data) == @d
  end
end