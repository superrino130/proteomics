require 'minitest/autorun'
require_relative '../rbteomics/mass/mass'

class TestMass < Minitest::Unit::TestCase
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

    @mass_H = mass.nist_mass['H'][0][0]
    @mass_O = mass.nist_mass['O'][0][0]
    @test_aa_mass = {'X' => 1.0, 'Y' => 2.0, 'Z' => 3.0}
    @random_peptides = (0...10).to_a.map{ (0...20).to_a.map{ 'XYZ'.split('')[rand(3)] } }.join('')

    @aa_comp = {
            'X' =>   mass.Composition({'A' => 1}, mass_data=@mass_data),
            'Y' =>   mass.Composition({'B' => 1}, mass_data=@mass_data),
            'Z' =>   mass.Composition({'C' => 1}, mass_data=@mass_data),
            'F' =>   mass.Composition({'F' => 1}, mass_data=@mass_data),
            'H-' =>  mass.Composition({'D' => 1}, mass_data=@mass_data),
            '-OH' => mass.Composition({'E' => 1}, mass_data=@mass_data),
        }

    @ion_comp = {
      'M' => mass.Composition({}, mass_data=@mass_data),
      'a' => mass.Composition({'A' => -1}, mass_data=@mass_data)
    }

    @mods = {'a' => mass.Composition({'A' => 1}), 'b' => mass.Composition({'B' => 1})}
    @d = {atom: 1 for atom in 'ABCDE'}
  end

  def test1

  end
end