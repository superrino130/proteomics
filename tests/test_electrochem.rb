require 'minitest/autorun'
require_relative '../rbteomics/electrochem'

class TestEctrochem < Minitest::Test
  def setup
     
  end
  
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
    # assert charge(
    #   ['A','A','A'],
    #   5.0,
    #   **{
    #     'pK' => {
    #       'H-' => [[9.0, 1]],
    #       '-OH' => [[8.0, -1]]
    #     },
    #     'pK_nterm' => {
    #       'H-' => {'A' => [[3.0, 1]]}
    #     }
    #   }
    # )
    assert (charge(['H-','A','A','A','-OH'], 0.0) - 1.0).abs < 0.01
    assert (charge(['H-','A','A','A','-OH'], 14.0) + 1.0).abs < 0.01
    assert (charge(['H-','A','A','A','-OH'], (2.34 + 9.69) / 2.0)).abs < 0.01
  end
end