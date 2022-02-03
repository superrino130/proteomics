require 'minitest/autorun'
require_relative '../rbteomics/electrochem'
require_relative '../rbteomics/tandem'
require_relative '../rbteomics/version'
require_relative '../rbteomics/auxiliary/structures'

require 'numpy'
require 'pandas'

PSMS = []
alph = ('A'..'Z').to_a + ('a'..'z').to_a
52.times do |i|
  PSMS << [i, alph[i], 0.01 + i * 0.001]
end

class TestAuxiliary < Minitest::Test
  @@key = lambda { |x| x[0] }
  @@is_decoy = lambda { |x| x[1] == x[1].downcase }
  @@pep = lambda { |x| x[2] }

  def setup
    @psms = PSMS
    Numpy.random.shuffle(@psms)
  end
  
  def _run_check(q, formula)
      assert Numpy.allclose(q['q'][0...26], 0)
    if formula == 2
      assert Numpy.allclose(q['q'][26..], 2 * Numpy.arange(1.0, 27.0) / (26 + Numpy.arange(1, 27)))
    else
      assert Numpy.allclose(q['q'][26..], Numpy.arange(1.0, 27.0) / 26)
    end
    assert Numpy.allclose(q['is decoy'][0...26], 0)
    assert Numpy.allclose(q['is decoy'][26..], 1)
    assert Numpy.allclose(q['score'], Numpy.arange(52))
    setup
    spsms = @psms.sort_by(&key)
    assert Numpy.allclose(spsms.map(&@@is_decoy), q['is decoy'])
    assert Numpy.allclose(spsms.map(&@@key), q['score'])
    setup
  end

  def _run_check_pep(q)
    assert Numpy.allclose(q['q'], Numpy.arange(0.01, 0.036, 0.0005))
    setup
  end

  def test_qvalues
    q = Target_decoy.qvalues(@psms, key: @@key, is_decoy: @@is_decoy, remove_decoy: true)
    assert Numpy.allclose(q['q'], 0)
    assert Numpy.allclose(q['is decoy'], 0)
    assert Numpy.allclose(q['score'], Numpy.arange(26))
  end
end