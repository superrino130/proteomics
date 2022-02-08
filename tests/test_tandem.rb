require 'minitest/autorun'
require_relative '../rbteomics/tandem'
require_relative 'data'

class TestTandem < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'test.t.xml'
  end

  def testReadPSM
    # Not started
  end

  def test_df
    df = Tandem.dataframe(@path)
    assert_equal df.shape, [1, 29]
  end
end