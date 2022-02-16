require 'minitest/autorun'
require_relative '../rbteomics/tandem'
require_relative 'data'

class TestTandem < Minitest::Test
  def setup
    @maxDiff = nil
    @path = 'tests/test.t.xml'
  end

  def test_ReadPSM
    # [tandem.TandemXML, tandem.read, tandem.chain,
    #   lambda x, **kw: tandem.chain.from_iterable([x], **kw),
    #   lambda x, **kw: tandem.filter(x, fdr=1, full_output=False),
    #   lambda x, **kw: tandem.filter.chain(x, fdr=1, full_output=False),
    #   lambda x, **kw: tandem.filter.chain.from_iterable([x], fdr=1, full_output=False)].each do |func|
    #   for it in range(2):
    #     with func(self.path, iterative=it) as r:
    #         self.assertEqual(list(r), tandem_spectra)
    #   end
    # end
    [0, 1].each do |it|
      r = Tandem::TandemXML.new(@path, 'iterative' => it)
      assert_equal r, Tandem_spectra
    end
  end

  def test_df
    df = Tandem.dataframe(@path)
    p [17, df, df.shape]
    assert_equal df.shape, [1, 29]
  end
end