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

class TestUseIndex < Minitest::Test
  class MockFile
    def initialize(seekable, mode)
      @mode = mode if mode.nil?.!
      @seekable = seekable if seekable.nil?.!
    end

    def seekable

    end
  end


  # Not started
end