require 'minitest/autorun'
require 'tempfile'
require 'stringio'
require_relative '../rbteomics/electrochem'
require_relative '../rbteomics/tandem'
require_relative '../rbteomics/version'
require_relative '../rbteomics/auxiliary/structures'
require_relative '../rbteomics/auxiliary/file_helpers'

require 'numpy'
require 'pandas'

PSMS = []
alph = ('A'..'Z').to_a + ('a'..'z').to_a
52.times do |i|
  PSMS << [i, alph[i], 0.01 + i * 0.001]
end

module Warning
  def self.warn(message, category: nil)
    @count << [category, message] if message.include?("$SAFE").!
    super("#{category} warning : #{message.chomp}\n")
  end

  def count
    @count
  end

  def reset_count
    @count = []
  end
end

class TestUseIndex < Minitest::Test
  class MockFile
    attr_reader :seekable, :mode
    def initialize(seekable, mode)
      @mode = mode if mode.nil?.!
      @seekable = seekable if seekable.nil?.!
    end

    def seekable

    end
  end


  def _check_file_object(fo, value)
    Warning.reset_count
    result = _check_use_index(fo, nil, nil)
    assert_equal Warning.count.size, 0
    assert_equal result, value
  end

  def test_str_name
    [false, true].each do |ui|
      [false, true].each do |default|
        assert_equal _check_use_index('test.mgf', ui, default), ui
      end
    end
  end

  def test_textfile
    f = './tests/test.fasta'
    _check_file_object(f, false)
  end

  # there is no binary in Ruby
  # def test_binfile
  #   f = './tests/test.mgf'
  #   # test.mgf is not binary file!
  #   # _check_file_object(f, true)
  #   _check_file_object(f, false)
  # end

  def test_tmpfile_text
    Tempfile.open("temp"){ |f|
      _check_file_object(f, false)
    }
  end

  # there is no binary in Ruby
  # def test_tmpfile_bin; end

  def test_stringio
    Warning.reset_count
    StringIO.open('test') do |f|
      _check_use_index(f, nil, nil)
      assert_equal Warning.count.size, 1
      assert_equal Warning.count[0][0], nil
    end
  end

  # there is no 'rb' in Ruby
  # def test_error_not_seekable
  #   source = MockFile.new(false, 'rb')
  #   Warning.reset_count
  #   _check_use_index(source, nil, nil)
  #   assert_equal Warning.count.size, 1
  #   assert_equal Warning.count[0][0], nil
  #   p [97,Warning.count[0]]
  # end

  def test_warning_not_seekable
    source = MockFile.new(false, 'r')
    Warning.reset_count
    _check_use_index(source, true, nil)
    assert_equal Warning.count.size, 1
    assert_equal Warning.count[0][0], nil

  end

  def test_warning_wrong_mode
    source = MockFile.new(true, 'r')
    Warning.reset_count
    _check_use_index(source, true, nil)
    assert_equal Warning.count.size, 1
    assert_equal Warning.count[0][0], nil
  end

  # def test_warning_no_mode
  #   source = MockFile.new(nil, nil)
  #   p [84, _check_use_index(source, nil, nil)]
  # end

  
end