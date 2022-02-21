require './rbteomics/auxiliary/file_helpers'
require 'set'
require 'bio'

protein = Struct.new(:protein, :description, :sequence)

module FASTABase
  def __init__(source, **kwargs)
    @_comments = '>;'.split('').to_set
    @_ignore_comments = kwargs.delete('ignore_comments') || false
    @parser = kwargs.delete('parser') || nil
    super
  end

  def _is_comment(line)
    @_comments.include?(line[0])
  end

  def get_entry(key)
    raise NotImplementedError.new
  end
end

class FASTA
  include FASTABase
  include FileReader

  def initialize(...)
    __init__(...)
  end

  def __init__(source, ignore_comments: false, parser: nil, encoding: nil)
    super(source, 'mode' => 'r', 'parser_func' => @_read, 'pass_file' => false, 'args' => [], 'kwargs' => {},
      'encoding' => encoding, 'ignore_comments' => ignore_comments, 'parser' => parser)
  end
end