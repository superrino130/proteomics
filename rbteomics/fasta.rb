require './rbteomics/auxiliary/file_helpers'
require 'set'

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

class FASTA < FileReader
  include FASTABase
end