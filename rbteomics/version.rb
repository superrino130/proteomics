versioninfo = Struct.new('VersionInfo', 'major', 'minor', 'micro', 'releaselevel', 'serial')

class VersionInfo
  def initialize(...)
    __new__(...)
  end

  def __new__(cls, version_str)
    if version_str.instance_of?(String)
      @groups = /(\d+)\.(\d+)(?:\.)?(\d+)?([a-zA-Z]+)?(\d+)?/.match(version_str)
      @inst = super(cls).__new__(cls, *groups)
    else
      @inst = super(cls).__new__(cls, version_str.map{ _1.nil?.! ? x.to_s : x })
    end
    @inst._version_str = version_str
    @inst._version_ints = @inst.map{ _1.instance_of?(String) && /^[0-9]*$/.match(_1) ? _1.to_i : _1 }
    @inst
  end

  def __str__
    "Rbteomics version #{self._version_str}"
  end

  def __lt__(other)
    other = self.new(other) if other.instance_of?(VersionInfo).!
    self._version_inits < other._version_ints
  end

  def __gt__(other)
    other = self.new(other) if other.instance_of?(VersionInfo).!
    self._version_inits > other._version_ints
  end

  def __le__(other)
    self == other || self < other
  end

  def __ge__(other)
    self == other || self > other
  end

  def __eq__(other)
    other = self.new(other) if other.instance_of?(VersionInfo).!
    super(self).__eq__(other)
  end
end

# version_info = VersionInfo.new(__version__)
# version = __version__