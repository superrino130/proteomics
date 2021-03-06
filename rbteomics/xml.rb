# import re
# import socket
# from traceback import format_exc
# import warnings
# from collections import OrderedDict, namedtuple
# from itertools import islice
# from lxml import etree
# import numpy as np
require 'numpy'

# from .auxiliary import FileReader, PyteomicsError, basestring, _file_obj, HierarchicalOffsetIndex
# from .auxiliary import unitint, unitfloat, unitstr, cvstr
# from .auxiliary import _keepstate_method as _keepstate
# from .auxiliary import BinaryDataArrayTransformer
# from .auxiliary import TaskMappingMixin, IndexedReaderMixin, IndexSavingMixin
require_relative 'auxiliary/structures'
require_relative 'auxiliary/file_helpers'
require_relative 'auxiliary/utils'
BaseString ||= String


require 'set'
require 'open-uri'
require 'rexml/document'

module Xml
  module_function
  def _local_name(element)
    begin
      return element.local_name
    rescue => exception
      return ''
    end
  end
  
  def xsd_parser(schema_url)
    ret = Hash.new{ |h, k| h[k] = Set.new }
    if (schema_url.start_with?('http://') || schema_url.start_with?('https://') || schema_url.start_with?('file://')).!
      schema_url = 'file://' + schema_url
    end
    doc = ''
    begin
      URI.open(schema_url) do |schema_file|
        doc = REXML::Document.new schema_file
      end
    rescue => exception
      begin
        File.open('xml_schema/' + schema_url.split('/')[-1]) do |local_file|
          doc = REXML::Document.new local_file.read.gsub("\n", "").gsub("\t", "")
        end
      rescue => exception
        raize "Error URI and FILE"
      end
    end
    types = {
      'ints' => Set.new(['int', 'long', 'nonNegativeInteger', 'positiveInt', 'integer', 'unsignedInt']),
      'floats' => Set.new(['float', 'double']),
      'bools' => Set.new(['boolean']),
      'intlists' => Set.new(['listOfIntegers']),
      'floatlists' => Set.new(['listOfFloats']),
      'charlists' => Set.new(['listOfChars', 'listOfCharsOrAny'])
    }
    trees = Set.new
    doc.children.each do |ch|
      if ch.respond_to?(:expanded_name)
        trees << ch.expanded_name
      end
    end

    while trees.empty?.!
      adtrees = Set.new()
      trees.each do |tree|
        doc.elements.each(tree) do |chs|
          chs.children.each do |ch|
            if ch.respond_to?(:expanded_name)
              adtrees << tree + '/' + ch.expanded_name if ch.name != 'attribute'
              if ch.name == 'element'
                if ch.attribute('maxOccurs') && ch.attributes['maxOccurs'] != '1'
                  ret['lists'] << ch.attributes['name'] if ch.attributes['name'].nil?.!
                end
              elsif ch.name == 'attribute'
                  types.each do |k, val|
                  if ch.attribute('type') && val.include?(ch.attributes['type'].sub(ch.prefix + ':', ''))
                    case ch.parent.local_name
                    when 'element'
                      ret[k] << [ch.parent.attributes['name'], ch.attributes['name']]
                      break
                    when 'complexType'
                      if ch.parent.parent.local_name == 'element'
                        ret[k] << [ch.parent.parent.attributes['name'], ch.attributes['name']]
                        break
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
      trees = adtrees
    end
    ret
  end
  
  module XMLValueConverter
    module_function
    # @@_duration_parser = Regexp.compile('(?P<sign>-?)P(?:(?P<years>\d+\.?\d*)Y)?(?:(?P<months>\d+\.?\d*)M)?(?:(?P<days>\d+\.?\d*)D)?(?:T(?:(?P<hours>\d+\.?\d*)H)?(?:(?P<minutes>\d+\.?\d*)M)?(?:(?P<seconds>\d+\.?\d*)S)?)?')
    @@_duration_parser = Regexp.compile('(?<sign>-?)\D*(?<years>\d+\.?\d*)\D*(?<months>\d+\.?\d*)\D*(?<days>\d+\.?\d*)\D*(?<hours>\d+\.?\d*)\D*(?<minutes>\d+\.?\d*)\D*(?<seconds>\d+\.?\d*)')
    
    Duration_str_to_float = lambda do |s|
      if s.start_with?('P').!
        begin
          return UnitFloat.new(s, 'duration')
        rescue => exception
          return Unitstr.new(s, 'duration')
        end
      end
      match = @@_duration_parser.search(s)
      if match.nil?.!
        matchdict = match.groupdict()
        hours = matchdict['hours'].to_f || 0
        minutes = matchdict['minutes'].to_f || 0
        seconds = matchdict['seconds'].to_f || 0
        minutes += hours * 60.0
        minutes += seconds / 60.0
        return UnitFloat.new(minutes, 'minute')
      else
        return Unitstr.new(s, 'duration')
      end
    end
  
    Str_to_bool = lambda do |s|
      return true if ['true', '1', 'y'].include?(s.downcase)
      return false if ['false', '0', 'n'].include?(s.downcase)
      raise PyteomicsError.new('Cannot convert string to bool: ' + s)
    end
  
    Str_to_num = lambda do |s, numtype|
      ['', nil].include?(s).! ? numtype(s) : nil
    end
  
    def self.to(t)
      convert_from = lambda do |s|
        return str_to_num(s, t)
      end
      convert_from
    end
  
    def self.converters
      return {
        'ints' => self.to(Unitint),
        'floats' => self.to(UnitFloat),
        'bools' => Str_to_bool,
        'intlists' => lambda { |x| Numpy.fromstring(x.replace('\n', ' '), dtype: int, sep: ' ') },
        'floatlists' => lambda { |x| Numpy.fromstring(x.replace('\n', ' '), sep: ' ') },
        'charlists' => [],
        'duration' => Duration_str_to_float
      }
    end
  end
  
  SXMLParam ||= Struct.new(:name, :value, :type)
  
  module XMLParam
    module_function
    def __init__(sxmlparam: nil)
      @sxmlparam = sxmlparam
    end
  end

  module XML
    include File_helpers
    include File_helpers::FileReader
    module_function
    
    def _get_info_smart(element, **kwargs)
      raise NotImplementedError
    end
  
    def __init__(source, **kwargs)
      @file_format ||= 'XML'
      @_root_element ||= nil
      @_default_schema ||= {}
      @_read_schema ||= false
      @_default_version ||= 0
      @_default_iter_tag ||= nil
      @_default_iter_path ||= nil
      @_structures_to_flatten ||= []
      @_schema_location_param ||= 'schemaLocation'
      @_default_id_attr ||= 'id'
      @_huge_tree ||= false
      @_retrieve_refs_enabled ||= nil  # only some subclasses implement this
      @_iterative ||= true
  
      # Configurable plugin logic
      @_converters ||= XMLValueConverter.converters()
      @_element_handlers ||= {}
  
      read_schema = kwargs['read_schema'] || nil
      iterative = kwargs['iterative'] || nil
      build_id_cache = kwargs['build_id_cache'] || false
  
      super(source, 'mode' => 'rb', 'parser_func' => Set_Iterfind, 'pass_file' => false,
        'args' => (@_default_iter_path || @_default_iter_tag), 'kwargs' => kwargs)
      iterative = @_iterative if iterative.nil?
      if iterative.nil?.!
        @_tree ||= nil
      else
        build_tree()
      end
      if build_id_cache
        build_id_cache()
      else
        @_id_dict ||= nil
      end
  
      @version_info ||= _get_version_info()
      if read_schema.nil?.!
        @_read_schema ||= read_schema
      end
      @schema_info ||= _get_schema_info(read_schema: read_schema)
  
      @_converters_items ||= @_converters
      @_huge_tree ||= kwargs['huge_tree'] || @_huge_tree
      @_retrieve_refs_enabled ||= kwargs['retrieve_refs']
    end
  
    def __reduce_ex__(protocol)
      [self.class, [@_source_init, @_read_schema, @_tree.nil?, false], __getstate__()]
    end
  
    def __getstate__
      state = super
      state['_huge_tree'] = @_huge_tree
      state['_retrieve_refs_enabled'] = @_retrieve_refs_enabled
      state['_id_dict'] = @_id_dict
      state
    end
  
    def __setstate__(state)
      super
      @_huge_tree = state['_huge_tree']
      @_retrieve_refs_enabled = state['_retrieve_refs_enabled']
      @_id_dict = state['_id_dict']
    end
  
    #@_keepstate
    def in_get_version_info
      version = '1.0'
      schemaLocation = ''
      @_source.each_with_index do |line, i|
        next if i == 0
        if line.include?('version')
          version = line.match(/version=\"([0-9]|\.)*\"/)[0].sub('version=', '').gsub('"', '')
        elsif line.include?(':schemaLocation')
          schemaLocation = line.match(/(http|ftp).+xsd/)[0]
          break
        end
      end
      [version, schemaLocation]

      # mzIdentML
      # return ['1.1.0', 'https://psidev.info/sites/default/files/2017-06/mzIdentML1.1.0.xsd']

      # pepxml
      # return ['1.0', 'http://regis-web.systemsbiology.net/pepXML http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v115.xsd']
    end
    def _get_version_info
      _keepstate_method(:in_get_version_info)
      in_get_version_info
    end
  
    #@_keepstate
    def in_get_schema_info(read_schema: true)
      return @_default_schema if read_schema.!
  
      version, schema = @version_info
      return @_default_schema if version == @_default_version
  
      ret = {}
      begin
        if ['', 0, nil, false, [], {}].include?(schema)
          schema_url = ''
          raise PyteomicsError.new("Schema information not found in #{@name}.")
        end
        schema_url = schema.split[-1]
        ret = xsd_parser(schema_url)
      rescue => e
        # if [OpenURI::HTTPError, socket.error, socket.timeout].include?(e.class)
        if e.instance_of?(OpenURI::HTTPError)
            warn "Can't get the #{self.class} schema for version `#{@version}` from <#{@schema_url}> at the moment.\nUsing defaults for {0._default_version}.\nYou can disable reading the schema by specifying `read_schema=False`."
        else
          warn "Unknown #{self.class} version `#{@version}`.\n" +
            "Attempt to use schema " +
            "information from <#{schema_url}> failed.\n"
        end
        ret = @_default_schema
      end
      ret
    end
    def _get_schema_info(read_schema: true)
      _keepstate_method(:in_get_schema_info)
      in_get_schema_info(read_schema: read_schema)
    end
  
    def _handle_param(element, **kwargs)
      types = {'int' => unitint, 'float' => unitfloat, 'string' => unitstr}
      attribs = element.attrib
      unit_info = nil
      unit_accesssion = nil
      if attribs.include?('unitCvRef') || attribs.include?('unitName')
        unit_accesssion = attribs['unitAccession']
        unit_name = attribs['unitName'] || unit_accesssion
        unit_info = unit_name
      end
      accession = attribs['accession']
      value = attribs['value'] || ''
      begin
        if types.include?(attribs['type'])
          value = types[attribs['type']].call(value, unit_info)
        else
          value = unitfloat(value, unit_info)
        end
      rescue => exception
        value = unitstr(value, unit_info)
      end
      return XMLParam.new(Cvstr.new(attribs['name'], accession, unit_accesssion), value, _local_name(element))
    end
  
    def _find_immediate_params(element, **kwargs)
      return element.xpath('./*[local-name()="cvParam" or local-name()="userParam" or local-name()="UserParam"]')
    end
  
    def _insert_param(info_dict, param)
      key = param.name
      if info_dict.include?(key)
        if info_dict[key].is_a?(Array)
          info_dict[key] << param.value
        else
          info_dict[key] = [info_dict[key], param.value]
        end
      else
        info_dict[key] = param.value
      end
    end
  
    Promote_empty_parameter_to_name = lambda do |info, params|
      empty_values = []
      not_empty_values = []
      params.each do |param|
        if param.is_empty()
          empty_values << param
        else
          not_empty_values << param
        end
      end
  
      if empty_values.size == 1 && info.include?('name').!
        info['name'] = empty_values[0].name
        return [info, not_empty_values]
      end
      [info, params]
    end
  
    def _get_info(element, **kwargs)
      begin
        name = kwargs.delete('ename')
      rescue => exception
        name = _local_name(element)
      end
      schema_info = @schema_info
      if ['cvParam', 'userParam', 'UserParam'].include?(name)
        return _handle_param(element, **kwargs)
      end
  
      info = element.attrib.to_h
      params = []
      if kwargs['recursive']
        element.iterchildren().each do |child|
          cname = _local_name(child)
          if ['cvParam', 'userParam', 'UserParam'].include?(cname)
            newinfo = _handle_param(child, **kwargs)
            params << newinfo
          else
            if schema_info['lists'].include?(cname).!
              info[cname] = self._get_info_smart(child, ename=cname, **kwargs)
            else
              info[cname] = [] if info.include?(cname).!
              info[cname] << _get_info_smart(child, 'ename' => cname, **kwargs)
            end
          end
        end
      else
        _find_immediate_params(element, **kwargs).each do |child|
          params << _handle_param(child, **kwargs)
        end
      end
  
      handler = @_element_handlers[name]
      if handler.nil?.!
        info, params = handler(info, params)
      end
  
      params.each do |param|
        _insert_param(info, param)
      end
  
      if element.text
        stext = element.text.strip
        if stext
          if info
            info[name] = stext
          else
            return stext
          end
        end
      end
  
      begin
        info.each do |k, v|
          @_converters_items.each do |t, a|
            if schema_info.include?(t) && schema_info[t].include?([name, k])
              info[k] = a(v)
            end
          end
        end
      rescue => e
        message = "Error when converting types: #{e.args}"
        if _read_schema.!
          message += '\nTry reading the file with read_schema=True'
        end
        raise PyteomicsError.new(message)
      end
  
      if kwargs['retrieve_refs'] || @_retrieve_refs_enabled
        _retrieve_refs(info, **kwargs)
      end
  
      info.to_h.each do |k, v|
        if @_structures_to_flatten.include?(k)
          if v.is_a?(Array)
            v.each do |vi|
              info.merge!(vi)
            end
          else
            info.merge!(v)
          end
          info.delete(k)
        end
      end
  
      info.to_h.each do |k, v|
        if v.is_a?(Hash) && v.include?('name') && v.size == 1
          info[k] = v['name']
        end
      end
      if info.size == 2 && info.include?('name') && (info.include?('value') || info.include?('values'))
        name = info.delete('name')
        info = {'name' => info.pop[1]}
      end
      info
    end
  
    #@_keepstate
    def in_build_tree
      ps = etree.XMLParser(remove_comments: true, huge_tree: true)
      @_tree = etree.parse(@_source, 'parser' => ps)
    end
    def build_tree
      _keepstate_method(:in_build_tree)
      in_build_tree
    end
  
    def clear_tree
      @_tree = nil
    end
  
    def _retrieve_refs(info, **kwargs)
      raise NotImplementedError("_retrieve_refs is not implemented for #{self.class}. Do not use `retrieve_refs=True`.")
    end
  
    Set_Iterfind = lambda do |path, **kwargs|
      Iterfind.__init__(path, **kwargs)
    end
  
    #@_keepstate
    def in_iterfind_impl(path, **kwargs)
      begin
        m = pattern_path.match(path)
        path = m[1]
        tail = m[0].sub(path,'')
      rescue => exception
        raise PyteomicsError.new('Invalid path: ' + path)
      end
      if path[0...2] == '//' || path[0] != '/'
        absolute = false
        if path[0...2] == '//'
          path = path[2..]
          if path[0] == '/' || path.include?('//')
            raise PyteomicsError.new("Too many /'s in a row.")
          end
        end
      else
        absolute = true
        path = path[1..]
      end
      nodes = path.rstrip('/').split('/')
      if nodes.!
        raise PyteomicsError.new('Invalid path: ' + path)
      end
  
      if @_tree.!
        if tail
          if tail[0] == '['
            tail = '(.)' + tail
          else
            raise PyteomicsError.new('Cannot parse path tail: ' + tail)
          end
          xpath = etree.XPath(tail)
        end
        localname = nodes[0]
        found = false
        etree.iterparse('events' => ['start', 'end'], 'remove_comments' => true, 'huge_tree' => @_huge_tree).each do |ev, elem|
          name_lc = _local_name(elem)
          if ev == 'start'
            if name_lc == localname || localname == '*'
              found += 1
            end
          else
            if name_lc == localname || localname == '*'
              if (absolute && elem.getparent().nil?) || absolute.!
                get_rel_path(elem, nodes[1..]).each do |child|
                  if tail
                    xpath(child).each do |elem|
                      info = _get_info_smart(elem, **kwargs)
                      yield info
                    end
                  else
                    info = _get_info_smart(child, **kwargs)
                    yield info
                  end
                end
              end
              if localname != '*'
                found -= 1
              end
            end
            if found.!
              elem.clear
            end
          end
        end
      else
        xpath = (absolute ? '/' : '//') + nodes.map{ |node| node != '*' ? "*[local-name()='#{node}']" : '*' }.join('/') + tail
        @_tree.xpath(xpath).each do |elem|
          info = _get_info_smart(elem, **kwargs)
          yield info
        end
      end
    end
    def _iterfind_impl(path, **kwargs)
      _keepstate_method(:in_iterfind_impl)
      in_iterfind_impl
    end
  
    #@_keepstate
    def in_build_id_cache
      stack = 0
      id_dict = {}
      etree.iterparse(@_source, 'events' => ['start', 'end'],
        'remove_comments' => true, 'huge_tree' => @_huge_tree).each do |event, elem|
        if event == 'start'
          if elem.attrib.include?('id')
            stack += 1
          end
        else
          if elem.attrib.include?('id')
            stack -= 1
            id_dict[elem.attrib['id']] = elem
          elsif stack == 0
            elem.clear
          end
        end
      end
      @_id_dict = id_dict
    end
    def build_id_cache
      _keepstate_method(:in_build_id_cache)
      in_build_id_cache
    end
  
    def clear_id_cache
      @_id_dict = {}
    end
  
    def _find_by_id_no_reset(elem_id, id_key: nil)
      found = false
      id_key = @_default_id_attr if id_key.nil?
      etree.iterparse(@_source, 'events' => ['start', 'end'], 'remove_comments' => true, 'huge_tree' => @_huge_tree).each do |event, elem|
        if event == 'start'
          if elem.attrib[id_key] == elem_id
            found = true
          end
        else
          return elem if elem.attrib[id_key] == elem_id
          elem.clear if found.!
        end
      end
      raise KeyError elem_id
    end
  
    #@_keepstate
    def in_get_by_id(elem_id, **kwargs)
      if @_id_dict.!
        elem = _find_by_id_no_reset(elem_id)
      else
        elem = @_id_dict[elem_id]
      end
      _get_info_smart(elem, **kwargs)
    end
    def get_by_id(elem_id, **kwargs)
      _keepstate_method(:in_get_by_id)
      in_get_by_id(elem_id, **kwargs)
    end

    def self._element_handlers
      @_element_handlers || {}
    end

    def version_info
      # PASS
    end
  end
  
  Pattern_path = Regexp.compile('([\w/*]*)(.*)')
  
  def get_rel_path(element, names)
    # Not started
  end
  
  def xpath(tree, path, ns: nil)
    # Not started
  end
  
  def _make_version_info(cls)
    version_info = lambda do |source|
      cls.new(source).version_info
    end
    # version_info.__doc__ = "
    # Provide version information about the #{@file_format} file.
  
    # .. note:: This function is provided for backward compatibility only.
    #     It simply creates an :py:class:`#{cls.class}` instance
    #     and returns its :py:data:`!version_info` attribute.
  
    # Parameters
    # ----------
    # source : str or file
    #     File name or file-like object.
  
    # Returns
    # -------
    # out : tuple
    #     A (version, schema URL) tuple, both elements are strings or None.
    # "
    version_info.call('Not started')
  end
  
  class ByteCountingXMLScanner < File_helpers::File_obj
    include File_helpers
    include Xml

    @@entities = {
      'quot' => '"',
      'amp' => '&',
      'apos' => "'",
      'lt' => '<',
      'gt' => '>'
    }
    @@xml_entity_pattern = Regexp.compile("&(#{@@entities.keys.join('|')});")

    def initialize(...)
      __init__(...)
    end

    def __init__(source, indexed_tags, block_size: 1000000)
      super(source, 'rb')
      @indexed_tags = ensure_bytes(indexed_tags)
      # @indexed_tags = Xml.ensure_bytes(indexed_tags)
      @block_size = block_size
    end

    def _chunk_iterator
      # n s
    end

    def _generate_offsets
      # n s
    end

    _entity_sub_cb = lambda do |match|
      # n s
    end

    def replace_entities(key)
      return @@xml_entity_pattern.sub(_entity_sub_cb, key)
    end

    def in_build_byte_index(lookup_id_key_mapping: nil)
      if lookup_id_key_mapping.nil?
        lookup_id_key_mapping = {}
      end
      lookup_id_key_mapping = lookup_id_key_mapping.map{ |key, value| [ensure_bytes_single(key), ensure_bytes_single(value)] }

      @indexed_tags.each do |name|
        bname = ensure_bytes_single(name)
        lookup_id_key_mapping[bname] = 'id' if lookup_id_key_mapping.include?(bname).!
        lookup_id_key_mapping[bname] = ensure_bytes_single(lookup_id_key_mapping[bname])
      end

      indices = HierarchicalOffsetIndex.new
      g = _generate_offsets()
      g.each do |offset, offset_type, attrs|
        k = attrs[lookup_id_key_mapping[offset_type]].decode('utf-8')
        if k.include?('&')
          k = replace_entities(k)
        end
        indices[offset_type.decode('utf-8')][k] = offset
      end
      return indices
    end
    def build_byte_index(lookup_id_key_mapping: nil)
      File_helpers._keepstate_method(:in_build_byte_index)
      in_build_byte_index(lookup_id_key_mapping: lookup_id_key_mapping)
    end

    def in_scan(cls, source, indexed_tags)
      inst = cls.new(source, indexed_tags)
      return inst.build_byte_index()
    end
    def scan(cls, source, indexed_tags)
      _keepstate_method(:in_scan)
      in_scan(cls, source, indexed_tags)
    end
  end
  
  class TagSpecificXMLByteIndex
    def initialize(...)
      @_default_indexed_tags = []
      @_default_keys = {}
      @_scanner_class = ByteCountingXMLScanner
      __init__(...)
    end
  
    def __init__(source, indexed_tags: nil, keys: nil)
      keys = @_default_keys.dup if keys.nil?
      indexed_tags = @_default_indexed_tags if indexed_tags.nil?
      @indexed_tags = indexed_tags
      @indexed_tag_keys = keys
      @source = source
      @offsets = File_helpers::HierarchicalOffsetIndex.__init__(base: nil)
      build_index()
    end
  
    def __getstate__
      # Not started
    end
  
    def __setstate__(state)
      # Not started
    end
  
    def __getitem__(key)
      @offsets[key]
    end
  
    def build_index
      scanner = @_scanner_class.new(@source, @indexed_tags)
      @offsets = scanner.build_byte_index(@indexed_tag_keys)
      @offsets
    end
  
    def items
      @offsets
    end
  
    def keys
      @offsets.keys
    end
  
    def __iter__
      @keys.to_enum
    end
  
    def __len__
      items().sum{ |_, group| group.size }
    end
  
    def self.build(source, indexed_tags: nil, keys: nil)
      indexer = self.new(source, indexed_tags: indexed_tags, keys: keys)
      indexer.offsets
    end
  end
  
  def ensure_bytes_single(string)
    begin
      return string.encode(Encoding::UTF_8)
    rescue => exception
      raise PyteomicsError.new("#{string.inspect} could not be encoded")
    end
  end
  
  def ensure_bytes(strings)
    if strings.instance_of?(BaseString)
      strings = [strings]
    end
    return strings.map{ |string| ensure_bytes_single(string) }

  end
  
  def _flatten_map(hierarchical_map)
    all_records = []
    hierarchical_map.each do |key, records|
      all_records.concat(records)
    end
    all_records.sort!{ |x| x[1] }
    OrderedDict(all_records)
  end
  
  module IndexedXML
    extend File_helpers
    extend File_helpers::IndexSavingMixin
    include XML
    include File_helpers::IndexedReaderMixin

    module_function
  
    def __init__(source, *args, **kwargs)
      @_indexed_tags ||= Set.new
      @_indexed_tag_keys ||= {}
      @_use_index ||= true
  
      kwargs['read_schema'] ||= false
      kwargs['iterative'] ||= true
      kwargs['build_id_cache'] ||= false
      build_id_cache = kwargs['build_id_cache']
      use_index = kwargs['use_index'] || nil
      tags = kwargs['indexed_tags']
      tag_index_keys = kwargs['indexed_tag_keys']
  
      @_indexed_tags = tags if tags.nil?.!
      @_indexed_tag_keys = tag_index_keys if tag_index_keys.nil?.!
  
      @_use_index = use_index if use_index.nil?.!
      if ['', 0, nil, false, [], {}].include?(use_index).!
        build_id_cache = false
        if ['', 0, nil, false, [], {}].include?(@_default_iter_path).! && @_default_iter_path != @_default_iter_tag
          warn('_default_iter_path differs from _default_iter_tag and index is enabled. _default_iter_tag will be used in the index, mind the consequences.')
        end
      end
      super(source, **kwargs)
  
      @_offset_index = nil
      self._build_index()
    end
  
    # @property
    def default_index
      @_offset_index[@_default_iter_tag]
    end
  
    def __reduce_ex__(protocol)
      reconstructor, args, state = XML.new.__reduce_ex__(protocol)
      args += [False]
      [reconstructor, args, state]
    end
  
    def __getstate__
      @k.each do |obj|
        if obj.respond_to?(:__getstate__)
          state = obj.__getstate__
          break
        end
      end
      state['_indexed_tags'] = @_indexed_tags
      state['_indexed_tag_keys'] = @_indexed_tag_keys
      state['_use_index'] = @_use_index
      state['_offset_index'] = @_offset_index
      state
    end
  
    def __setstate__(state)
      @k.each do |obj|
        if obj.respond_to?(:__setstate__)
          obj.__setstate__
          break
        end
      end
      @_indexed_tags = state['_indexed_tags']
      @_indexed_tag_keys = state['_indexed_tag_keys']
      @_use_index = state['_use_index']
      @_offset_index = state['_offset_index']
    end
  
    # @_keepstate
    def in_build_index
      if ['', 0, nil, false, [], {}].include?(@_indexed_tags) || ['', 0, nil, false, [], {}].include?(@_use_index)
        return
      end
      @_offset_index = TagSpecificXMLByteIndex.build(
        @_source, indexed_tags: @_indexed_tags, keys: @_indexed_tag_keys)
    end
    def self._build_index
      _keepstate_method(:in_build_index)
      in_build_index
    end
  
    #@_keepstate
    def in_find_by_id_reset(elem_id, id_key: nil)
      _find_by_id_no_reset(elem_id, id_key: id_key)
    end
    def _find_by_id_reset(elem_id, id_key: nil)
      _keepstate_method(:in_find_by_id_reset)
      in_find_by_id_reset(elem_idk id_key: id_key)
    end
  
    #@_keepstate
    def in_get_by_id(elem_id, **kwargs)
      id_key = kwargs['id_key'] || nil
      element_type = kwargs['element_type'] || nil
      begin
        index = @_offset_index
        if element_type.nil?
          offset, element_type = index.find_no_type(elem_id)
        else
          offset = index.find(elem_id, element_type)
        end
        @_source.seek(offset)
        if id_key.nil?
          id_key = @_indexed_tag_keys[element_type]
        end
        elem = _find_by_id_no_reset(elem_id, id_key: id_key)        
      rescue => exception
        elem = _find_by_id_reset(elem_id, id_key: id_key)      
      end
      data = _get_info_smart(elem, **kwargs)
    end
    def get_by_id(elem_id, **kwargs)
      _keepstate_method(:in_get_by_id)
      in_get_by_id(elem_id, **kwargs)
    end
  
    def __contains__(key)
      @_offset_index[@_default_iter_tag].include?(key)
    end
  
    def size
      @_offset_index[@_default_iter_tag].size
    end
  
    def iterfind(path, **kwargs)
      if @_indexed_tags.nil?.! && @_indexed_tags.include?(path) && @_use_index
        return IndexedIterfind.new(path, **kwargs)
      end
      Iterfind.new(path, **kwargs)
    end
  end
  
  module MultiProcessingXML
    include File_helpers::TaskMappingMixin
    include IndexedXML
    module_function
  
    def _task_map_iterator
      @_offset_index[self._default_iter_tag].to_enum
    end  
  end
  
  module IndexSavingXML
    include Xml
    include IndexedXML
    include File_helpers::IndexSavingMixin
    module_function
    
    def __init__(...)
      if @_index_class.nil? || @_index_class.empty?
        @_index_class = File_helpers::HierarchicalOffsetIndex.__init__
      end
      super
    end
    
    def _read_byte_offsets
      File.open(_byte_offset_filename, 'w') do |f|
        index = @_index_class.load(f)
        if index.schema_version.nil?
          raise TypeError("Legacy Offset Index!")
        end
        @_offset_index = index
      end
    end
  end
  
  class ArrayConversionMixin < BinaryDataArrayTransformer
    @@_dtype_dict = {}
    @@_array_keys = ['m/z array', 'intensity array']
    
    def initialize(...)
      __init__(...)
    end
  
    def __init__(*args, **kwargs)
      @@_dtype_dict = {'None' => nil}
      dtype = kwargs.delete('dtype') || nil
      if dtype.is_a?(Hash)
        @@_dtype_dict.merge!(dtype)
      elsif dtype.nil?.!
        @@_dtype_dict = @@_array_keys.map{ |k| [k, dtype] }.to_h
        @@_dtype_dict['None'] = dtype
      end
      # super
    end
  
    def __getstate__
      state = super
      state['_dtype_dict'] = @@_dtype_dict
      state
    end
  
    def __setstate__(state)
      super
      @@_dtype_dict = state['_dtype_dict']
    end
  
    def _convert_array(k, array)
      dtype = @@_dtype_dict[k]
      if dtype.nil?.!
        return array.astype(dtype)
      end
      array
    end
  
    def _finalize_record_conversion(array, record)
      key = record.key
      return _convert_array(key, array)
    end
  end
  
  module Iterfind
    @@md = Struct.new(:start, :stop, :step)
    module_function
  
    def __init__(parser, tag_name, **kwargs)
      @parser ||= parser
      @tag_name ||= tag_name
      @config ||= kwargs
      @_iterator ||= nil
    end
  
    def __repr__
      if @config != ''
        @config = ", " + repr(@config)
      else
        @config = ''
      end
      template = "#{self.class}(#{__repr__(@tag_name)}#{@config})"
    end
  
    def __iter__
      self
    end
  
    def _make_iterator
      @parser._iterfind_impl(@tag_name, **@config)
    end
  
    def __next__
      if @_iterator.nil?
        @_iterator = _make_iterator()
      end
      return _next(@_iterator)
    end
  
    def _next
      __next__()
    end
  
    #@property
    def is_indexed
      false
    end
  
    def reset
      @_iterator = nil
      @parser.reset()
    end
  
    def __enter__
      self
    end
  
    def __exit__(*args, **kwargs)
      reset()
    end
  
    def map(*args,**kwargs)
      raise NotImplementedError.new("This query isn't indexed, it cannot be mapped with multiprocessing")
    end
  
    def _get_by_index(idx)
      self.reset()
      value = _next(enum_for(:islice, self, idx, idx + 1))
    end
  
    def _get_by_slice(*slc)
      self.reset()
      value = enum_for(:islice, self, *slc)
    end
  
    def __getitem__(i)
      if i.is_a?(Struct)
          return _get_by_slice(i)
      end
      return _get_by_index(i)
    end
  end
  
  class IndexedIterfind
    include Iterfind
    include File_helpers::TaskMappingMixin

    def _task_map_iterator
      _index.to_enum
    end
  
    def _index
      @parser.index(@tag_name)
    end
  
    def _get_reader_for_worker_spec
      @parser
    end
  
    def _yield_from_index
      _task_map_iterator.each do |key, _|
        yield @parser.get_by_id(key, @config)
      end
    end
  
    def _make_iterator
      if is_indexed
        return _yield_from_index()
      end
      warn "Non-indexed iterator created from #{self}"
      super
    end
  
    #@property
    def is_indexed
      if @parser.respond_to?(:_index)
        if @parser._index.nil?.!
          index = @parser._index
          if index.is_a?(HierarchicalOffsetIndex)
            return index.include?(@tag_name) && ['', 0, nil, false, [], {}].include?(index[@tag_name]).!
          end
        end
      end
      false
    end
  
    def _get_by_index(idx)
      key = @index.form_index(idx)
      @parser.get_by_id(key)
    end
  
    def _get_by_slice(slc)
      keys = @index.form_slice(slc)
      @parser.get_by_ids(keys)
    end
  
    def __len__
      @index.size
    end  
  end
end
