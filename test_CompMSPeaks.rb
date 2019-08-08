#!/usr/bin/ruby

require 'test/unit'
require './CompMSPeaks.rb'

class ClassCompMSPeaks
  include CompMSPeaks
end

class TC_CompMSPeaks < Test::Unit::TestCase

  def setup3
    @msps = ClassCompMSPeaks::MSPeaks.new( 5, 100, 300, [File.new("./peak.test1.peak.table"), File.new("./peak.test2.peak.table")] ) 
  end

  def test_mspeaks
    setup3
#    assert_equal( CompMSPeaks::MSPeaks, @msps.class )
#    assert_equal( 2, @msps.mspeaks.keys.size )
#    assert_equal( [:n, :p], @msps.mspeaks.keys.sort )
#    assert_equal( 283, @msps.mspeaks[:p].size )
#    assert_equal( 152, @msps.mspeaks[:n].size )
#    assert_equal( 478, @msps.mspeaks[:p].values.flatten.size )
#    assert_equal( 266, @msps.mspeaks[:n].values.flatten.size )
#    assert_equal( 478, @msps.pos_array.size )
#    assert_equal( 266, @msps.neg_array.size )
#    assert_equal(@msps.pos_array.size, @msps.mspeaks[:p].values.flatten.size)
#    assert_equal( [], @msps.diff_two_peaks( @msps.pos_array, @msps.mspeaks[:p].values.flatten ) )
  end
=begin
=end

end

