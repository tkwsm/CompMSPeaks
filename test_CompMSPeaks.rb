#!/usr/bin/ruby

require 'test/unit'
require './CompMSPeaks.rb'

class ClassCompMSPeaks
  include CompMSPeaks
end

class TC_CompMSPeaks < Test::Unit::TestCase

  def setup1
    @cmsps  = ClassCompMSPeaks.new
  end

  def setup2
    sample_id = "0100"
    etype     = "p"
    mol_mass  = 200.0
    ret       = 18.5
    adducts    = [ "[M-H2O+H]+", "[M+H]+" ]
    @msp  = ClassCompMSPeaks::MSPeak.new( sample_id, etype, mol_mass, ret, adducts )
  end

  def setup3
    @msps = ClassCompMSPeaks::MSPeaks.new( 5, 100, 300, [File.new("./peak.test1.peak.table"), File.new("./peak.test2.peak.table")] ) 
  end

  def test_calc_ppm
    setup1
    assert_equal( ClassCompMSPeaks, @cmsps.class )
    assert_equal( 5.0e-4, @cmsps.calc_ppm( 100.0, 5) )
    assert_equal( 1.0e-3, @cmsps.calc_ppm( 100.0, 10) )
  end

  def test_mspeak
    setup2
    assert_equal( CompMSPeaks::MSPeak, @msp.class )
    assert_equal( "p", @msp.etype )
    assert_equal( [ "[M-H2O+H]+", "[M+H]+" ], @msp.adducts )
  end

  def test_mspeaks2
    setup1
    orig_a = [100.0000, 100.0001, 100.1005, 100.1006, 100.1007, 101.0001 ]
    cent_a = [100.0000,           100.1005,                     101.0001 ]
    assert_equal( cent_a, @cmsps.create_center_peaks( orig_a, [], 5 ) )
  end

  def test_mspeaks
    setup3
    assert_equal( CompMSPeaks::MSPeaks, @msps.class )
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

