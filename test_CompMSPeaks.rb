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
    etype     = "pos"
    mol_mass  = 200.0
    aducts    = [ "[M-H2O+H]+", "[M+H]+" ]
    @msp  = ClassCompMSPeaks::MSPeak.new( sample_id, etype, mol_mass, aducts )
  end

  def setup3
    @msps = ClassCompMSPeaks::MSPeaks.new( 5, File.open( "./peak.test1.table"), File.open( "./peak.test2.table" ) )
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
    assert_equal( "pos", @msp.etype )
    assert_equal( [ "[M-H2O+H]+", "[M+H]+" ], @msp.aducts )
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
    assert_equal( 2, @msps.mspeaks.keys.size )
    assert_equal( ["neg", "pos"], @msps.mspeaks.keys.sort )
    assert_equal( 169, @msps.mspeaks["pos"].size )
    assert_equal(  81, @msps.mspeaks["neg"].size )
    assert_equal( 304, @msps.mspeaks["pos"].values.flatten.size )
    assert_equal( 175, @msps.mspeaks["neg"].values.flatten.size )
    assert_equal( 304, @msps.pos_array.size )
    assert_equal( 175, @msps.neg_array.size )
    assert_equal(@msps.pos_array.size, @msps.mspeaks["pos"].values.flatten.size)
    assert_equal( [], @msps.diff_two_peaks( @msps.pos_array, @msps.mspeaks["pos"].values.flatten ) )
  end

end

=begin
=end
