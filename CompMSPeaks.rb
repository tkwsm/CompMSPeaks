#!/usr/local/bin/ruby

# = comp_peaks.rb - Comparing MSMS peak data for metabolome.
#
# Copyright::   Copyright (C) 2018 
#               Takeshi Kawashima
#
# $ID: $
#
# Comparing several MSMS peak data for metabolome.
#
# USAGE: ruby comp_peaks.rb <peak-table-1> <peak-table-2> ...
#
# Format of Peak-Table(s)
# Sample-ID  etype  mol_mass  aduct
# 
# type can take  "pos" or "neg".
# aduct (ion) can take multiple data separated by ", ".
# 
# Example of Peak-Table(s)
# 01111   pos     0.0     [M-2(H2O)+Na]+, [M-H2O+Na]+
# 01111   pos     0.0     [M-H2O+Na]+, [M+Na]+
# 01112   pos    12.0     [M-H2O+Na]+, [M+Na]+
# 01113   pos   330.0     [2M-H2O+Na]+, [M+Na]+


# Initialize Variables


###############################################
# Main Part : 0.
###############################################

module CompMSPeaks

  def calc_ppm( mol_mass, ppm )
    digit_ppm = 0.000001
    range     = mol_mass.to_f * digit_ppm * ppm.to_f
    return range
  end

  def diff_two_peaks( apeaks, bpeaks)
    ( apeaks - bpeaks )
  end

  def create_center_peaks( original_peaks, center_peaks, ppm )
    original_peaks.sort!
    original_peaks.each do |val|
      mol_mass = val
      range    = calc_ppm( mol_mass, ppm )
      if center_peaks.size == 0
        center_peaks << mol_mass
      else
        t_mol_mass = center_peaks[-1] 
        t_range    = calc_ppm( t_mol_mass, ppm )
        if t_mol_mass + t_range <= mol_mass 
          center_peaks << mol_mass
        end
      end
    end
    center_peaks.sort!
    return center_peaks
  end

  class MSPeak
    include CompMSPeaks

    def initialize( sid, etype, mol_mass, aducts )
      @sid = sid
      @etype    = etype
      @mol_mass = mol_mass
      @aducts   = aducts
    end

    attr_reader :sid, :etype, :mol_mass, :aducts

  end

  class MSPeaks
    include CompMSPeaks

    def initialize( ppm, *peak_table_files )
      @min_search_mass = 100.0
      @max_search_mass = 2000.0
      @ppm         = ppm
      @pos_array   = []
      @neg_array   = []
      @c_pos_peaks = []
      @c_neg_peaks = []
      @mspeaks = { "pos" => Hash.new, "neg" => Hash.new };
      add_table_files( peak_table_files )                   # @pos/@neg_array 
      create_center_peaks( @pos_array.collect{|item| item.mol_mass }, @c_pos_peaks, @ppm ) # @pos_center_peaks
      create_center_peaks( @neg_array.collect{|item| item.mol_mass }, @c_neg_peaks, @ppm ) # @neg_center_peaks
      init_mspeaks( "pos" )                   # @mspeaks["pos"]
      init_mspeaks( "neg" )                   # @mspeaks["neg"]
      fill_mspeaks( "pos" )
      fill_mspeaks( "neg" )
    end

    attr_reader :ppm, :pos_array, :neg_array, :c_pos_peaks, :c_neg_peaks,
                :mspeaks

    def fill_mspeaks( etype )
      center_peaks = @c_pos_peaks if etype == "pos"
      center_peaks = @c_neg_peaks if etype == "neg"
      orig_array   = @pos_array   if etype == "pos"
      orig_array   = @neg_array   if etype == "neg"
      orig_array.each do |a_mspeak|
        center_peaks.each do |mspeak|
          range = calc_ppm( mspeak, @ppm )
          if mspeak - range <= a_mspeak.mol_mass and 
             mspeak + range >  a_mspeak.mol_mass
            @mspeaks[ etype ][ mspeak ] << a_mspeak
            break
          end
        end
      end
    end

    def init_mspeaks( etype )
      center_peaks = @c_pos_peaks if etype == "pos"
      center_peaks = @c_neg_peaks if etype == "neg"
      center_peaks.each{ |masspeak| @mspeaks[ etype ][ masspeak ] = [] }
    end

    def add_original_peaks( sid, etype, mol_mass, aducts )
      @pos_array << MSPeak.new( sid, etype, mol_mass, aducts ) if etype == "pos"
      @neg_array << MSPeak.new( sid, etype, mol_mass, aducts ) if etype == "neg"
    end

    def add_table_data( peak_table_file )
      a = []
      open( peak_table_file ) .each do |x|
        a = x.chomp.split("\s")
        sid = a[0]
        etype    =  a[1]
        mol_mass =  a[2].to_f
        aducts   =  a[3..-1]
        next if mol_mass < @min_search_mass
        next if mol_mass > @max_search_mass
        add_original_peaks( sid, etype, mol_mass, aducts )
      end
    end

    def add_table_files( peak_table_files )
      peak_table_files.flatten.each{ |afile| add_table_data( File.expand_path(afile) ) }
    end

    def list_peak_cluster( etype )
      @mspeaks[etype].each_key do |mspeak|
        print @mspeaks[etype][mspeak].collect{|x| x.sid }.sort.uniq.join("\t"), "\n" if @mspeaks[etype][mspeak].size > 2
      end
    end

    def table_peak_cluster2( etype )
      h = {}
      @mspeaks[etype].each_key do |mspeak|
        p @mspeaks[etype][mspeak]
        p @mspeaks[etype][mspeak].collect{|x| x.sid }
        p @mspeaks[etype][mspeak].collect{|x| x.sid }.uniq 
      end
    end

    def table_peak_cluster( etype )
      h = {}
      @mspeaks[etype].each_key do |mspeak|
        if @mspeaks[etype][mspeak].size > 2
          k = @mspeaks[etype][mspeak].collect{|x| x.sid }.sort.uniq
          h[ k ] = 0 if h[ k ] == nil
          h[ k ] += 1
        end
      end
      flatten_keys = h.keys.flatten
      flatten_keys.sort!
      flatten_keys.uniq!

      print "num"; flatten_keys.each{|v| print "\t#{v}"}; print "\n"
      h_keys = h.keys.sort{|x, y| h[y] <=> h[x]}
      h_keys.each do |k|
        print "#{h[k]}"
        flatten_keys.each do |v|
          if k.include?( v )
            print "\t#{v}" 
          else
            print "\t" 
          end
        end
        print "\n"
      end

    end
  end

end

if $0 == __FILE__ 

  msps = CompMSPeaks::MSPeaks.new( 5, ARGV )
#  msps.list_peak_cluster( "pos" )
#  msps.table_peak_cluster( "pos" )
  msps.table_peak_cluster2( "neg" )

end

