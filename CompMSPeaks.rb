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
# Sample-ID  etype  mol_mass  ret aduct
# 
# type can take  "pos" or "neg".
# aduct (ion) can take multiple data separated by ", ".
# 
# Example of Peak-Table(s)
# 01111   pos     0.0     18.1	[M-2(H2O)+Na]+, [M-H2O+Na]+
# 01111   pos     0.0     10.1	[M-H2O+Na]+, [M+Na]+
# 01112   pos    12.0     11.1	[M-H2O+Na]+, [M+Na]+
# 01113   pos   330.0     15.2	[2M-H2O+Na]+, [M+Na]+

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

    def initialize( sid, etype, mol_mass, ret, adducts )
      @sid = sid
      @etype    = etype
      @mol_mass = mol_mass
      @ret      = ret
      @adducts   = adducts
    end

    attr_reader :sid, :etype, :mol_mass, :ret, :adducts

  end

  class MSPeaks
    include CompMSPeaks

    def initialize( ppm, min_search_mass, max_search_mass, peak_table_files )
      @min_search_mass = min_search_mass.to_i
      @max_search_mass = max_search_mass.to_i
      @ppm           = ppm
      @pos_array     = []
      @neg_array     = []
      @c_pos_peaks   = []
      @c_neg_peaks   = []
      @samples       = {}
      @samples_order = []
      @mspeaks = { :p => Hash.new, :n => Hash.new };
      STDERR.puts "The loading files:  #{peak_table_files.join(":")}"
      add_table_files( peak_table_files )                   # @pos/@neg_array 
      STDERR.puts "pos stored peaks total #{@pos_array.size}"
      STDERR.puts "neg stored peaks total #{@neg_array.size}"
      create_center_peaks( @pos_array.collect{|item| item.mol_mass }, @c_pos_peaks, @ppm ) # @pos_center_peaks
      STDERR.puts "pos-peak keys size #{@c_pos_peaks.size}"
      create_center_peaks( @neg_array.collect{|item| item.mol_mass }, @c_neg_peaks, @ppm ) # @neg_center_peaks
      STDERR.puts "neg-peak keys size #{@c_neg_peaks.size}"
      init_mspeaks( :p )                   # @mspeaks[ :p ]
      init_mspeaks( :n )                   # @mspeaks[ :n ]
      fill_mspeaks( :p )
      fill_mspeaks( :n )
    end

    attr_reader :ppm, :pos_array, :neg_array, :c_pos_peaks, :c_neg_peaks,
                :mspeaks

    def fill_mspeaks( etype )
      center_peaks = @c_pos_peaks if etype == :p
      center_peaks = @c_neg_peaks if etype == :n
      orig_array   = @pos_array   if etype == :p
      orig_array   = @neg_array   if etype == :n
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
      center_peaks = @c_pos_peaks if etype == :p
      center_peaks = @c_neg_peaks if etype == :n
      center_peaks.each{ |masspeak| @mspeaks[ etype ][ masspeak ] = [] }
    end

    def add_original_peaks( sid, etype, mol_mass, ret, adducts )
      @pos_array << MSPeak.new(sid, etype, mol_mass, ret, adducts) if etype== :p
      @neg_array << MSPeak.new(sid, etype, mol_mass, ret, adducts) if etype== :n
    end

    def add_table_data( peak_table_file )
      a = []
      open( peak_table_file ) .each do |x|
        a = x.chomp.split("\s")
        sid = a[0]
        etype    =  :p if a[1] == "pos"
        etype    =  :n if a[1] == "neg"
        mol_mass =  a[2].to_f
        ret      =  a[3].to_f
        adducts   =  a[4..-1]
        @samples[ sid ] = ""
        @samples_order = @samples.keys.sort
        next if mol_mass < @min_search_mass
        next if mol_mass > @max_search_mass
        add_original_peaks( sid, etype, mol_mass, ret, adducts )
      end
    end

    def add_table_files( peak_table_files )
      peak_table_files.flatten.each{ |afile| add_table_data( File.expand_path(afile) ) }
    end

    def list_sid_by_peak_cluster( etype )
      @mspeaks[etype].each_key do |mspeak|
        print @mspeaks[etype][mspeak].collect{|x| x.sid }.sort.uniq.join("\t"), "\n" if @mspeaks[etype][mspeak].size > 2
      end
    end

    def table_spotcounts_by_peak_cluster( etype )
      print "mspeak"
      @samples_order.each{|sid| print "\t", sid }
      print "\n"
      @mspeaks[ etype ].each_key do |mspeak|
        print "#{mspeak}"
        tmp_samples = @mspeaks[ etype ][ mspeak ].collect{|x| x.sid } 
        @samples_order.each do |sid|
           print "\t", tmp_samples.count( sid )
        end
         print "\n"
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
        k = @mspeaks[etype][mspeak].collect{|x| x.sid }.sort.uniq
        if k.size > 1
          h[ k ] = 0 if h[ k ] == nil
          h[ k ] += 1
        end
      end

      print "num"; @samples_order.each{|v| print "\t#{v}"}; print "\n"
      h_keys = h.keys.sort{|x, y| h[y] <=> h[x]}
      h_keys.each do |k|
        print "#{h[k]}"
        @samples_order.each do |v|
          if k.include?( v )
            print "\t#{v}" 
          else
            print "\t" 
          end
        end
        print "\n"
      end

    end

    def list_clade_specific_peaks( etype, sid_list )
      h = {}
      @mspeaks[etype].each_key do |mspeak|
        k = @mspeaks[etype][mspeak].collect{|x| x.sid }.sort.uniq
        next if k != sid_list.sort
        next if @mspeaks[etype][mspeak].size < k.size
        @mspeaks[etype][mspeak].each do |apeak|
          print "#{etype}\t#{mspeak}\t#{apeak.mol_mass}\t#{apeak.ret}\n"
        end
      end
    end

  end


end

if $0 == __FILE__ 

  if ARGV[0] !~ /\d/ or ARGV[1] !~ /\d/ or ARGV[0] =~ /\D/
    STDERR.print "Please check the USAGE\n or \n"
    STDERR.print "ruby CompMSPeaks.rb 100, 2000, <peaktable1> ...\n\n"
    exit
  end

  STDERR.puts "START"

  min_peak       = ARGV.shift
  max_peak       = ARGV.shift
  sample_id_list = ARGV.shift
  peak_table_dir = ARGV.shift

  check_list_hash = {}
  File.open( sample_id_list ).each do |afile_line|
    STDERR.puts afile_line
    file_name = afile_line.chomp.split("\t")[0]
    STDERR.puts file_name
    check_list_hash[ file_name ] = ""
  end
  
  STDERR.puts  check_list_hash
  STDERR.puts "READFILES"
  peak_tables    = []
  Dir.open( peak_table_dir ).each do |afile|

    next if afile =~ /^\./
    next if afile !~ /\S+.peak.table/
    STDERR.puts afile
    file_name = afile.slice(/^(\S+).peak.table/, 1)
    STDERR.puts file_name
    next unless check_list_hash[ file_name ]
    STDERR.puts "set #{afile} for survey"

    peak_tables << "#{peak_table_dir}/#{afile}"
  end

  STDERR.puts "#{peak_tables.join(":")}"

  msps = CompMSPeaks::MSPeaks.new( 5, min_peak, max_peak, peak_tables )
#  msps.list_sid_by_peak_cluster( "pos" )
#  msps.table_spotcounts_by_peak_cluster( etype )
#  msps.table_peak_cluster( :p )
#  msps.table_peak_cluster( :n )
  msps.table_spotcounts_by_peak_cluster( :p )
#  msps.table_peak_cluster2( :n )
#   sid_list = ["01026", "01038"]
#   msps.list_clade_specific_peaks( :p, sid_list )
#   msps.list_clade_specific_peaks( :n, sid_list )

  STDERR.puts "DONE"

end

