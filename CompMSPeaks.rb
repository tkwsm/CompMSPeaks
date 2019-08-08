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
# USAGE: ruby comp_peaks.rb ms 100 700 <peak-table-1> <peak-table-2> ...
# Ex: ruby comp_peaks.rb bin-size-of-ms min-ms max-ms <peak-table> <peak-table> 
#
# Format of Peak-Table(s)
# Sample-ID  charge  mol_mass  ret aduct
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

require '../MolCalc/molcalc.rb'

module CompMSPeaks

  def calc_ppm( mol_mass, ppm )
    digit_ppm = 0.000001
    range     = mol_mass.to_f * digit_ppm * ppm.to_f
    return range
  end

  def diff_two_peaks( apeaks, bpeaks)
    ( apeaks - bpeaks )
  end

  def create_raito_bins( original_peaks, ratio )
    original_peaks.sort!
  end

  class MSPeak
    include CompMSPeaks

    def initialize( sid, charge, mol_mass, ret, adduct )
      @sid = sid              # S12002
      @charge   = charge      # pos or neg
      @mol_mass = mol_mass    # 101.00120
      @ret      = ret         # 37.0 etc.   # Retention Time
      @adduct   = adduct      # [H]+ etc.
      @mc       = MolCalc.new
      @ms_without_adduct = calc_ms_without_adduct # 101.00120
    end

    attr_reader :sid, :charge, :mol_mass, :ret, :adducts, :ms_without_adduct 

    def calc_ms_without_adduct
      adduct_mw = @mc.adduct_mw( @adduct )
      ( @mol_mass - adduct_mw.adduct_mw_divided_by_charge ) / adduct_mw.mult
    end

  end

  class MSPeaks
    include CompMSPeaks

    def initialize(bin_ms, min_search_mass, max_search_mass, peak_table_files)
      @min_search_mass = min_search_mass.to_f.round(6)
      @max_search_mass = max_search_mass.to_f.round(6)
      @bin_ms        = bin_ms # 0.001
      @pos_array     = []
      @neg_array     = []
      @c_pos_peaks   = []
      @c_neg_peaks   = []
      @samples       = {}
      @samples_order = []
      @mspeaks = { :p => Hash.new, :n => Hash.new };

      STDERR.puts "The loading files:  #{peak_table_files.join(":")}"
      add_table_files( peak_table_files ) # @pos/@neg_array 
      STDERR.puts "pos stored peaks total #{@pos_array.size}"
      STDERR.puts "neg stored peaks total #{@neg_array.size}"

      create_mspeaks_for_initialize
      create_mspeaks

    end

    attr_reader :bin_ms, :pos_array, :neg_array, :c_pos_peaks, :c_neg_peaks,
                :mspeaks
    def change_samples_order( sample_id_list )
       @samples_order = []
       ff = File.open( sample_id_list )
       ff.rewind
       ff.each do |x|
         next if x !~ /\S/
         next if x =~ /^\#/
         @samples_order << x.slice(/^(\S+).peak.table/, 1)
       end
    end

    def create_mspeaks_for_initialize
      current_bin_max = 0.0
      current_bin_min = @min_search_mass.round(6)
      until current_bin_max > @max_search_mass
        current_bin_max = ( current_bin_min + @bin_ms ).round(6)
        for charge in [:p, :n]
          @mspeaks[charge] = {} if @mspeaks[charge] == nil
          bin_ms_key = [ current_bin_min, current_bin_max ]
          @mspeaks[charge][ bin_ms_key ] = [] if @mspeaks[charge][ bin_ms_key ] == nil
        end
        current_bin_min = current_bin_max
      end
    end

    def create_mspeaks
      @pos_array.each do |item|
        @mspeaks.each_key do |charge| 
          @mspeaks[charge].each_key do |bin_ms_key|
            current_bin_min, current_bin_max = bin_ms_key
            next if charge != :p
            next if item.ms_without_adduct < current_bin_min
            next if item.ms_without_adduct > current_bin_max
            @mspeaks[ charge ][ bin_ms_key ] << item
          end
        end
      end
      @neg_array.each do |item|
        @mspeaks.each_key do |charge| 
          @mspeaks[charge].each_key do |bin_ms_key|
            current_bin_min, current_bin_max = bin_ms_key
            next if charge != :n
            next if item.ms_without_adduct < current_bin_min
            next if item.ms_without_adduct > current_bin_max
            @mspeaks[ charge ][ bin_ms_key ] << item
          end
        end
      end
    end

    def add_original_peaks( sid, charge, mol_mass, ret, adduct )
      @pos_array << MSPeak.new(sid,charge,mol_mass,ret,adduct) if charge == :p
      @neg_array << MSPeak.new(sid,charge,mol_mass,ret,adduct) if charge == :n
    end

    def add_table_data( peak_table_file )
      a = []
      open( peak_table_file ) .each do |x|
        a = x.chomp.split("\s")
        sid = a[0]
        charge   =  :p if a[1] == "pos"
        charge   =  :n if a[1] == "neg"
        mol_mass =  a[2].to_f
        ret      =  a[3].to_f
        adducts   =  a[5..-1]
        @samples[ sid ] = ""
        @samples_order = @samples.keys.sort
        next if mol_mass < @min_search_mass
        next if mol_mass > @max_search_mass
        adducts.each do |an_adduct|
          add_original_peaks( sid, charge, mol_mass, ret, an_adduct )
        end
      end
    end

    def add_table_files( peak_table_files )
      peak_table_files.flatten.each{ |afile| add_table_data( File.expand_path(afile) ) }
    end

    def list_sid_by_peak_cluster( charge )
      charge = :p if charge == "pos"
      charge = :n if charge == "neg"
      @mspeaks[charge].each_key do |mspeak|
        print @mspeaks[charge][mspeak].collect{|x| x.sid }.sort.uniq.join("\t"), "\n" if @mspeaks[charge][mspeak].size > 2
      end
    end

    def table_spotcounts_by_peak_cluster( charge )
      print "mspeak"
      @samples_order.each{|sid| print "\t", sid }
      print "\n"
      charge = :p if charge == "pos"
      charge = :n if charge == "neg"
      @mspeaks[ charge ].each_key do |mspeak|
        print "#{mspeak}"
        tmp_samples = @mspeaks[ charge ][ mspeak ].collect{|item| item.sid } 
        @samples_order.each do |sid|
          print "\t", tmp_samples.count( sid )
        end
         print "\n"
      end
    end

    def table_headers( charge )
      print "mspeak"
      @samples_order.each{|sid| print "\t#{sid}" }
      print "\n"
    end

    def table_peak_cluster( charge )

      @mspeaks[charge].each_key do |mspeak|
        next if @mspeaks[charge][mspeak].size <= 0
        print mspeak.join("-")
        @samples_order.each do |sid|
          if @mspeaks[charge][mspeak].collect{|item| item.sid }.include?( sid )
            countv = @mspeaks[charge][mspeak].count{|item| item.sid == sid }
            print "\t#{countv}"
          else
            print "\t0" 
          end
        end
        print "\n" 
      end

    end

    def list_clade_specific_peaks( charge, sid_list )
      h = {}
      @mspeaks[charge].each_key do |mspeak|
        k = @mspeaks[charge][mspeak].collect{|x| x.sid }.sort.uniq
        next if k != sid_list.sort
        next if @mspeaks[charge][mspeak].size < k.size
        @mspeaks[charge][mspeak].each do |apeak|
          print "#{charge}\t#{mspeak}\t#{apeak.mol_mass}\t#{apeak.ret}\n"
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
    file_name = afile_line.chomp.split("\t")[0]
    next if file_name !~ /\S+.peak.table/
    file_name = file_name.slice(/^(\S+).peak.table/, 1)
    check_list_hash[ file_name ] = ""
  end
  
  STDERR.puts  check_list_hash
  STDERR.puts "READFILES"
  peak_tables    = []
  Dir.open( peak_table_dir ).each do |afile|

    next if afile =~ /^\./
    next if afile !~ /\S+.peak.table/
    file_name = afile.slice(/^(\S+).peak.table/, 1)
    next unless check_list_hash[ file_name ]
    STDERR.puts "set #{afile} for survey"

    peak_tables << "#{peak_table_dir}/#{afile}"
  end

  STDERR.puts "PEAK TABLES Loaded"
  STDERR.puts "#{peak_tables.join(":")}"

##                                   bin, min, max, peak_tables ##
##  msps = CompMSPeaks::MSPeaks.new( 5, min_peak, max_peak, peak_tables )
  msps = CompMSPeaks::MSPeaks.new( 0.001, min_peak, max_peak, peak_tables )
  msps.change_samples_order( sample_id_list )
#  msps.list_sid_by_peak_cluster( "pos" )
#  msps.table_spotcounts_by_peak_cluster( "pos" )
  msps.table_headers( :p )
  msps.table_peak_cluster( :p )
  msps.table_peak_cluster( :n )
#  msps.table_spotcounts_by_peak_cluster( :p )
#  msps.table_peak_cluster2( :n )
#   sid_list = ["10086", "10199"]
#   msps.list_clade_specific_peaks( :p, sid_list )
#   msps.list_clade_specific_peaks( :n, sid_list )

  STDERR.puts "DONE"

end

