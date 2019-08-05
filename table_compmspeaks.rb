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


# require "/home/takeshik/scripts/Tkwsm/CompMSPeaks/CompMSPeaks.rb"
require "./CompMSPeaks.rb"

if ARGV[0] !~ /\d/ or ARGV[1] !~ /\d/ or ARGV[0] =~ /\D/
  STDERR.print "Please check the USAGE\n or \n"
  STDERR.print "ruby CompMSPeaks.rb 100, 2000, <peaktable1> ...\n\n"
  exit
end

STDERR.puts "START"

min_peak       = ARGV.shift
max_peak       = ARGV.shift
sample_id_list = ARGV.shift   # 11011 9913 meat  natural  Eumetazoa  Bos taurus
elec_type      = ARGV.shift   # p or n # p: positive, n: negative
peak_table_dir = ARGV.shift

if    elec_type == "p"
  etype = :p 
elsif elec_type == "n"
  etype = :n 
else
  STDERR.puts "something error"
  exit
end

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

msps.table_spotcounts_by_peak_cluster( etype )

STDERR.puts "DONE"


