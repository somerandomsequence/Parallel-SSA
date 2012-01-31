#!/usr/bin/env ruby

FloatRE = /\d+\.\d+/

dh = Dir.open(".")
dh.each{ |l|
  next unless l =~ /sa_slave_(\d+)_(#{FloatRE})\.RData/
  nid = $1.to_i
  ts = $2.to_f
  akv_before = nil
  akv_after = nil
  wpe_before = nil
  wpe_after = nil
  points = nil
  ph = File.popen("~/R-2.14.0/bin/R --vanilla --args #{l} < fitness.R","r")
  ph.each_line{ |l2|
    if l2 =~ /FITNESS (#{FloatRE}) (#{FloatRE}) (#{FloatRE}) (#{FloatRE})/
      wpe_before = $1.to_f
      wpe_after = $2.to_f
      akv_before = $3.to_f
      akv_after = $4.to_f
    elsif l2 =~ /SAMPLE (.*)/
      points = $1.strip.split(/\s+/).collect{ |e| e.to_f }
    end
  }
  ph.close
  puts "#{ts} #{nid} #{wpe_before} #{wpe_after} #{akv_before} #{akv_after} #{points.join(",")}"
}
dh.close
   
