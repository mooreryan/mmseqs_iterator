#!/usr/bin/env ruby
Signal.trap("PIPE", "EXIT")

require "set"

iter = ARGV[0]
btab = ARGV[1]

connections = Set.new
File.open(btab, "rt").each_line do |line|
  source, target, *rest = line.chomp.split "\t"

  connections << [source, target]
end

connections.each do |(source, target)|
  puts [source, target, iter].join "\t"
end
