#! /usr/bin/env ruby

def convert_eda(in_file, out_file)

	def hex_string_to_ascii(str)
		str = "0#{str}" if str.size.odd?
		h = str.scan(/../).map { |pair| pair.hex.chr }.join
	end
	
	def parse_signed(bytes)
		# Get the sign
		if bytes[0] == "1"
			sgn = "-"
		else
			sgn = ""
		end
	
		# Get the value	
		value = sprintf '%03d', bytes[3..16].to_i(2)
	
		# Get the floating point
		case bytes[1..2]
			when "00"
				offset = "0.xxx"
				value_off = sgn + "0." + value.to_s
			when "01"
				offset = "x.xx"
				value_off = sgn + value.to_s[0] + "." + value.to_s[1..2]
			when "10"
				offset = "xx.x"
				value_off = sgn + value.to_s[0..1] + "." + value.to_s[2]			
			when "11"
				offset = "xxx."
				value_off = sgn + value.to_s[0..2] + "."			
		end	
		return value_off
	end
	
	def parse_unsigned(bytes)
		# Get the value
		value = sprintf '%04d', bytes[2..16].to_i(2)
		
		# Get the floating point
		case bytes[0..1]
			when "00"
				offset = "0.xxxx"
				value_off = "0." + value.to_s			
			when "01"
				value_off = value.to_s[0] + "." + value.to_s[1..3]
			when "10"
				offset = "xx.xx"
				value_off = value.to_s[0..1] + "." + value.to_s[2..3]			
			when "11"
				offset = "xxx.x"
				value_off = value.to_s[0..2] + "." + value.to_s[3]						
		end	
		return value_off
	end

	f = File.open(in_file, "rb").read
	puts "start " + in_file
	# Convert the file to hex for easier processing
	g = f.unpack('H*')
	
	# Split the file into a header chunk (array) and a data chunk (long string)
	str = g[0].split("0d0a")
	
	## get the running variables that we'll use throughout ##
	file_data = in_file.split("/")[2]
	puts file_data
	file_number = file_data[3..4]
	bracelet_id = file_data[6..9]
	bracelet_start_time = hex_string_to_ascii(str[5])[12..30]

	# Parse the data
	samps = str[8..str.length].join("0d0a")
	pos = 0
	sample_number = 1
	begin
		sample = samps[pos..(pos+27)]
		if sample == "e3e7e3e7e3e7e70fe3e7e70fe5e2"
			z_val = ""
			y_val = ""
			x_val = ""
			b_val = ""
			t_val = ""
			e_val = ""
		else
			# Z acceleration
			z = sample[0..3]
			z_b = "0"*(16-z.to_i(16).to_s(2).length).to_i + z.to_i(16).to_s(2)
			z_val = parse_signed(z_b)
				
			# Y acceleration
			y = sample[4..7]
			y_b = "0"*(16-y.to_i(16).to_s(2).length).to_i + y.to_i(16).to_s(2)	
			y_val = parse_signed(y_b)	
			
			# X acceleration
			x = sample[8..11]	
			x_b = "0"*(16-x.to_i(16).to_s(2).length).to_i + x.to_i(16).to_s(2)	
			x_val = parse_signed(x_b)
			
			# Battery voltage
			b = sample[12..15]	
			b_b = "0"*(16-b.to_i(16).to_s(2).length).to_i + b.to_i(16).to_s(2)	
			b_val = parse_unsigned(b_b)	
			
			# sensor temperature
			t = sample[16..19]	
			t_b = "0"*(16-t.to_i(16).to_s(2).length).to_i + t.to_i(16).to_s(2)	
			t_val = parse_signed(t_b)	
			
			# EDA
			e = sample[20..23]	
			e_b = "0"*(16-e.to_i(16).to_s(2).length).to_i + e.to_i(16).to_s(2)	
			e_val = parse_unsigned(e_b)	
		end				
		val_arr = Array[bracelet_id, file_number, bracelet_start_time, sample_number, z_val, y_val, x_val, b_val, t_val, e_val]
		write_line = val_arr.join(",") + "\n"		
		File.open(out_file, 'a') { |f| f.write(write_line) }				
		pos += 28
		sample_number += 1
	end while pos < samps.length
	puts "end" + in_file	
end
