def to_matlab_array(s)
	 s.scan( /[[:digit:]],[[:digit:]]/ ).join(';')
end
# dann als matlab array und zb scatter(a[:,1],a[:,2])
