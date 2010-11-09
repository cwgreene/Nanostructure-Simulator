def unroll_attr(complex,unrolled):
	result = []
	for attr in complex:
		if attr[0] == "_":
			unroll_attr(complex[attr],unrolled)
		else:
			unrolled[attr]=types[complex[attr]]
	return unrolled

