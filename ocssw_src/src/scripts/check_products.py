#!/usr/bin/env python3

import xml.etree.ElementTree as et
import argparse

def float_or_double(prod_type: str) -> bool:
	return prod_type == 'float' or prod_type == 'double'

def set_integer_bounds(integer_type: str):
	if integer_type == 'short':
		return -32768, 32767
	elif integer_type == 'int':
		return  -2147483648, 2147483647
	elif integer_type == 'byte':
		return -128, 127
	else:
		return 0, 0

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

product_string     = '{http://oceancolor.gsfc.nasa.gov}product'
type_string        = '{http://oceancolor.gsfc.nasa.gov}type'
range_string       = '{http://oceancolor.gsfc.nasa.gov}range'
validMax_string    = '{http://oceancolor.gsfc.nasa.gov}validMax'
validMin_string    = '{http://oceancolor.gsfc.nasa.gov}validMin'
displayMin_string  = '{http://oceancolor.gsfc.nasa.gov}displayMin'
displayMax_string  = '{http://oceancolor.gsfc.nasa.gov}displayMax'
addOffset_string   = '{http://oceancolor.gsfc.nasa.gov}addOffset'
scaleFactor_string = '{http://oceancolor.gsfc.nasa.gov}scaleFactor'
algorithm_string   = '{http://oceancolor.gsfc.nasa.gov}algorithm'

xml_tree      = et.parse(args.filename)
xml_tree_root = xml_tree.getroot()
for product in xml_tree_root.findall(product_string):
	product_name 			 = product.attrib['name']
	default_range            = None
	product_has_error 		 = False	
	product_invalid_validMin = False
	product_invalid_validMax = False
	product_missing_type     = False
	product_missing_validMax_or_validMin    = False	
	product_missing_default_range           = False
	product_validMax, product_validMin 	= 0xdeadbeef, 0xdeadbeef
	product_addOffset, product_scaleFactor 	= None, None

	if args.verbose:
		print(product_name)
	else:
		(error_string := f'{product_name}')
	
	try:
		product_type = product.find(type_string).text
	except: #TODO: Default to float
		product_type = 'float'
		# product_has_error = True
		# product_missing_type = True
		# (product_type := 'Missing')
		# if args.verbose:
		# 	print("    Missing type")
		# else:
		# 	error_string += "\n    Missing type"

	try:
		(default_range := product.find(range_string))
		_ = default_range.text # Try to induce an error for missing default range
	except:
		product_has_error = True
		product_missing_default_range = True
		if args.verbose:
			print("    Missing default range")
		else:
			error_string += "\n    Missing default range"

	try:
		product_validMax = default_range.find(validMax_string).text
	except:
		product_has_error = True
		product_missing_validMax_or_validMin = True
		if product_missing_type is False and product_missing_default_range is False:
			if args.verbose:
				print("    Missing validMax")
			else:
				error_string += "\n    Missing validMax"
			
	try:
		product_validMin = default_range.find(validMin_string).text
	except:
		product_has_error 	     = True
		product_missing_validMax_or_validMin = True
		if product_missing_type is False and product_missing_default_range is False:
			if args.verbose:
				print("    Missing validMin")
			else:
				error_string += "\n    Missing validMin" 
	

	# Try to assign scaleFactor/addOffset. If any errors, they fall back to 1 and 0, respectively
	try:
		if not float_or_double(product_type):
			product_addOffset = float(default_range.find(addOffset_string).text)
			product_scaleFactor = float(default_range.find(scaleFactor_string).text)
		else:
			product_addOffset = 0
			product_scaleFactor = 1	
			
		if ((product_addOffset is not None and product_scaleFactor is not None) 
			and (product_addOffset != 0 and product_scaleFactor != 1)
			and float_or_double(product_type)):
			# Getting here means that a float/double has a scaleFactor/addOffset
			raise TypeError # Technically not the right kind of exception; used just to differentiate
	except AttributeError as ae:
		if str(ae) == "'NoneType' object has no attribute 'text'":
			pass
		else:
			product_has_error = True
			if product_missing_default_range is False and not float_or_double(product_type):
				if args.verbose:
					print("    If scaleFactor/addOffset is defined, both must be")
				else:
					error_string += "\n    If scaleFactor/addOffset is defined, both must be"
		product_addOffset = 0
		product_scaleFactor = 1
	except TypeError as te:
		product_has_error = True
		if (product_missing_default_range is False and float_or_double(product_type) and
				product_addOffset is not None and product_scaleFactor is not None):
			if args.verbose:
				print(f"    {product_type} with a scaleFactor/addOffset")
			else:
				error_string += f"\n    {product_type} with scaleFactor/addOffset"
		product_addOffset = 0
		product_scaleFactor = 1

	try:
		_ = default_range.find(displayMin_string).text # Attempt to induce an exception 
		_ = default_range.find(displayMax_string).text # Attempt to induce an exception
	except: 
		product_has_error = True
		if product_missing_type is False and product_missing_default_range is False:
			if args.verbose:
				print("    Missing displayMin/displayMax")
			else:
				error_string += "\n    Missing displayMin/displayMax"				

	###### validMin/validMax checking ######
	int_min, int_max = set_integer_bounds(product_type)
	min_possible = (int_min * product_scaleFactor) + product_addOffset
	max_possible = (int_max * product_scaleFactor) + product_addOffset
	# min_is_valid = min_possible <= float(product_validMin)
	try:	
		min_is_valid = float(product_validMin) >= min_possible
	except:
		product_has_error = True
		if args.verbose:
			print(f"    Issue parsing validMin, {product_validMin} cannot be parsed as a number")
		else:
			error_string += f"\n    Issue parsing validMin of {product_validMin}"
	try:	
		max_is_valid = float(product_validMax) <= max_possible
	except ValueError:
		product_has_error = True
		if args.verbose:
			print(f"    Issue parsing validMax, {product_validMax} cannot be parsed as a number")
		else:
			error_string += f"\n    Issue parsing validMax of {product_validMax}"

	if float_or_double(product_type): # Assume scaleFactor = 1, addOffset = 0
		min_is_valid = True
		max_is_valid = True
	else:
		try:
			if (float(product_validMin) > float(product_validMax) and 
				product_missing_validMax_or_validMin is False):
				product_has_error = True
				if args.verbose:
					print(f"    validMin of {product_validMin} is greater than validMax of {product_validMax}")
				else:
					error_string += f"\n    validMin of {product_validMin} is greater than validMax of {product_validMax}"
		except:
			pass

	if not min_is_valid and product_missing_validMax_or_validMin is False and product_missing_type is False:
		product_has_error 		 = True
		product_invalid_validMin = True
		if args.verbose:
			print(f"    Invalid validMin:", end= ' ') # to include ↓ would make too long a string in the editor
			print(f"Minimum possible value is {min_possible} but validMin is {product_validMin}")
		else:
			error_string += f"\n    Invalid validMin, Min possible value is {min_possible}"

	if not max_is_valid and product_missing_validMax_or_validMin is False and product_missing_type is False:
		product_has_error 		 = True
		product_invalid_validMax = True
		if args.verbose:
			print(f"    Invalid validMax:", end= ' ') # to include ↓ would make too long a string in the editor
			print(f"Max possible value is {max_possible} but validMax is {product_validMax}")
		else:
			error_string += f"\n    Invalid validMax, Max possible value is {max_possible}"


	### Algo checking ###
	for algo in product.findall(algorithm_string):
		algo_validMin 	  = product_validMin
		algo_validMax 	  = product_validMax
		algo_scaleFactor  = product_scaleFactor
		algo_addOffset 	  = product_addOffset
		algo_name         = product_name
		algo_error_string = ''		
		algo_has_error    = False
	
		try:
			algo_name = algo.attrib['name']	
		except: # Default algo
			algo_name = product_name
		algo_error_string += f"\n     {algo_name}"

		if args.verbose:
			print(f"       {algo_name}")

		try:
			(algo_range := algo.find(range_string))
		except: # Algo ranges default to product range
			(algo_range := default_range)

		try:
			algo_validMax    = algo_range.find(validMax_string).text
			algo_validMin    = algo_range.find(validMin_string).text
			algo_addOffset   = float(algo_range.find(addOffset_string).text)
			algo_scaleFactor = float(algo_range.find(scaleFactor_string).text)
		except: # Algos default to product's values.
			algo_validMax 	 = product_validMax
			algo_validMin 	 = product_validMin
			algo_addOffset 	 = product_addOffset
			algo_scaleFactor = product_scaleFactor

		### Algo validity checking ###
		algo_min_possible     = (int_min * product_scaleFactor) + product_addOffset
		algo_max_possible     = (int_max * product_scaleFactor) + product_addOffset
		try:
			algo_invalid_validMin = (min_possible > float(algo_validMin))
		except:
			if product_has_error: #Assume it would be the same error
				pass
			else:
				if args.verbose:
					print(f"    Issue parsing validMin, {algo_validMin} cannot be parsed as a number")
				else:
					error_string += f"\n    Issue parsing validMin {algo_validMin}"
		
		try:
			algo_invalid_validMax = (float(algo_validMax) > max_possible)
		except:
			if product_has_error: #Assume it would be the same error
				pass
			else:
				if args.verbose:
					print(f"    Issue parsing validMax, {algo_validMin} cannot be parsed as a number")
				else:
					error_string += f"\n    Issue parsing validMax: {algo_validMax}"

		if float_or_double(product_type): # min and max should be all valid b/c sF of 1, aO of 0?
			algo_invalid_validMin = False
			algo_invalid_validMax = False

		if algo_invalid_validMin and not product_invalid_validMin:
			if product_missing_type is False and product_missing_validMax_or_validMin is False:
				product_has_error = True
				algo_has_error 	  = True
				if args.verbose: 
					print(f"    Invalid validMin:", end= ' ') # to include ↓ would make too long a string in the editor
					print(f"Minimum possible value is {min_possible} but validMin is {algo_validMin}")
				else:
					algo_error_string += f"\n        Invalid validMin, Min possible value is {min_possible}"
			else:
				pass
				
		if algo_invalid_validMax and not product_invalid_validMax: 
			if product_missing_type is False and product_missing_validMax_or_validMin is False:
				product_has_error = True
				algo_has_error 	  = True
				if args.verbose: 
					print(f"    Invalid validMax:", end= ' ') # to include ↓ would make too long a string in the editor
					print(f"Maximum possible value is {max_possible} but validMax is {algo_validMax}")
				else:
					algo_error_string += f"\n        Invalid validMax, Max possible value is {max_possible}"
			else:
				pass

		if algo_has_error:
			error_string += algo_error_string
	product_type = 'none'

	# If args.verbose, we've been printing this whole time, no need to print the error string
	if not args.verbose and product_has_error:
		print(error_string)	
