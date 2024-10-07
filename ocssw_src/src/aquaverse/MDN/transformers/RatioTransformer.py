from ._CustomTransformer import _CustomTransformer
from .utils import spearmanr, get_unique_features

from itertools import combinations
from functools import wraps, partial
from pathlib import Path 

import pandas as pd
import numpy as np 
import math 

import pandas as pd
import numpy as np


def add_labels(n_bands, label_function):
	''' Decorator which adds labels to the RatioTransformer.labels
		list, using the given label_function to combine lists of
		wavelengths that are passed to the decorated function. '''
	def decorator(function):
		setattr(function, 'n_bands', n_bands)

		@wraps(function)
		def wrapper(Rrs, labels, *args, **kwargs):
			# Cast to a series to handle single values, and ensure they are strings
			as_str_series = lambda v: pd.Series(v, dtype=str) 
			labels.extend( map(label_function, *map(as_str_series, args)) )
			return function(Rrs, *args, **kwargs)
		return wrapper
	return decorator


#=========================================
# Below functions represent commonly used
# product formulations found in literature
#=========================================  

@add_labels(2, lambda W1, W2: f'{W2}/{W1}')
def BR2(Rrs, W1, W2):
	''' 2-Band ratio '''
	return Rrs(W2) / Rrs(W1)


@add_labels(3, lambda W1, W2, W3: f'{W3}/{W1}-{W3}/{W2}')
def BR3(Rrs, W1, W2, W3):
	''' 3-Band ratio '''
	return Rrs(W3)/Rrs(W1) - Rrs(W3)/Rrs(W2)


@add_labels(3, lambda W1, W2, W3: f'{W1}|{W2}|{W3}')
def LH(Rrs, W1, W2, W3):
	''' Line height '''
	c = (W3 - W2) / (W3 - W1)
	return Rrs(W2) - c*Rrs(W1) - (1-c)*Rrs(W3)


@add_labels(2, lambda W1, W2: f'{W2}-{W1}/{W2}+{W1}')
def ND(Rrs, W1, W2):
	''' Normalized difference '''
	return (Rrs(W2)-Rrs(W1)) / (Rrs(W2)+Rrs(W1))


@add_labels(3, lambda W1, W2, W3: f'({W3}+{W1})/2-{W2}')		
def AVG(Rrs, W1, W2, W3):
	''' Average difference '''
	return (Rrs(W3)+Rrs(W1))/2 - Rrs(W2)


@add_labels(2, lambda W1, W2: f'peak|{W1}-{W2}') 
def PEAK(Rrs, min_wvl, max_wvl, wavelengths=None):
	''' Peak location in a range of wavelengths '''
	wave = np.array(wavelengths)
	area = (wave[:, None] >= min_wvl) & (wave[:, None] <= max_wvl)
	idxs = [Rrs( wave[a] ).argmax(axis=-1) for a in area.T]
	peak = [wave[a][i] - wave[a].min() for a,i in zip(area.T, idxs)]
	return np.stack(peak, axis=1)

#=========================================  



class RatioTransformer(_CustomTransformer):	
	''' Add ratio features '''
	functions = [BR2, BR3, LH, ND, AVG]

	def __init__(self, wavelengths, sensor, *args, excl_Rrs=False, all_ratio=False, cutoff=0.5, **kwargs):
		self.wavelengths = list(wavelengths)
		self.all_ratio   = all_ratio
		self.excl_Rrs    = excl_Rrs
		self.sensor      = sensor
		self.cutoff      = cutoff


	@staticmethod
	def config_info(*args, **kwargs):
		transformer = RatioTransformer(*args, **kwargs)
		return transformer.get_n_features()


	def get_n_features(self):
		n_wavelengths  = len(self.wavelengths)
		n_combinations = lambda n, r, f=math.factorial: f(n) // f(r) // f(n-r) 
		function_names = ', '.join([f.__name__ for f in self.functions])
		function_bands = [f.n_bands for f in self.functions]
		total_number   = sum(map(partial(n_combinations, n_wavelengths), function_bands))
		return f'[{function_names}] applied to {n_wavelengths} wavelengths yields up to {total_number} features'


	def _inverse_transform(self, X, *args, **kwargs): 
		if self.excl_Rrs: raise Exception('No inverse when Rrs is excluded from features')
		return np.array(X)[..., :self.shape]


	def _fit(self, X, y, *args, **kwargs):
		from ..benchmarks.utils import has_band, closest_wavelength

		self.features = []
		self.labels   = []
		self.shape    = X.shape[-1]

		Rrs_df = pd.DataFrame(X, columns=self.wavelengths)
		Rrs    = lambda cols: Rrs_df[cols].values

		# Steps:
		# 	1. Enumerate all wavelength combinations for the selected feature functions
		# 	2. Calculate correlations between features and target variables, filtering those below cutoff
		# 	3. Filter remaining combinations by removing those which are too similar to others
		# 	4. Store the final remaining wavelength combinations for use in the _transform function
		for function in self.functions:
			wavelengths = list(combinations(self.wavelengths, function.n_bands))
			values      = function(Rrs, [], *map(np.array, zip(*wavelengths)))
			correlation = [spearmanr(yy[np.isfinite(yy), None].T, values[np.isfinite(yy)].T) for yy in y.T]
			mask        = np.max(correlation, axis=0) > self.cutoff
			mask[mask] &= get_unique_features(values[:, mask], threshold=0.05)

			mask_index  = pd.MultiIndex.from_tuples(wavelengths)
			mask_frame  = pd.DataFrame(mask, index=mask_index)
			wavelengths = mask_frame[mask_frame].dropna().index.to_frame()
			self.features.append( (function, wavelengths.values.T) )

		# Add any additional manually specified features to the list
		additional = [
			(partial(PEAK, wavelengths=self.wavelengths), [
				(450, 520), # Peak in blue
				(550, 600), # Peak in green
				(620, 670), # Peak in red
				(680, 720), # Peak in nir
			])
		]

		# Ensure required wavelengths are available, and use the exact values available
		for function, wavelengths in additional:
			valid       = lambda wvls: all([has_band(w, self.wavelengths) for w in wvls])
			select_wvls = closest = np.array(wavelengths)[list(map(valid, wavelengths))]

			if len(select_wvls): 
				closest = closest_wavelength(select_wvls.flatten(), self.wavelengths)
			wavelengths = closest.reshape(select_wvls.shape).T
			self.features.append( (function, wavelengths) )

		# Run the transformation to add the labels being used
		self._transform(X, labels=self.labels)
		print(f'\nFit RatioTransformer with {len(self.labels)} features')


	def _transform(self, X, *args, labels=None, **kwargs):		 
		Rrs_df = pd.DataFrame(X, columns=self.wavelengths)
		Rrs    = lambda cols: Rrs_df[cols].values 
		x_new  = [] if self.excl_Rrs else [Rrs_df.values]

		if labels is None:
			labels = []

		if not self.excl_Rrs:
			labels.extend( list(map(str, self.wavelengths)) )

		for function, wavelengths in self.features:
			x_new.append( function(Rrs, labels, *wavelengths) )
		x_new = np.concatenate(x_new, axis=-1)
		x_new = np.maximum(-1e8, np.minimum(1e8, x_new))

		# Sanity checks
		assert(x_new.shape[-1] == len(labels)), f'Mismatch between features and labels: {x_new.shape[-1]} vs {len(labels)}'
		assert(len(self.labels)== len(labels)), f'Mismatch between old and new labels: {len(self.labels)} vs {len(labels)}'
		return x_new
