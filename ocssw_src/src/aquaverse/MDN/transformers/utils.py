from scipy.spatial.distance import pdist, squareform

import pandas as pd
import numpy as np 


def get_unique_features(features, threshold=0.025, chunksize=10, _bottom=False):
	''' Mask out features which are near duplicates of some other 
		feature(s), based on correlation. The threshold parameter
		determines how similar features can be to keep them.
		Roughly based on https://stackoverflow.com/a/66326102
	'''  
	# To make this more efficient, break into chunks and filter those chunks first
	if (features.shape[1] > chunksize) and not _bottom:
		''' 
		There are two options, which boils down a question of how information 
		propagates between chunks:
			1. Binary-split recursion down to the chunksize. This will filter 
				features more towards the bottom of the tree, assuming it's 
				more likely for similar features to be near each other

			2. Single level tree, with (nearly) all nodes having the same number
				of features (i.e. chunksize). This performs most filtering at the
				root of the tree, when all nodes are combined after the initial filter

			Paradoxically, option 1 appears to result in a higher number of final 
			features returned in spite of applying more filters. This could be a 
			consequence of chained features, e.g. features A and B have a distance
			of 0.015, B and C have a distance of 0.015, and A and C have 0.03; with
			a threshold of 0.025, B and C would both be removed if all three were within 
			the same filter. If instead A and B were in one filter (with B then 
			removed), but C was in another, then the final result would be A and C
			rather than just A. 

			These chains (referred to as "Hamming chains" in the stackoverflow post)
			can be arbitrarily long, and breaking them should intuitively offer 
			a better representation for ML model features. Option 1 therefore seems 
			to be the better choice, as well as using the smallest chunksize possible -
			but note that this hasn't been rigorously tested. 
		'''
		recursion_size = features.shape[-1] // 2  # Option 1 - binary recursion
		# recursion_size = chunksize                # Option 2 - flattened recursion

		mask = np.concatenate([
			get_unique_features(features[:, i: i+recursion_size], threshold, chunksize) 
			for i in range(0, features.shape[1], recursion_size)
		])
		mask[mask] = get_unique_features(features[:, mask], threshold, chunksize, True)
		return mask

	dist = squareform( pdist(features.T, metric='correlation') )
	dist[np.triu_indices_from(dist)] = -1
	mask = (0 <= dist) & (dist <= threshold)
	return ~np.any(mask, axis=1)



def spearmanr(x, y):
	''' Scipy's provided spearmanr is too slow for large matrices, and so we instead use a custom implementation.
		Source: https://stackoverflow.com/questions/52371329/fast-spearman-correlation-between-two-pandas-dataframes/59072032#59072032
	'''
	def rankdata_average(data):
		''' Row-based rankdata using method=mean '''
		dc = np.asarray(data).copy()
		sorter = np.apply_along_axis(np.argsort, 1, data)
		
		inv   = np.empty(data.shape, np.intp)
		ranks = np.tile(np.arange(data.shape[1]), (len(data), 1))
		np.put_along_axis(inv, sorter, ranks, axis=1)

		dc  = np.take_along_axis(dc, sorter, 1)
		res = np.apply_along_axis(lambda r: r[1:] != r[:-1], 1, dc)
		obs = np.column_stack([np.ones(len(res), dtype=bool), res])

		dense = np.take_along_axis(np.apply_along_axis(np.cumsum, 1, obs), inv, 1)
		len_r = obs.shape[1]

		nonzero = np.count_nonzero(obs, axis=1)
		nonzero = pd.Series(nonzero)
		dense = pd.DataFrame(dense)
		obs = pd.DataFrame(obs)

		ranks = []
		for _nonzero, nzdf in obs.groupby(nonzero, sort=False):
			nz = np.apply_along_axis(lambda r: np.nonzero(r)[0], 1, nzdf)

			_count = np.column_stack([nz, np.ones(len(nz)) * len_r])
			_dense = dense.reindex(nzdf.index).values

			_result = 0.5 * (np.take_along_axis(_count, _dense, 1) + np.take_along_axis(_count, _dense - 1, 1) + 1)
			result  = pd.DataFrame(_result, index=nzdf.index)
			ranks.append(result)
		return pd.concat(ranks).sort_index()


	def compute_corr(x, y):
		# Thanks to https://github.com/dengemann
		def ss(a, axis):
			return np.sum(a * a, axis=axis)

		x = np.asarray(x)
		y = np.asarray(y)

		mx = x.mean(axis=-1)
		my = y.mean(axis=-1)

		xm, ym = x - mx[..., None], y - my[..., None]

		r_num = np.add.reduce(xm * ym, axis=-1)
		r_den = np.sqrt(ss(xm, axis=-1) * ss(ym, axis=-1))

		with np.errstate(divide='ignore', invalid="ignore"):
			r = r_num / r_den
		return r


	x = np.asarray(x)
	y = np.asarray(y)

	rx = rankdata_average(x)
	ry = rankdata_average(y)
	return compute_corr(rx, ry)
