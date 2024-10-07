'''
Mixture Density Network (MDN) ocean color water quality product retrieval 

MDN is a machine learning algorithm, trained to use remote sensing reflectance (Rrs) 
to estimate various water quality products. This package includes the model which
retrieves Chlorophyll-a (chl), Total Suspended Solids (tss), Colored Dissolved Organic 
Matter at 440nm (cdom), and Phycocyanin Concentration (pc) from HICO imagery. 

To utilize this package, activate the provided virtual environment and call the script:
	$ source venv/Scripts/activate
	$ python create_images.py

Code base can additionally be found at:
	https://github.com/BrandonSmithJ/MDN
	https://github.com/STREAM-RS/STREAM-RS

Brandon Smith, NASA Goddard Space Flight Center, October 2021
'''

from netCDF4 import Dataset
from pathlib import Path 

from MDN.parameters import get_args
from MDN.utils import get_sensor_bands, closest_wavelength
from MDN import image_estimates

import matplotlib.colors as colors
import matplotlib.pyplot as plt 
import numpy as np
import time 


def gamma_stretch(data, gamma=2): 
	''' Apply gamma stretching to brighten imagery '''
	return (255. * data ** 0.5).astype(np.uint8) 



def extract_data(image, avail_bands, req_bands, allow_neg=False, key='Rrs'):
	''' Extract the requested bands from a given NetCDF object '''

	def extract(requested):
		bands = [closest_wavelength(band, avail_bands) for band in requested]
		return np.ma.stack([image[f'{key}_{band}'][:] for band in bands], axis=-1)
	
	# Extract the requested bands from the image object	 
	extracted = extract(req_bands)

	# Set any values <= 0 to nan if we disallow negatives
	if not allow_neg: 
		extracted[extracted <= 0] = np.nan
	
	# Return the data, filling any masked values with nan
	return extracted.filled(fill_value=np.nan)



def plot_product(ax, title, product, rgb, vmin, vmax):
	''' Plot a given product on the axis using vmin/vmax as the 
		colorbar min/max, and rgb as the visible background '''
	ax.imshow( gamma_stretch(rgb) )
	ax.axis('off')
	ax.set_title(key.upper())

	norm = colors.LogNorm(vmin=vmin, vmax=vmax)
	img  = ax.imshow(np.squeeze(product), norm=norm, cmap='turbo')
	plt.colorbar(img, ax=ax)




if __name__ == '__main__':
	sensor = 'HICO'
	kwargs = {
		'sensor'        : sensor,
		'product'       : 'chl,tss,cdom,pc',
		'sat_bands'     : True,
		'use_ratio'     : True,
		'use_excl_Rrs'  : True,
	}

	# Load the bands required for the given sensor
	req_bands = get_sensor_bands(sensor, get_args(**kwargs))
	rgb_bands = [660, 550, 440]

	for location in Path(f'{sensor}-imagery').glob('*'):
		time_start = time.time()

		# Load HICO data, using rhos as the visible background
		image = Dataset(location.joinpath('l2gen.nc'))['geophysical_data']
		bands = sorted([int(k.replace('Rrs_', '')) for k in image.variables.keys() if 'Rrs_' in k])
		Rrs   = extract_data(image, bands, req_bands)
		rgb   = extract_data(image, bands, rgb_bands, key='rhos')

		# Generate product estimates - 'slices' contains the index of each product within 'products'
		products, slices = image_estimates(Rrs, **kwargs)
		
		# Create plot for each product, bounding the colorbar per product
		f, axes = plt.subplots(1, len(slices), figsize=(4*len(slices), 8))
		bounds  = {
			'chl' : (1,  100),
			'tss' : (1,  100),
			'pc'  : (1,  100),
			'cdom': (0.1, 10),
		}
		for i, (key, idx) in enumerate(slices.items()):
			plot_product(np.atleast_1d(axes)[i], key, products[..., idx], rgb, *bounds[key])
		plt.tight_layout()
		plt.savefig(f'{sensor}_{location.stem}.png')
		plt.clf()

		print(f'Generated {sensor}_{location.stem}.png in {time.time()-time_start:.1f} seconds')

