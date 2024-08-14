import ee
import Geo_assets as Ga  # Custom module for geo assets
import geemap  # Library for interactive mapping
import helper as HE  # Custom helper functions
import math  # Standard math library

# Initialize the map
Map = geemap.Map()

def change(pre, post, index): 
    """
    Compute the specified index for both pre and post images and return the quality mosaics 
    and relative difference index.

    Parameters:
    pre (ee.ImageCollection): Pre-event image collection.
    post (ee.ImageCollection): Post-event image collection.
    index (str): Index to be calculated.

    Returns:
    list: List containing quality mosaics for pre and post images and the relative difference index.
    """
    # Calculate the specified index for pre and post images
    index_pre = pre.map(globals()[index])
    index_post = post.map(globals()[index])
    
    # Create quality mosaics for the specified index
    index_pre = index_pre.qualityMosaic(index).rename(index + '_pre')
    index_post = index_post.qualityMosaic(index).rename(index + '_post')

    rd_index = ee.Image().expression(
        '((post - pre) / sqrt((post + pre)))*100',
        {
        'pre': index_pre,
        'post': index_post,
        },).rename('rd' + index)
    
    return [index_pre, index_post, rd_index]

def RVI_V(img):
    """
    Calculate Radar Vegetation Index (RVI) for vertical polarization.

    Parameters:
    img (ee.Image): Sentinel-1 image.

    Returns:
    ee.Image: Image with the RVI_V band.
    """
    RVI_V = img.expression(
    '((4*VH) / (VV + VH))',
    {
        'VH': img.select('VH'),
        'VV': img.select('VV'),
    },
    ).rename('RVI_V')
    return RVI_V

def RFDI(img):
    """
    Calculate Radar Forest Degradation Index (RFDI).

    Parameters:
    img (ee.Image): Sentinel-1 image.

    Returns:
    ee.Image: Image with the RFDI band.
    """
    RFDI = img.normalizedDifference(['VV', 'VH']).rename('RFDI')
    return RFDI

def RVI4S1(img):
    """
    Calculate a specific Radar Vegetation Index for Sentinel-1.

    Parameters:
    img (ee.Image): Sentinel-1 image.

    Returns:
    ee.Image: Image with the RVI4S1 band.
    """
    RVI4S1 = img.expression(
    '(q*(q+3))/((q+1)*(q+1))',
    {
        'q' : img.select('VH').divide(img.select('VV')),
    },
    ).rename('RVI4S1')
    return RVI4S1
