import geemap  # Library for interactive mapping with Google Earth Engine
import os  # Standard library for interacting with the operating system
import ee  # Google Earth Engine library
import SAR_indices  # Custom module for SAR indices
import Geo_assets as Ga  # Custom module for geographic assets

import border_noise_correction as bnc  # Custom module for border noise correction
import speckle_filter as sf  # Custom module for speckle filtering
import terrain_flattening as tf  # Custom module for terrain flattening
import helper as help  # Custom helper functions

###########################################
# DO THE JOB
###########################################

def s1_preproc(params):
    """
    Preprocess Sentinel-1 data for the specified indices and parameters.

    Parameters:
    params (dict): Dictionary containing preprocessing parameters.

    Returns:
    geemap.Map: Map object with preprocessed data layers.
    """
    # Extract parameters from the dictionary
    APPLY_BORDER_NOISE_CORRECTION = params['APPLY_BORDER_NOISE_CORRECTION']
    APPLY_TERRAIN_FLATTENING = params['APPLY_TERRAIN_FLATTENING']
    APPLY_SPECKLE_FILTERING = params['APPLY_SPECKLE_FILTERING']
    POLARIZATION = params['POLARIZATION']
    PLATFORM_NUMBER = params['PLATFORM_NUMBER']
    ORBIT = params['ORBIT']
    ORBIT_NUM = params['ORBIT_NUM']
    SPECKLE_FILTER_FRAMEWORK = params['SPECKLE_FILTER_FRAMEWORK']
    SPECKLE_FILTER = params['SPECKLE_FILTER']
    SPECKLE_FILTER_KERNEL_SIZE = params['SPECKLE_FILTER_KERNEL_SIZE']
    SPECKLE_FILTER_NR_OF_IMAGES = params['SPECKLE_FILTER_NR_OF_IMAGES']
    TERRAIN_FLATTENING_MODEL = params['TERRAIN_FLATTENING_MODEL']
    DEM = params['DEM']
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = params['TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER']
    FORMAT = params['FORMAT']
    START_DATE = params['START_DATE']
    EVENT_DATE = params['EVENT_DATE']
    STOP_DATE = params['STOP_DATE']
    INDEX = params['INDEX']
    ROI = params['ROI']
    CLIP_TO_ROI = params['CLIP_TO_ROI']
    MAKE_MAP = params['MAKE_MAP']
    SAVE_ASSET = params['SAVE_ASSET']
    Filename = params['FILENAME']

    ###########################################
    # 0. CHECK PARAMETERS
    ###########################################

    # Set default values for parameters if not specified
    if APPLY_BORDER_NOISE_CORRECTION is None:
        APPLY_BORDER_NOISE_CORRECTION = True
    if APPLY_TERRAIN_FLATTENING is None:
        APPLY_TERRAIN_FLATTENING = True
    if APPLY_SPECKLE_FILTERING is None:
        APPLY_SPECKLE_FILTERING = True
    if POLARIZATION is None:
        POLARIZATION = 'VVVH'
    if ORBIT is None:
        ORBIT = 'BOTH'
    if SPECKLE_FILTER_FRAMEWORK is None:
        SPECKLE_FILTER_FRAMEWORK = 'MULTI'
    if SPECKLE_FILTER is None:
        SPECKLE_FILTER = 'GAMMA MAP'
    if SPECKLE_FILTER_KERNEL_SIZE is None:
        SPECKLE_FILTER_KERNEL_SIZE = 7
    if SPECKLE_FILTER_NR_OF_IMAGES is None:
        SPECKLE_FILTER_NR_OF_IMAGES = 10
    if TERRAIN_FLATTENING_MODEL is None:
        TERRAIN_FLATTENING_MODEL = 'VOLUME'
    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER is None:
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0
    if FORMAT is None:
        FORMAT = 'DB'

    # Validate polarization parameter
    pol_required = ['VV', 'VH', 'VVVH']
    if POLARIZATION not in pol_required:
        raise ValueError("ERROR!!! Parameter POLARIZATION not correctly defined")

    # Validate orbit parameter
    orbit_required = ['ASCENDING', 'DESCENDING', 'BOTH']
    if ORBIT not in orbit_required:
        raise ValueError("ERROR!!! Parameter ORBIT not correctly defined")

    # Validate terrain flattening model parameter
    model_required = ['DIRECT', 'VOLUME']
    if TERRAIN_FLATTENING_MODEL not in model_required:
        raise ValueError("ERROR!!! Parameter TERRAIN_FLATTENING_MODEL not correctly defined")

    # Validate format parameter
    format_required = ['LINEAR', 'DB']
    if FORMAT not in format_required:
        raise ValueError("ERROR!!! FORMAT not correctly defined")

    # Validate speckle filter framework parameter
    frame_needed = ['MONO', 'MULTI']
    if SPECKLE_FILTER_FRAMEWORK not in frame_needed:
        raise ValueError("ERROR!!! SPECKLE_FILTER_FRAMEWORK not correctly defined")

    # Validate speckle filter parameter
    format_sfilter = ['BOXCAR', 'LEE', 'GAMMA MAP', 'REFINED LEE', 'LEE SIGMA']
    if SPECKLE_FILTER not in format_sfilter:
        raise ValueError("ERROR!!! SPECKLE_FILTER not correctly defined")

    # Validate additional layover/shadow buffer parameter
    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0:
        raise ValueError("ERROR!!! TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER not correctly defined")

    # Validate speckle filter kernel size parameter
    if SPECKLE_FILTER_KERNEL_SIZE <= 0:
        raise ValueError("ERROR!!! SPECKLE_FILTER_KERNEL_SIZE not correctly defined")
    
    # Validate index parameter
    index_required = ['RVI_V', 'RFDI', 'RVI4S1']
    if INDEX not in index_required:
        raise ValueError("ERROR!!! Parameter INDEX not correctly defined")

    ###########################################
    # 1. DATA SELECTION
    ###########################################

    # Select pre-event and post-event Sentinel-1 image collections
    pre = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.eq('resolution_meters', 10)) \
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
        .filterDate(START_DATE, EVENT_DATE) \
        .filterBounds(ROI)\
        .select(['VV', 'VH', 'angle'])
    
    post = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')\
        .filter(ee.Filter.eq('instrumentMode', 'IW'))\
        .filter(ee.Filter.eq('resolution_meters', 10)) \
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
        .filterDate(EVENT_DATE, STOP_DATE) \
        .filterBounds(ROI)\
        .select(['VV', 'VH', 'angle'])

    # Filter by platform number if specified
    if PLATFORM_NUMBER in ['A', 'B']:
        pre = pre.filter(ee.Filter.eq('platform_number', PLATFORM_NUMBER))
        post = post.filter(ee.Filter.eq('platform_number', PLATFORM_NUMBER))
   
    # Filter by orbit number if specified
    if ORBIT_NUM is not None:
        pre = pre.filter(ee.Filter.eq('relativeOrbitNumber_start', ORBIT_NUM))
        post = post.filter(ee.Filter.eq('relativeOrbitNumber_start', ORBIT_NUM))

    # Select orbit type
    if ORBIT != 'BOTH':
        pre = pre.filter(ee.Filter.eq('orbitProperties_pass', ORBIT))
        post = post.filter(ee.Filter.eq('orbitProperties_pass', ORBIT))

    # Print the number of images in each collection
    print('Number of images in collection: ', pre.size().getInfo())
    print('Number of images in collection: ', post.size().getInfo())

    ###########################################
    # 2. ADDITIONAL BORDER NOISE CORRECTION
    ###########################################

    if APPLY_BORDER_NOISE_CORRECTION:
        pre = pre.map(bnc.f_mask_edges)
        post = post.map(bnc.f_mask_edges)
        print('Additional border noise correction is completed')

    ###########################################
    # 3. SPECKLE FILTERING
    ###########################################

    if APPLY_SPECKLE_FILTERING:
        if SPECKLE_FILTER_FRAMEWORK == 'MONO':
            pre = ee.ImageCollection(sf.MonoTemporal_Filter(pre, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER))
            post = ee.ImageCollection(sf.MonoTemporal_Filter(post, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER))
            print('Mono-temporal speckle filtering is completed')
        else:
            pre = ee.ImageCollection(sf.MultiTemporal_Filter(pre, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER, SPECKLE_FILTER_NR_OF_IMAGES))
            post = ee.ImageCollection(sf.MultiTemporal_Filter(post, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER, SPECKLE_FILTER_NR_OF_IMAGES))
            print('Multi-temporal speckle filtering is completed')

    ###########################################
    # 4. TERRAIN CORRECTION
    ###########################################

    if APPLY_TERRAIN_FLATTENING:
        pre = tf.slope_correction(pre, TERRAIN_FLATTENING_MODEL, DEM, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)
        post = tf.slope_correction(post, TERRAIN_FLATTENING_MODEL, DEM, TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER)
        print('Radiometric terrain normalization is completed')

    ###########################################
    # 5. OUTPUT
    ###########################################

    # Convert to dB if specified
    if FORMAT == 'DB':
        pre = pre.map(help.lin_to_db)
        post = post.map(help.lin_to_db)
     
    # Clip to ROI if specified
    if CLIP_TO_ROI:
        pre = pre.map(lambda image: image.clip(ROI))
        post = post.map(lambda image: image.clip(ROI))

    # Calculate SAR index changes
    collection = SAR_indices.change(pre, post, INDEX)
    index_pre = collection[0]
    index_post = collection[1]
    change = collection[2]

    # Create a map if specified
    if MAKE_MAP:
        Map = geemap.Map()
        palette = ['#800000', '#FF0000', '#FFA500', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF']
        vis_params = {
            'bands': [change.bandNames().getInfo()[0]],
            'min': -50,
            'max': 50,
            'palette': palette
        }
        Map.addLayer(index_pre, {}, name=index_pre.bandNames().getInfo()[0])
        Map.addLayer(index_post, {}, name=index_post.bandNames().getInfo()[0])
        Map.addLayer(change, vis_params, name=change.bandNames().getInfo()[0])

    # Save the results to assets if specified
    if SAVE_ASSET:
        names = [change.bandNames().getInfo()[0], index_post.bandNames().getInfo()[0], index_pre.bandNames().getInfo()[0]]
        number = 0
        for i in collection:
            description = Filename
            assetId = 'projects/ee-thvdamme/assets/' + Filename
            task = ee.batch.Export.image.toDrive(image=i,
                                                 folder='RVI',
                                                 description=description + '_' + names[number],
                                                 region=ROI.getInfo()['coordinates'],
                                                 scale=10,
                                                 maxPixels=1e13)
            task.start()
            number += 1
            print('Exporting ' + names[number - 1])

    return Map
