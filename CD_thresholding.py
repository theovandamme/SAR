import ee
import Geo_assets as Ga
import geemap
ee.Initialize()

Map = geemap.Map()
merit = ee.Image("MERIT/Hydro/v1_0_1")

visparam = {'bands':['VV','VH','VV/VH'],"min":[-25,-25,0],"max":[0,0,2]}
RgbVis = {"min":[-25,-25,0],"max":[0,0,2]}
Red = {'min': 0, 'max': 1, 'palette': ['FF0000']}
Yellow = {'min': 0, 'max': 1, 'palette': ['FFFF00']}
def add_VV_VH_band(image):
    ratioBand = image.select('VV').divide(image.select('VH')).rename('VV/VH')
    Rgb_ref = image.addBands(ratioBand)
    return Rgb_ref

def difference(Rgb_ref, Rgb_cris):
    diff_S = Rgb_cris.subtract(Rgb_ref)
    return diff_S

def calculate_stat(diff_S):
    reducers = ee.Reducer.mean().combine(
        reducer2 = ee.Reducer.stdDev(),sharedInputs= True)
    stats_S = diff_S.reduceRegion(
        reducer = reducers,
        scale= 10,
        bestEffort= True,
        tileScale= 16)
    return stats_S

def define_threshold(stats_S, diff_S, hand, hand_threshold, CD): 
    mu_VH = ee.Number(stats_S.get('VH_mean'))
    stddev_VH = ee.Number(stats_S.get('VH_stdDev'))
    mu_VV = ee.Number(stats_S.get('VV_mean'))
    stddev_VV = ee.Number(stats_S.get('VV_stdDev'))
    variablesVH = {'x': mu_VH, 'y': stddev_VH}
    #threshold_VH = ee.Number.expression('x - (1.5 * y)', variablesVH)
    threshold_VH = ee.Number(mu_VH).subtract(ee.Number(stddev_VH).multiply(CD))
    variablesVV = {'x': mu_VV, 'y': stddev_VV}
    #threshold_VV = ee.Number.expression('x - (1.5 * y)', variablesVV)
    threshold_VV = ee.Number(mu_VV).subtract(ee.Number(stddev_VV).multiply(CD))
    flood_VH = diff_S.select('VH').lt(threshold_VH).selfMask().rename("VH_CDAT")
    flood_VV = diff_S.select('VV').lt(threshold_VV).selfMask().rename("VV_CDAT")
    flood_VV = flood_VV.updateMask(hand.lt(hand_threshold))
    flood_VH = flood_VH.updateMask(hand.lt(hand_threshold))
    return (flood_VH, flood_VV)

def select_pixels(water_cris, Rgb_cris, flood_VH, flood_VV, hand, hand_threshold):
    LowBS_cris_VH = Rgb_cris.select('VH').lt(ee.Number(water_cris.get('threshold_VH')).add(ee.Number(1.5))).selfMask().rename("VH_water")
    LowBS_cris_VV = Rgb_cris.select('VV').lt(ee.Number(water_cris.get('threshold_VV')).add(ee.Number(1.5))).selfMask().rename("VV_water")
    
    # Retain only the pixels which fulfill both conditions/thresholds
    flood_VH = flood_VH.multiply(LowBS_cris_VH)
    flood_VV = flood_VV.multiply(LowBS_cris_VV)
    
    # Post-processing 
    
    # Post processing: mask out areas with too high HAND
    flood_VH = flood_VH.updateMask(hand.lt(hand_threshold))
    flood_VV = flood_VV.updateMask(hand.lt(hand_threshold))
    return (flood_VH, flood_VV)


def label_and_eliminate(flood_VH, flood_VV, MMU):
    # Eliminate small groups of pixels (MMU >= 10)  
# Uniquely label the image objects.
    objectId_VH = flood_VH.connectedComponents(
      connectedness=ee.Kernel.square(1),
      maxSize=50
    )
    objectId_VV = flood_VV.connectedComponents(
      connectedness=ee.Kernel.square(1),
      maxSize=50
    )
    objectSize_VH = objectId_VH.select('labels').connectedPixelCount(
      maxSize=50, eightConnected=True
    )
    objectSize_VV = objectId_VV.select('labels').connectedPixelCount(
      maxSize=50, eightConnected=True
    )
    sizeMaskVH = objectSize_VH.lt(MMU)
    sizeMaskVV = objectSize_VV.lt(MMU)
    flood_VH = flood_VH.updateMask(sizeMaskVH.unmask().Not())
    flood_VV = flood_VV.updateMask(sizeMaskVV.unmask().Not())
    return (flood_VH, flood_VV)

def ComputeWater(image, hand, hand_threshold):
    histogram = image.reduceRegion(
        reducer = ee.Reducer.histogram(255, 0.1).combine(reducer2 = 'mean',outputPrefix = None, sharedInputs=True).combine(reducer2 = 'variance', outputPrefix = None, sharedInputs=True),
        geometry = Ga.Geometry_otsu,
        scale =10,
        bestEffort = True)
    threshold_VV = otsu(histogram.get('VV_histogram')).add(ee.Number(0))
    threshold_VH = otsu(histogram.get('VH_histogram')).add(ee.Number(0))
 
    waterVV = image.select('VV').lt(threshold_VV).selfMask().rename("VV_water")
    waterVH = image.select('VH').lt(threshold_VH).selfMask().rename("VH_water")
  
    waterVV = waterVV.updateMask(hand.lt(hand_threshold))
    waterVH = waterVH.updateMask(hand.lt(hand_threshold))

    waterLayer = waterVV.addBands(waterVH)
  
    waterLayer = waterLayer.set('threshold_VV', threshold_VV)
    waterLayer = waterLayer.set('threshold_VH', threshold_VH)
    return waterLayer

def otsu(histogram):
    counts = ee.Array(ee.Dictionary(histogram).get("histogram"))
    means = ee.Array(ee.Dictionary(histogram).get("bucketMeans"))
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean = sum.divide(total)

    indices = ee.List.sequence(1, size)

    # Compute between sum of squares, where each mean partitions the data.

    def func_xxx(i):
        aCounts = counts.slice(0, 0, i)
        aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
        aMeans = means.slice(0, 0, i)
        aMean = aMeans.multiply(aCounts).reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount)
        bCount = total.subtract(aCount)
        bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount)
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)))

    bss = indices.map(func_xxx)

    # Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1])

def water_ref_otsu(water_ref):
    objectId_refVV = water_ref.select('VV_water').connectedComponents(
      connectedness=ee.Kernel.square(1),
      maxSize=50
    )
    objectId_refVH = water_ref.select('VH_water').connectedComponents(
      connectedness=ee.Kernel.square(1),
      maxSize=50
    )
    objectSize_refVV = objectId_refVV.select('labels').connectedPixelCount(
      maxSize=50, eightConnected=True
    )
    objectSize_refVH = objectId_refVH.select('labels').connectedPixelCount(
      maxSize=50, eightConnected=True
    )
    sizeMaskVV = objectSize_refVV.lt(50)
    sizeMaskVH = objectSize_refVH.lt(50)
    water_ref_maskedVV = water_ref.select('VV_water').updateMask(sizeMaskVV.unmask().Not())
    water_ref_maskedVH = water_ref.select('VH_water').updateMask(sizeMaskVH.unmask().Not())
    water_ref_masked = water_ref_maskedVV.addBands(water_ref_maskedVH)
    return water_ref_masked 
  
def CD_Thresholding(img_ref, img_crisis, date, hand, MMU, CD):
    hand_threshold = hand
    hand = merit.select('hnd').clip(Ga.BujxUvira)
    hand_masked = hand.updateMask(hand.lt(hand_threshold))
    Rgb_ref = add_VV_VH_band(img_ref)
    Rgb_cris = add_VV_VH_band(img_crisis)
    diff_S = difference(Rgb_ref, Rgb_cris)
    stats_S = calculate_stat(diff_S)
    flood_VH = define_threshold(stats_S, diff_S, hand, hand_threshold, CD)[0]
    flood_VV = define_threshold(stats_S, diff_S, hand, hand_threshold, CD)[1]
    water_cris = ComputeWater(Rgb_cris, hand, hand_threshold)
    flood_VH = select_pixels(water_cris, Rgb_cris, flood_VH, flood_VV, hand, hand_threshold)[0]
    flood_VV = select_pixels(water_cris, Rgb_cris, flood_VH, flood_VV, hand, hand_threshold)[1]
    flood_VH = label_and_eliminate(flood_VH, flood_VV, MMU)[0]
    flood_VV = label_and_eliminate(flood_VH, flood_VV, MMU)[1]
    water_ref = ComputeWater(Rgb_ref, hand, hand_threshold)
    water_ref_masked = water_ref_otsu(water_ref)
    flood_VH = flood_VH.updateMask(water_ref_masked.select('VH_water').unmask(0).Not()) #toegevoegde lijn
    flood_VV = flood_VV.updateMask(water_ref_masked.select('VV_water').unmask(0).Not())
    flood = flood_VH.addBands(flood_VV)
    task = ee.batch.Export.image.toDrive(image=Rgb_cris,
                                                 description= 'rgv' + date,
                                                 
                                                 region=Ga.BujxUvira,
                                                 scale=10,
                                                 maxPixels=1e13)
    #task.start()
    Map.centerObject(Ga.BujxUvira, 10)
    Map.addLayer(hand_masked, name = 'hand', shown = False)
    Map.addLayer(Rgb_ref, RgbVis, name ='RGB ref',shown = False)
    Map.addLayer(Rgb_cris, RgbVis, name ='Crisis RGB',shown = False)
    Map.addLayer(diff_S.select('VV'),{min: -12, max: 8}, name ='difference subtraction VV',shown = False)
    Map.addLayer(water_ref_masked.select('VV_water'), Yellow, name ='water layer ref image VV',shown = False)
    Map.addLayer(water_ref_masked.select('VH_water'), Yellow, name ='water layer ref image VH', shown =True)
    Map.addLayer(flood_VH, Red, name ='flooded CDAT VH', shown =True)
    Map.addLayer(flood_VV, Red, name ='flooded CDAT VV', shown =False)
    
    
    return Map
