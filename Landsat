var faith =  0.3
var ndvith = 0.15
var shp = table//矢量
var oeel=require('users/OEEL/lib:loadAll')
;//指定目标范围
////////////////////裁剪功能
var imgname = "hehu210623"
Map.centerObject(shp,10);
//////////////////////////////////////去云
//////////////////////////////////////
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;

  // Get the pixel QA band.
  var qa = image.select('pixel_qa');

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));

  // Return the masked image, scaled to reflectance, without the QA bands.
  return image
      .updateMask(mask)//去云选项，不用就注释掉
      .divide(10000)
      .select("B[0-9]*")
      .copyProperties(image, ["system:time_start"]);
}
//////////////////
// 影像数据集
var collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    //.filterDate('2019-09-20', '2019-09-22')
    //.filterDate('2020-08-20', '2020-08-25')
    //.filterDate('2020-09-07', '2020-09-10')
    //.filterDate('2018-07-25', '2018-08-01')
    .filterDate('2019-10-21', '2019-10-25')
    //.filter(ee.Filter.lte('CLOUD_COVER',100))
////////////////////////////////////////////
var expt='2018_nyh';//导出文件名称
////////////////////////////影像合成与显示
var image = collection
    .map(maskL8sr)
    .mosaic()
    .clip(shp);
print(image)
Map.addLayer(image.select(['B5', 'B4', 'B3', 'B2']),visParamimage,"image");
///////////////////////////////////////////////////

function bandmathfai(img) {
 var NIR = img.select("B5");
 var RED = img.select("B4");
 var BLUE = img.select("B2");
 var SNIR = img.select("B6");
 var result = img.expression(
   "A-B-(C-B)*((0.855-0.655)/(1.61-0.655))",
   {
     "A": NIR,
     "B": RED,
     "C": SNIR,
   }
 );
 return result;
}//

/////////////////////////////////////////
var svsimath=function(image){
var coefficients = ee.Array([
  [0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872],
  [-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608],
  [0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559],

]);
// Make an Array Image, with a 1-D Array per pixel.
var arrayImage1D = image.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7']).toArray();

// Make an Array Image with a 2-D Array per pixel, 6x1.
var arrayImage2D = arrayImage1D.toArray(1);

// Do a matrix multiplication: 6x6 times 6x1.
var componentsImage = ee.Image(coefficients)
  .matrixMultiply(arrayImage2D)
  // Get rid of the extra dimensions.
  .arrayProject([0])
  .arrayFlatten(
    [['brightness', 'greenness', 'wetness']]);
return componentsImage;
};
//TC变化

///////////////////////////
////////////////////////////////////影像显示参数
var visParamimage = {
 min: 0.0012,
 max: 0.12,
};
//////////////////////////////////////

Map.addLayer(shp,{palette: "ff0000"},"region",false);
/////////////////////////////////////////////////

var visParamndvi = {min: -0.36, max: 0.52, palette: [
"060eff","0a8eff","0cffb0","00ff1f","36ff04","ceff08",
"ffeb00","ff0000"
]};
var visParamsnir = {min: 0.022, max: 0.028, palette: [
"00ffff","0008ff"
]};
var visParamfai = {min: -0.03, max: 0.05, palette: [
"060eff","0a8eff","0cffb0","00ff1f","36ff04","ceff08",
"ffeb00","ff0000"
]};
var visParamatc3 = {min: 0.028, max: 0.09, palette: [
"060eff","0a8eff","0cffb0","00ff1f","36ff04","ceff08",
"ffeb00","ff0000"
]};
var visParamtcart = {
 min: 0,
 max: 3,
 palette: '100cff,#ff0000,10ff00,#f3ff35'
};
////////////////////
var find_peak = function(image,name1,name2){
    var histogram = image.reduceRegion({
    reducer: ee.Reducer.histogram(), 
    geometry: shp, 
    scale: 10,
    maxPixels: 1e13,
    tileScale: 10
    });
    //print(histogram)
    var threshold = otsu(histogram.get(name1));
    var mask = image.gte(threshold);
    var bufen = mask.rename("tiqu");

    var prue_water = image.updateMask(mask)
    var peak = prue_water.reduceRegion({
    reducer: ee.Reducer.mode(),
    geometry: shp,
    scale: 10,
    maxPixels: 1e13
     });
    var peak_value1 = peak.get(name1)
    var peak_value2 = ee.Number(peak_value1)
  //print(peak_value2 )
    var mask_pure = image.gte(peak_value2);
    var bufen2 = mask_pure.rename("tiqu2");//.updateMask(mask)
    return  peak_value2
  }
  var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);

  var indices = ee.List.sequence(1, size);

  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });

  return means.sort(bss).get([-1]);
};
/////////////////////////////////////////////////////////
var snir = image.select("B6");
Map.addLayer(snir,visParamsnir,"snir");
var TCT = svsimath(image);
var TC3 = TCT.select(['wetness']);
var ndvi = image.normalizedDifference(["B5","B4"]);
var ndwi = image.normalizedDifference(["B3","B5"]);
var fai=bandmathfai(image);
var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(), 
    geometry: shp, 
    scale: 30,
    maxPixels: 1e13,
    tileScale: 9
  });
print(histogram)
var ndwichart = ui.Chart.image.histogram(ndwi,shp,30).setChartType('LineChart')
print(ndwichart)
Map.addLayer(ndvi,visParamndvi,"ndvi");
Map.addLayer(fai,visParamfai,"fai");
Map.addLayer(ndwi,visParamndvi,"ndwi");
////////////////////////////////////////////////////////////////
var ndwi_find_peak = find_peak(ndwi,"nd",'no_fliter')
var ndwi_threshold = ee.Number(ndwi_find_peak)
var ndwimask = ndwi.select(["nd"]).gte(ndwi_threshold);//NDWI提取纯水
var faicut = fai.select(["B5"]).updateMask(ndwimask);
Map.addLayer(faicut,visParamfai,"faicut");
var faimean = faicut.reduceRegion({
    reducer: ee.Reducer.mean(), 
    geometry: shp, 
    scale: 30,
    maxPixels: 1e13,
    tileScale: 9
  }).get("B5");
var purewater_fai = ee.Number(faimean).multiply(0.95);
var vegetataion_fai = ee.Number(faith).multiply(0.05);
var faimean_num = purewater_fai.add(vegetataion_fai)
print(faimean_num)
var faimask = fai.select(["B5"]).lte(faimean_num);
var faimask1 = fai.select(["B5"]).gte(faimean_num);
var snircut = snir.select(["B6"]).updateMask(faimask);
Map.addLayer(snircut,visParamsnir,"snircut");
///////////////////////////////////////////////////
var bufen = faimask1.rename("nd");
var FEAV =ndvi.gte(ndvith);
var AVG  =bufen.add(FEAV);
var AVG  =AVG.where(AVG.eq(1),0);
var AVG  =AVG.where(AVG.eq(2),1);
var AVG  =bufen.add(AVG);
Map.addLayer(AVG,visParamtcart,"AVG");
//Export.image.toDrive({
//  image: image.select(['B5', 'B4', 'B3', 'B2']),
//  description: imgname+"image",//更改名字
//  fileNamePrefix: imgname+"image",//更改名字
//  region: cutregion,
//  scale: 30,
//  maxPixels: 1e13
//});
Export.image.toDrive({
  image: AVG.add(1),
  description: imgname+"_feng",//更改名字
  fileNamePrefix: imgname+"feng",//更改名字
  region: shp,
  scale: 30,
  maxPixels: 1e13
});

