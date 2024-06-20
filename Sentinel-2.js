Landsat 版本
var ndvith = 0.25
var faith = 0.25
var ID = table //矢量
var oeel=require('users/OEEL/lib:loadAll')
//var cutregion=geometry;//指定目标范围
////////////////////裁剪功能
var cutimg = function(image) {
  var imagecut =image.clip(cutregion);
 return imagecut;
};
var imgname = "hehu210623"
Map.centerObject(ID,12);
//////////////////////////////////////去云
//////////////////////////////////////

function maskS2clouds(image) {
//选择有关云掩膜的波段
  var qa = image.select('QA60')
  //位10和11分别代表云和卷云。
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // 将有关云的像元置为0
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0))
  // 掩模数据并选择多光谱波段，将反射率归为0-1复制时间属性
  return image.updateMask(mask)
      .divide(10000)
      .select("B.*")
      .copyProperties(image, ["system:time_start"])
}
//////////////////
// 影像数据集
var collection  = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
//var collection  = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
   //.filterDate('2021-11-13', '2021-11-16')
   //.filterDate('2020-09-01', '2020-09-06')
   .filterDate('2019-10-19', '2019-10-25')
    //.filter(ee.Filter.lte('CLOUD_COVER',100))
    //.map(cutimg);
////////////////////////////////////////////
var expt='2018_nyh';//导出文件名称
////////////////////////////影像合成与显示
var image = collection
    .map(maskS2clouds)
    .mosaic()
    .clip(ID);
///////////////////////////////////////////////////

function bandmathfai(img) {
 var NIR = img.select("B8");
 var RED = img.select("B4");
 var BLUE = img.select("B2");
 var SNIR = img.select("B11");
 var result = img.expression(
   "A-B-(C-B)*((0.855-0.655)/(1.61-0.655))",
   {
     "A": NIR,
     "B": RED,
     "C": SNIR,
   }
 ).rename('B5');
 return result;
}//

/////////////////////////////////////////
//TC变化

///////////////////////////
////////////////////////////////////影像显示参数
var visParamimage = {
 min: 0.0012,
 max: 0.12,
};
//////////////////////////////////////

Map.addLayer(image.select(['B8', 'B4', 'B3', 'B2']),visParamimage,"image",false);
Map.addLayer(ID,{palette: "ff0000"},"region",false);
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
/////////////////////////////////////////////////////////
//var snir = image.select("B11");
//Map.addLayer(snir,visParamsnir,"snir");
var ndvi = image.normalizedDifference(["B8","B4"]);
var ndwi = image.normalizedDifference(["B3","B8"]);
var fai=bandmathfai(image);
var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(), 
    geometry: ID, 
    scale: 10,
    maxPixels: 1e13,
    tileScale: 9
  });
print(histogram)
var ndwichart = ui.Chart.image.histogram(ndwi,ID,30).setChartType('LineChart')
print(ndwichart)
//Map.addLayer(ndvi,visParamndvi,"ndvi");
//Map.addLayer(fai,visParamfai,"fai");
//Map.addLayer(ndwi,visParamndvi,"ndwi");
////////////////////////////////////////////////////////////////
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
var find_peak = function(image,name1,name2){
    var histogram = image.reduceRegion({
    reducer: ee.Reducer.histogram(), 
    geometry: ID, 
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
    geometry: ID,
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

//var ndwimask = ndwi.select(["nd"]).gte(ndwith);
var ndwi_find_peak = find_peak(ndwi,"nd",'no_fliter')
var ndwi_threshold = ee.Number(ndwi_find_peak)
var ndwimask = ndwi.select(["nd"]).gte(ndwi_threshold);//NDWI提取纯水
//////////////////
var faicut = fai.select(["B5"]).updateMask(ndwimask);
Map.addLayer(faicut,visParamfai,"faicut");
var faimean = faicut.reduceRegion({
    reducer: ee.Reducer.mean(), 
    geometry: ID, 
    scale: 10,
    maxPixels: 1e13,
    tileScale: 9
  }).get("B5");
var purewater_fai = ee.Number(faimean).multiply(0.95);
var vegetataion_fai = ee.Number(faith).multiply(0.05);
//var vegetataion_fai = ee.Number(0.2535).multiply(0.2);
var faimean_num = purewater_fai.add(vegetataion_fai)
print(faimean_num)
var faimask = fai.select(["B5"]).lte(faimean_num);
var faimask1 = fai.select(["B5"]).gte(faimean_num);
//var snircut = snir.select(["B6"]).updateMask(faimask);
//Map.addLayer(snircut,visParamsnir,"snircut");
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
  region: ID,
  scale: 30,
  maxPixels: 1e13
});
