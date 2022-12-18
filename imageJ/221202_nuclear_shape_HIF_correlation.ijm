// This does NOT correct for background signal!

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".oir") suffix

// Main
processFolder(input);
print("FINISHED");

// FUNCTIONS ----

function processFolder(input) {
  list = getFileList(input);
  list = Array.sort(list);
  for (i = 0; i < list.length; i++) {
    if(File.isDirectory(input + File.separator + list[i]))
      processFolder(input + File.separator + list[i]);
    if(endsWith(list[i], suffix))
      processFile(input, list[i]);
  }
}

function processFile(input, file) { 
  // Open image and do Z-projection
  BFopen(input, file);
  imageTitle = getTitle();
  imageName = File.nameWithoutExtension;
  selectWindow(imageTitle);
  run("Z Project...", "projection=[Max Intensity]");
  close(imageTitle);
  rename(imageTitle);
  
  // Create RGB duplicate for validation of segmentation
  Stack.setChannel(1);
  run("Duplicate...", "title=validation");
  run("Enhance Contrast", "saturated=0.35");
  setOption("ScaleConversions", true);
  run("RGB Color");
  
  // Segment nuclei and filter by size
  selectWindow(imageTitle);
  Stack.setChannel(1);
  run("Duplicate...", "title=nuclei");
  
  run("8-bit");
  run("Gaussian Blur...", "sigma=2");
  run("adaptiveThr ", "using=Mean from=45 then=-10");
  run("Convert to Mask");
  run("Distance Transform Watershed", "distances=[City-Block (1,2)] output=[32 bits] normalize dynamic=1 connectivity=8");
  run("8-bit");
  setThreshold(1, 255, "raw");
  run("Convert to Mask");
  close("nuclei");
  selectWindow("nuclei-dist-watershed");
  
  run("Analyze Particles...", "size=20-Infinity show=Overlay add");
  close("nuclei-dist-watershed");
  
  // Measure nuclei and get shape descriptors, then save results
  selectWindow(imageTitle);
  Stack.setChannel(2);
  run("Duplicate...", "title=hif");
  selectWindow("hif");
  //run("8-bit");
  run("Set Measurements...", "area mean modal shape feret's integrated median redirect=None decimal=2");
  roiManager("Measure");
  close("hif");
  
  saveAs("Results", output + File.separator + imageName + ".csv");
  
  // Draw outlines on validation image and save 
  selectWindow("validation");
  setForegroundColor(255, 0, 255);
  setLineWidth(1);
  roiManager("Draw");
  
  saveAs("Tiff", output + File.separator + imageName + ".tif");
  
  // Clean up
  roiManager("reset");
  run("Clear Results");
  close("*");

}

function BFopen(input, file) { 
  // Open input using the bioformats importer
  run("Bio-Formats Importer", 
  "open=[" + input + File.separator + file + 
  "] autoscale color_mode=Default rois_import=[ROI manager]" +
  " view=Hyperstack stack_order=XYCZT");
}
