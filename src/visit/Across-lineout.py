# Script for running velocity profiles in time. 
#
# Note: you must create destination folders before beginning! 
#
import sys
import gc
gc.enable()
WDin=sys.argv[1]
WDout=sys.argv[2]
startnum=int(float(sys.argv[3]))
endnum=int(float(sys.argv[4]))

for num in range(startnum,endnum):

	OpenDatabase(str(WDin)+"/viz_IB2d"+str(num)+"/dumps.visit", 0)
	AddPlot("Pseudocolor", "U_magnitude", 1, 1)
	DrawPlots()
	SetActivePlots((0, 1))
	SetActivePlots(0)
	DeleteActivePlots()
	SetTimeSliderState(5)
	SetQueryFloatFormat("%g")
	Query("Lineout", end_point=(1e-4, 20e-6, 0), num_samples=50, start_point=(-1e-4, 20e-6, 0), use_sampling=0)
	SetActiveWindow(2)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = str(WDout)+"/sim"+str(num)
	SaveWindowAtts.fileName = "Um_line20_"
	SaveWindowAtts.family = 1
	SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.screenCapture = 0
	SaveWindowAtts.saveTiled = 0
	SaveWindowAtts.quality = 80
	SaveWindowAtts.progressive = 0
	SaveWindowAtts.binary = 0
	SaveWindowAtts.stereo = 0
	SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
	SaveWindowAtts.forceMerge = 0
	SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.advancedMultiWindowSave = 0
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()
	DeleteActivePlots()
	
	SetActiveWindow(1)
	Query("Lineout", end_point=(1e-4, 30e-6, 0), num_samples=50, start_point=(-1e-4, 30e-6, 0), use_sampling=0)
	SetActiveWindow(2)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = str(WDout)+"/sim"+str(num)
	SaveWindowAtts.fileName = "Um_line30_"
	SaveWindowAtts.family = 1
	SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.screenCapture = 0
	SaveWindowAtts.saveTiled = 0
	SaveWindowAtts.quality = 80
	SaveWindowAtts.progressive = 0
	SaveWindowAtts.binary = 0
	SaveWindowAtts.stereo = 0
	SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
	SaveWindowAtts.forceMerge = 0
	SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.advancedMultiWindowSave = 0
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()
	DeleteActivePlots()
	
	SetActiveWindow(1)
	Query("Lineout", end_point=(1e-4, 50e-6, 0), num_samples=50, start_point=(-1e-4, 50e-6, 0), use_sampling=0)
	SetActiveWindow(2)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = str(WDout)+"/sim"+str(num)
	SaveWindowAtts.fileName = "Um_line50_"
	SaveWindowAtts.family = 1
	SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.screenCapture = 0
	SaveWindowAtts.saveTiled = 0
	SaveWindowAtts.quality = 80
	SaveWindowAtts.progressive = 0
	SaveWindowAtts.binary = 0
	SaveWindowAtts.stereo = 0
	SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
	SaveWindowAtts.forceMerge = 0
	SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.advancedMultiWindowSave = 0
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()
	DeleteActivePlots()
	SetActiveWindow(1)
	DeleteActivePlots()
	DeleteActivePlots()
	
exit()
