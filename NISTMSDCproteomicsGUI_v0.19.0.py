import subprocess, os, re, time
from multiprocessing import Queue, Process, cpu_count
from sqlitedict import SqliteDict
import tkinter as Tk
import tkinter.ttk as ttk
from tkinter import filedialog

#import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

#print("Loading data...")

#msgpack_file = open("tandemNIST17_gui.msgpack", 'rb')
#NIST17data = msgpack.unpack(msgpack_file, use_list=False, encoding='utf-8')
#print("Complete!")

#print(list(NISTlibdata))

#precPatt = re.compile(r'PrecursorMZ:\s(\d+\.\d+)')
#spectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)')

# this is MSP format stuff
namePatt = re.compile(r'Name:\s(.+)')
precPatt = re.compile(r'PrecursorMZ:\s(\d+\.\d+)')
spectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.?\d*)')

# this is for the MSConvert MGF file format
MGFnamePatt = re.compile(r'TITLE=(.+)')
MGFprecPatt = re.compile(r'PEPMASS=(\d+\.\d+)\s')
MGFspectrumPatt = re.compile(r'(\d{2,4}\.\d{3,8})\s(\d+\.\d+)')


processedResults = dict()
theOverPlot = []
theRange = 0

# this is the path for the mspepsearch library
#libPath = os.path.join("C:/", "LIBS", "human_hair_selected_keratins_and_related_updated")
zppm = 20

mspepsearchInputType = "msp" # input defaults to MSP


f = Figure(figsize=(10, 6), dpi=100, facecolor='white')
a = f.add_subplot(111)
f.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
a.stem([0], [0], markerfmt=' ', use_line_collection=True)
a.spines['top'].set_visible(False)
a.spines['right'].set_visible(False)


def batchSearcherMGF(inputFileName, libPath, libType, aQu):

	outFileName = inputFileName.split(".")[0] + "Results.tsv"

	if libType == "hiRes":
		resParameters = ["/MPPM", "50"]
	elif libType == "loRes":
		resParameters = ["/ZPPM", "20", "/M", "0.4"]
	else:
		print("LIBTYPE IMPROPERLY SET!!!!!!") # this is really only for diagnostics

	theCallBATCH = ["MSPepSearch64.exe", "iPval", # ydiPva
		"/ZI", "2.5"] + resParameters + ["/MzLimits", "50", "-1", 
		"/HITS", "1", 
		"/OutBestHitsOnly"]

	for aLib in libPath:
		theCallBATCH.extend(["/LIB", aLib])

	theCallBATCH.extend([ "/INP", # /LibInMem
		inputFileName, # this is the input file
		"/OUTTAB", outFileName, "/OutSpecNum", "/OutPrecursorMZ",
		"/OutNumComp", "/OutNumMP","/OutMaxScore", #"/OutPepLocalization",
		"/OutCE", "/MinMF", "200"])

	subprocess.run(theCallBATCH)

	batchResults = open(outFileName, "r")
	print("...batch processing complete!!!")

	qualifiedIDset = []
	lineCount = 0
	for line in batchResults:
		parsedstuff = line.split("\t")

		if (lineCount >= 4) & (len(parsedstuff) > 5):

			theSeq = parsedstuff[14] # 19
			theMod = parsedstuff[18] if parsedstuff[18] != "" else "0" # 23
			theCharge = parsedstuff[15] # 20
			theCollision = parsedstuff[16] # 21

			uniqueID = theSeq + "_" + theMod + "_" + theCharge + "_" + theCollision

			theTitle = parsedstuff[1].split()[0]
			libraryID = parsedstuff[4]
			MFscore = parsedstuff[9] # this needs to be double checked


			qualifiedIDset.append([theTitle, libraryID, MFscore, uniqueID])

		lineCount += 1


	batchResults.close()

	aQu.put(qualifiedIDset)


def rawProcessor(_TEMPinputSpectralData, _inputSpectralData, libPath, libType):
	print("Starting batch processing for %s spectra..." % (len(_TEMPinputSpectralData)))

	letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"]
	#SMPcount = len(letters)
	SMPcount = cpu_count() // 2 # this will always take half the number of logical cores available in the machine, rounded down if for some reason it's odd

	'''
	batchInputDataA = open('inputFileBATCH_A.mgf', 'w')
	batchInputDataB = open('inputFileBATCH_B.mgf', 'w')
	batchInputDataC = open('inputFileBATCH_C.mgf', 'w')
	batchInputDataD = open('inputFileBATCH_D.mgf', 'w')
	'''

	batchInputData = [[open('inputFileBATCH_'+aLetter+'.mgf', 'w'), 'inputFileBATCH_'+aLetter+'.mgf'] for aLetter in letters]

	'''
	batchInputData = [[batchInputDataA, 'inputFileBATCH_A.mgf'], [batchInputDataB, 'inputFileBATCH_B.mgf'], 
						[batchInputDataC, 'inputFileBATCH_C.mgf'], [batchInputDataD, 'inputFileBATCH_D.mgf']]
	'''

	for i, selectedQuery in enumerate(list(_TEMPinputSpectralData)):

		queryData = _TEMPinputSpectralData[selectedQuery]

		inputText = queryData[2]
		inputSpectrum = queryData[0]
		precursorMZ = queryData[1]

		batchInputData[i % SMPcount][0].write(inputText)

	#batchInputDataA.close(); batchInputDataB.close()
	[_batchinput[0].close() for _batchinput in batchInputData]

	resultsQu = Queue()

	#procA = Process(target = batchSearcherMGF, args=('inputFileBATCH_A.mgf', libPath, resultsQu))
	#procB = Process(target = batchSearcherMGF, args=('inputFileBATCH_B.mgf', libPath, resultsQu))
	#theProcs = [procA, procB]
	theProcs = [Process(target = batchSearcherMGF, args=(batchInputData[i][1], libPath, libType, resultsQu)) for i in range(0, SMPcount)]

	[aProc.start() for aProc in theProcs]

	itsAlive = any([aProc.is_alive() for aProc in theProcs])

	qualifiedIDset = []

	while (itsAlive | (not resultsQu.empty())):
		if not resultsQu.empty():
			partialResult = resultsQu.get()

			qualifiedIDset.extend(partialResult)

		time.sleep(0.2)
		itsAlive = any([aProc.is_alive() for aProc in theProcs])

	[aProc.join() for aProc in theProcs]
	resultsQu.close()

	print("%s spectra found to have search results" % (len(qualifiedIDset)))
	
	for qualifiedData in qualifiedIDset:
		qualifiedTitle = qualifiedData[0]
		libID = qualifiedData[1]
		MFscore = qualifiedData[2]
		uniqueID = qualifiedData[3]

		_inputSpectralData[qualifiedTitle] = _TEMPinputSpectralData[qualifiedTitle] + [libID, int(MFscore), uniqueID]


def MSCONVERTextraction(filename, libPath, libType, qu):
	theCommand = ["pwiz-bin-windows-x86_64\msconvert.exe", filename,
				"--filter", "msLevel 2",
				"--filter", "zeroSamples removeExtra", "--mgf"]

	#result = subprocess.run(theCommand, stdout=subprocess.PIPE)
	result = subprocess.run(theCommand)

	justFileName = filename.split("/")[-1].split(".")[0]
	inputSpectrumReader = open(justFileName+".mgf", "r")

	_TEMPinputSpectralData = dict()
	_inputSpectralData = dict()

	inputSpectrum = []; inputAllText = ""
	precursorMZ = None; precursorAb = 0
	theMaxAb = 0; totalAb = 0

	_count = 0; _validSpecCount = 0
	for line in inputSpectrumReader:
		if line.startswith("TITLE"):
			theName = re.search(MGFnamePatt, line).group(1)

		spectrumData = re.match(MGFspectrumPatt, line)
		if spectrumData is not None:
			theMZ = float(spectrumData.group(1))
			theAb = float(spectrumData.group(2))
			if theAb > 0:
				totalAb += theAb
				inputSpectrum.append([theMZ, theAb])

				if (theMZ < (precursorMZ + 1)):
					findThePrecursor_ppm = (abs(theMZ - precursorMZ) / precursorMZ) * 1000000

					if findThePrecursor_ppm < 20:
						precursorAb += theAb

				if theAb > theMaxAb:
					theMaxAb = theAb

			inputAllText += str(round(theMZ, 4)) + " " + str(round(theAb, 2)) + "\n"

		else:
			precData = re.match(MGFprecPatt, line)
			if precData is not None:
				precursorMZ = float(precData.group(1))

			inputAllText += line

		if line.rstrip()=="END IONS":

			if totalAb > 0:
				if precursorAb / totalAb <= 0.50:

					inputSpectrum = [[peak[0], 1000*(peak[1]/theMaxAb)] for peak in inputSpectrum]
					_TEMPinputSpectralData[theName] = [inputSpectrum, precursorMZ, inputAllText]

					_validSpecCount += 1
				
			inputSpectrum = []; inputAllText = ""; theMaxAb = 0; totalAb = 0
			precursorMZ = None; precursorAb = 0

			_count += 1
			if _count % 50 == 0: qu.put(["justCount", _count, _validSpecCount, len(_inputSpectralData)])



			if len(_TEMPinputSpectralData) >= 8000:
				#time.sleep(30)
				
				rawProcessor(_TEMPinputSpectralData, _inputSpectralData, libPath, libType)
				_TEMPinputSpectralData.clear()

				qu.put(["justCount", _count, _validSpecCount, len(_inputSpectralData)])
				
	inputSpectrumReader.close()

	rawProcessor(_TEMPinputSpectralData, _inputSpectralData, libPath, libType)
	_TEMPinputSpectralData.clear()

	qu.put(["justCount", _count, _validSpecCount, len(_inputSpectralData)])
	qu.put(["finalData", justFileName, _inputSpectralData, _count])


def dereplicator(inputSpectralData_raw):
	compileByID = dict()

	#print(inputSpectralData_raw)

	for theName in inputSpectralData_raw:
		#_prec = inputSpectralData_raw[theName][1]
		_libID = inputSpectralData_raw[theName][3]
		_MFscore = inputSpectralData_raw[theName][4]
		_uniqID = inputSpectralData_raw[theName][5]

		trueUniq = (_uniqID, _libID)
		if trueUniq in compileByID:
			compileByID[trueUniq].append([_MFscore, theName])
		else:
			compileByID[trueUniq] = [[_MFscore, theName]]

	inputSpectralData_derep = dict()
	for trueUniq in compileByID:
		currentSpecs = compileByID[trueUniq]

		if len(currentSpecs) > 1:
			maxSpec = max(currentSpecs)
		else:
			maxSpec = currentSpecs[0]

		inputSpectralData_derep[maxSpec[1]] = inputSpectralData_raw[maxSpec[1]]


	print("Original spectra: %s -> Dereplicated: %s" % (len(inputSpectralData_raw), len(inputSpectralData_derep)))
	return inputSpectralData_derep


def calculateFDR(inputSpectralData):
	# do the FDR calculation here
	compiledData = dict()
	for theName in inputSpectralData:
		_prec = inputSpectralData[theName][1]
		_libID = inputSpectralData[theName][3]
		_MFscore = inputSpectralData[theName][4]

		#compiledData.append([_MFscore, _libID])
		if _MFscore in compiledData:
			compiledData[_MFscore].append([_libID, _prec])
		else:
			compiledData[_MFscore] = [[_libID, _prec]]

	#compiledData = sorted(compiledData, key=lambda psig: psig[0], reverse=True) # sort the data high to low
	sortedMFs = sorted(list(compiledData), reverse=True)

	decoyCount = 0
	totalCount = 0
	finalMFcutoff = 200
	for i,anMF in enumerate(sortedMFs):
		totalCount += len(compiledData[anMF])
		for aResult in compiledData[anMF]:

			_lib = aResult[0]
			if "decoy" in _lib.lower():
				decoyCount += 1

		FDRcalculation = min((2*decoyCount)/(totalCount) , 1)

		if FDRcalculation < 0.01:
			print("FDR: %s MF: %s" % (FDRcalculation, anMF))
		else:
			print("FDR: %s MF: %s" % (FDRcalculation, anMF))
			finalMFcutoff = sortedMFs[i-1]
			break



	return finalMFcutoff # this is a dummy number for now


def rawExtract():
	global inputSpectralData, mspepsearchInputType, libPath, FDRbasedMFmin, libType

	rawFilePath = filedialog.askopenfilename()

	RAWwindow = Tk.Toplevel(root)
	RAWwindow.wm_title("Thermo RAW file extraction")  # Makes the title that will appear in the top left
	RAWwindow.geometry("750x400")
	RAWwindow.config(background="#FFFFFF")  # sets background color to white

	l = Tk.Label(RAWwindow, text="Extracting data from RAW file...")
	l.pack(side="top", fill="both", expand=True, padx=50, pady=50)

	#l.config(text = "Wait till I'm done...")

	l.update_idletasks()

	theQu = Queue()

	thread = Process(target = MSCONVERTextraction, args=(rawFilePath, libPath, libType, theQu))
	startTime = time.time()
	thread.start()

	tripwire = True

	cumulativeText = "Initial RAW extraction complete!!\n"

	while (thread.is_alive() | (not theQu.empty())):
		
		elapsedTime = round(time.time() - startTime)

		if tripwire:
			l.config(text = "MSConvert RAW extraction in progress (%ss elapsed)..." % elapsedTime)
			
		if not theQu.empty():
			#print("Something's in the line!!!")
			tripwire = False
			signals = theQu.get()

			if signals[0] == "justCount":
				rawSpecCount = signals[1]
				suffFragCount = signals[2]
				prelimCount = signals[3]
				l.config(text = cumulativeText + "Processing results...\n >> %s spectra found, %s sufficiently fragmented, %s with prelim. results" % (rawSpecCount, suffFragCount, prelimCount))
				
			elif signals[0] == "finalData":
				#inputSpectralData = signals[1]
				justFileName = signals[1]
				inputSpectralData_raw = signals[2]
				finalCount = signals[3]

		RAWwindow.update()
		time.sleep(0.05)

	thread.join()
	theQu.close()
	cumulativeText += "First-pass spectral search complete! <%s raw spectra -> %s sufficiently fragmented -> %s with hits>\n" % (rawSpecCount, suffFragCount, prelimCount)
	l.config(text = cumulativeText)

	inputSpectralData = dereplicator(inputSpectralData_raw)
	cumulativeText += "Spectral dereplication complete! <%s initial results -> %s dereplicated results>\n" % (len(inputSpectralData_raw), len(inputSpectralData))
	l.config(text = cumulativeText)

	FDRbasedMFmin = calculateFDR(inputSpectralData)
	cumulativeText += "FDR calculations complete! <FDR MF cutoff: %s>\n" % FDRbasedMFmin
	l.config(text = cumulativeText)

	toRemove = []
	for theName in inputSpectralData:
		_libID = inputSpectralData[theName][3]
		_MFscore = inputSpectralData[theName][4]

		if (_MFscore < FDRbasedMFmin) | (if "decoy" in _libID.lower()):
			toRemove.append(theName)

	for deleteName in toRemove:
		del inputSpectralData[deleteName]
	cumulativeText += "Removal of %s spectra w/ hits below FDR threshold or only decoy complete!\n" % (len(toRemove))
	l.config(text = cumulativeText)

	cumulativeText += "%s file extraction and analysis complete! %s valid MS2 spectra found!\n" % (justFileName+".raw", len(inputSpectralData))
	l.config(text = cumulativeText)
	#l.config(text = "%s extracted, %s valid MS2 spectra!" % (justFileName+".raw", len(inputSpectralData)))
	RAWwindow.update()

	right_tree.delete(*right_tree.get_children())
	for theName in inputSpectralData:
		_prec = inputSpectralData[theName][1]

		right_tree.insert("", "end", theName, text=theName, values=(_prec, "???"))

	mspepsearchInputType = "mgf"


def openMSPfile():
	global processedResults, theRange, inputSpectralData, mspepsearchInputType, FDRbasedMFmin
	#print("Hello!")
	file_path = filedialog.askopenfilename()
	inputSpectrumReader = open(file_path, "r")
	
	inputSpectralData = dict()
	inputSpectrum = []; inputAllText = ""
	precursorMZ = None

	_count = 0
	for line in inputSpectrumReader:
		if line.startswith("Name"):
			theName = re.search(namePatt, line).group(1)

		spectrumData = re.match(spectrumPatt, line)
		if spectrumData is not None:
			inputSpectrum.append([float(spectrumData.group(1)), float(spectrumData.group(2))])
			#_textOut += theMZ + " " + str(round(theAb, 2))
		else:
			precData = re.match(precPatt, line)
			if precData is not None:
				precursorMZ = float(precData.group(1))
		inputAllText += line

		if line.rstrip()=="":

			inputSpectralData[theName] = [inputSpectrum, precursorMZ, inputAllText]

			inputSpectrum = []; inputAllText = ""
			precursorMZ = None
			
			_count += 1

	# record the last entry---------
	if len(inputSpectrum) > 1:
		inputSpectralData[theName] = [inputSpectrum, precursorMZ, inputAllText]

		inputSpectrum = []; inputAllText = ""
		precursorMZ = None
		
		_count += 1
	# -----------------------------

	inputSpectrumReader.close()

	right_tree.delete(*right_tree.get_children())
	for theName in inputSpectralData:
		_prec = inputSpectralData[theName][1]

		right_tree.insert("", "end", theName, text=theName, values=(_prec, "???"))

	mspepsearchInputType = "msp"
	FDRbasedMFmin = 200


def itemClickedHitsTree(event):
	uniqueID = str(left_tree.focus())

	global processedResults, theOverPlot, theRange, xData, yData, sortedSpectrum

	resultAllText = ""
	#print(processedResults)
	deltaMass = processedResults[uniqueID][6]
	deltaMZ = float(deltaMass)

	nistSpectrum = NISTlibdata[uniqueID][3]
	theMW = NISTlibdata[uniqueID][0]
	theComment = NISTlibdata[uniqueID][2]
	nistPrecursorMZ = NISTlibdata[uniqueID][1]

	resultAllText += "%s\nPrecursor m/z: %s\nDelta m/z: %s\n\n" % (theComment, nistPrecursorMZ, deltaMass)

	#print(nistSpectrum)

	matchSpectrum = []
	maxIntensity = 0
	resultAllText += "m/z\tIntensity\n---------------\n"
	for peak in nistSpectrum:
		resultAllText += "%s\t%s\n" % (peak[0], peak[1])
		matchSpectrum.append([float(peak[0]), peak[1]])

		if peak[1] > maxIntensity: maxIntensity = peak[1]



	specDataR.delete(1.0, Tk.END)
	specDataR.insert(0.0, resultAllText)


	unmatchedSubsetSpectrum = []
	for queryPeak in sortedSpectrum:
		trueMZquery = queryPeak[0]

		matchBool = False
		for nistPeak in matchSpectrum:
			trueMZnist = nistPeak[0]
			check_ppmTrue = (abs(trueMZnist - trueMZquery) / min(trueMZnist, trueMZquery)) * 1000000

			if check_ppmTrue < zppm:
				matchBool = True
		
		if matchBool is False:
			unmatchedSubsetSpectrum.append(queryPeak)

	#print(unmatchedSubsetSpectrum)
	#print(matchSpectrum)

	shiftedSpectrum = [[(round(peak[0] + deltaMZ, 4), peak[0]), peak[1]] for peak in matchSpectrum]

	shiftedMatches = []; unshiftedSpectrum = []
	for nistPeak in shiftedSpectrum:
		trueMZnist = nistPeak[0][1]
		shiftMZnist = nistPeak[0][0]

		shiftMatchBool = False
		for queryPeak in unmatchedSubsetSpectrum:
			trueMZquery = queryPeak[0]

			check_ppmShifted = (abs(shiftMZnist - trueMZquery) / min(trueMZnist, trueMZquery)) * 1000000

			if check_ppmShifted < zppm:
				shiftMatchBool = True

				shiftedMatches.append([shiftMZnist, nistPeak[1]])

		if shiftMatchBool is False:
			unshiftedSpectrum.append([trueMZnist, nistPeak[1]])

	

	xData2 = []; yData2 = []
	for anAttr in unshiftedSpectrum:
		xData2.append(anAttr[0])
		yData2.append(-1 * (anAttr[1] / maxIntensity)*100)

	theRange = (max(xData + xData2) + 10) - (min(xData + xData2) - 10)
	a.set_xlim(min(xData + xData2) - 10, max(xData + xData2) + 10)

	xDataShift = []; yDataShift = []
	for anAttr in shiftedMatches:
		xDataShift.append(anAttr[0])
		yDataShift.append(-1 * (anAttr[1] / maxIntensity)*100)

	# --code for removing the old plots------
	if len(theOverPlot) > 0:

		oldPlotData = theOverPlot.pop()

		for oldPlotted in oldPlotData:
			if isinstance(oldPlotted, tuple):
				[subplot.remove() for subplot in oldPlotted]
			else:
				oldPlotted.remove()

		#del oldPlot
		del oldPlotData
	# --------------------------------------
	currentDrawnCollection = []
	if len(xData2) > 0:
		currentPlot = a.stem(xData2, yData2, 'red', markerfmt=' ', use_line_collection=True)
		currentPlot[2].remove()  # the baseline
		currentDrawnCollection.append((currentPlot[0], currentPlot[1]))

	if len(xDataShift) > 0:
		shiftPlot = a.stem(xDataShift, yDataShift, 'purple', markerfmt=' ', use_line_collection=True)
		shiftPlot[2].remove()
		currentDrawnCollection.append((shiftPlot[0], shiftPlot[1]))

	sortedSpectrum2 = sorted(unshiftedSpectrum + shiftedMatches)

	whatToPlot2 = []
	theRecords2 = []  # you use this to check if there is any peak that is in the range of your current peak that would get in the way of the text (i.e. it's a higher intensity)
	for mz_intensity in sortedSpectrum2:
		# this code will only allow text to be plotted either if the preceding plotted text is far enough away ( > theRange*0.05), or if it has a higher intensity than the preceding plotted text (which it subsequently removes)
		_xcoor = mz_intensity[0]
		_ycoor = (mz_intensity[1] / maxIntensity)*100


		if len(whatToPlot2) > 0:
			last_xcoor = whatToPlot2[-1][0]
			last_ycoor = whatToPlot2[-1][1]

			#print "%s --> %s" % (_xcoor, _xcoor - last_xcoor)

			if (_xcoor - last_xcoor < theRange*0.05) & (_ycoor > last_ycoor) & (_ycoor > 1):  # if it is NOT far enough away BUT it is taller than the last text then go ahead
				whatToPlot2.pop()
				whatToPlot2.append([_xcoor,_ycoor])
				#print "popped! %s" % (last_xcoor)
			elif (_xcoor - last_xcoor > theRange*0.05) & (_ycoor > 1):  # if it is far enough away from the last text then go ahead
				goAhead = True
				for checkPrevious in theRecords2:
					if (_xcoor - checkPrevious[0] < theRange*0.05) & (_ycoor < checkPrevious[1]):  # check to see if anything in range that is also greater than the current intensity
						goAhead = False
				if goAhead is True: whatToPlot2.append([_xcoor, _ycoor])
		else:
			if (_ycoor > 1): whatToPlot2.append([_xcoor, _ycoor])

		theRecords2.append([_xcoor, _ycoor])

	currentArrow = a.arrow(float(nistPrecursorMZ), -3, 0.0, 1, fc="k", ec="k", head_width=theRange * 0.005, head_length=2)
	currentDrawnCollection.append(currentArrow)


	for dataToPlot in whatToPlot2:
		aWord = a.annotate(str(dataToPlot[0]), xy=(dataToPlot[0], -1 * dataToPlot[1] - 5), ha='center', size='x-small')
		currentDrawnCollection.append(aWord)

	theOverPlot.append(currentDrawnCollection)
	canvas.draw()


def itemClickedQueryTree(event):
	global inputSpectralData, processedResults, xData, yData, sortedSpectrum, mspepsearchInputType, FDRbasedMFmin

	selectedQuery = str(right_tree.focus())

	queryData = inputSpectralData[selectedQuery]

	inputAllText = queryData[2]
	inputSpectrum = queryData[0]
	precursorMZ = queryData[1]

	specData.delete(1.0, Tk.END)
	specData.insert(0.0, inputAllText)

	specDataR.delete(1.0, Tk.END)
	specDataR.insert(0.0, "")

	inputSpectrum = sorted(inputSpectrum)

	sortedSpectrum = sorted(inputSpectrum)
	xData = []; yData = []

	for anAttr in sortedSpectrum:
		xData.append(anAttr[0])
		yData.append(anAttr[1] / 10)

	a.clear()
	lines1, stems1, baseline1 = a.stem(xData, yData, 'blue', markerfmt=' ', use_line_collection=True)
	baseline1.remove()

	theRange = (max(xData) + 10) - (min(xData) - 10)

	a.set_xlim(min(xData) - 10, max(xData) + 10)

	a.set_ylim(-110, 110)
	a.axvline(0, color='black', linewidth=2)
	a.axhline(0, color='black', linewidth=1)

	a.grid(color='lightgray', alpha=0.7)

	a.set_xlabel("m/z", fontsize=12)
	a.set_ylabel("Relative Intensity", fontsize=12)


	whatToPlot = []
	theRecords = []  # you use this to check if there is any peak that is in the range of your current peak that would get in the way of the text (i.e. it's a higher intensity)
	for mz_intensity in sortedSpectrum:
		# this code will only allow text to be plotted either if the preceding plotted text is far enough away ( > theRange*0.05), or if it has a higher intensity than the preceding plotted text (which it subsequently removes)
		_xcoor = mz_intensity[0]
		_ycoor = mz_intensity[1] / 10.0


		if len(whatToPlot) > 0:
			last_xcoor = whatToPlot[-1][0]
			last_ycoor = whatToPlot[-1][1]

			#print "%s --> %s" % (_xcoor, _xcoor - last_xcoor)

			if (_xcoor - last_xcoor < theRange * 0.05) & (_ycoor > last_ycoor) & (_ycoor > 1):  # if it is NOT far enough away BUT it is taller than the last text then go ahead
				whatToPlot.pop()
				whatToPlot.append([_xcoor, _ycoor])
				#print "popped! %s" % (last_xcoor)
			elif (_xcoor - last_xcoor > theRange * 0.05) & (_ycoor > 1):  # if it is far enough away from the last text then go ahead
				goAhead = True
				for checkPrevious in theRecords:
					if (_xcoor - checkPrevious[0] < theRange * 0.05) & (_ycoor < checkPrevious[1]):  # check to see if anything in range that is also greater than the current intensity
						goAhead = False
				if goAhead is True: whatToPlot.append([_xcoor, _ycoor])
		else:
			if (_ycoor > 1): whatToPlot.append([_xcoor, _ycoor])

		theRecords.append([_xcoor, _ycoor])

	#a.arrow(float(precursorMZ), 3, 0.0, -1, fc="k", ec="k", head_width=theRange*0.05, head_length=5)
	a.arrow(float(precursorMZ), 3, 0.0, -1, fc="k", ec="k", head_width=theRange * 0.005, head_length=2)

	for dataToPlot in whatToPlot:
		a.annotate(str(dataToPlot[0]),xy=(dataToPlot[0],dataToPlot[1]),ha='center',size='x-small')

	canvas.draw()

	# ---------------- mspepsearch processing---------------------------------
	#print("Running MSpepsearch...")
	inputData = open('inputFile.' + mspepsearchInputType, 'w')
	inputData.write(inputAllText)
	inputData.close()

	if libType == "hiRes":
		resParameters = ["/MPPM", "50"]
	elif libType == "loRes":
		resParameters = ["/ZPPM", "20", "/M", "0.4"]
	else:
		print("LIBTYPE IMPROPERLY SET!!!!!!") # this is really only for diagnostics

	theCall = ["MSPepSearch64.exe", "ydiPva", 
		"/ZI", "2.5"] + resParameters + ["/MPPM", "50", 
		"/MzLimits", "50", "-1", 
		"/HITS", "100", 
		"/OutBestHitsOnly"]

	for aLib in libPath:
		if "decoy" not in aLib.lower():
			theCall.extend(["/LIB", aLib])

	theCall.extend(["/LibInMem","/HiPri", "/INP",
		"inputFile." + mspepsearchInputType, # this is the input file
		"/OUTTAB", "outfile.tsv", "/OutSpecNum", "/OutPrecursorMZ",
		"/OutNumComp", "/OutNumMP","/OutMaxScore", "/OutPepLocalization",
		"/OutCE", "/MinMF", str(FDRbasedMFmin)]) # the MinMF should be set by FDR calculation

	subprocess.call(theCall)

	theResults = open("outfile.tsv", "r")

	processedResults = dict()
	lineCount = 0
	for line in theResults:
		parsedstuff = line.split("\t")

		if (lineCount >= 4) & (len(parsedstuff) > 5):


			theSeq = parsedstuff[19]
			theMod = parsedstuff[23] if parsedstuff[23] != "" else "0"
			theCharge = parsedstuff[20]
			theCollision = parsedstuff[21]

			uniqueID = theSeq + "_" + theMod + "_" + theCharge + "_" + theCollision

			MFscore = parsedstuff[9] # this needs to be double checked
			dotProduct = parsedstuff[10]
			deltaMass = parsedstuff[17]

			finalOutput = [theSeq, theMod, theCharge, theCollision, MFscore, dotProduct, deltaMass]

			processedResults[uniqueID] = finalOutput

		lineCount += 1

	#print(processedResults)
	left_tree.delete(*left_tree.get_children())
	for aName in processedResults:
		aResult = processedResults[aName]
		theSeq = aResult[0]
		theCharge = aResult[2]
		MFscore = aResult[4]
		deltaMass = aResult[6]
		#print(NIST17data[theNISTNO])
		left_tree.insert("", "end", aName, text=aName, values=(MFscore, deltaMass))


def changeLib():
	global libPath, NISTlibdata, libType

	prelibPath = filedialog.askdirectory(title="Select NIST formatted spectral library").replace("/", '\\')
	print(prelibPath)

	theRoot, theDirs, theFiles = next(os.walk(prelibPath))

	if len(theDirs) > 0:
		libPath = []

		for aDir in theDirs:
			libPath.append(os.path.join(theRoot, aDir))
	else:
		libPath = [prelibPath]

	sqlName = ""
	for aFile in theFiles:
		if ".sqlite" in aFile:
			sqlName = os.path.join(theRoot, aFile)

	print(sqlName)
	if sqlName == "":
		print("Missing SQLITE file!")


	root.wm_title("NIST MSDC Proteomics Spectral Searcher [%s]" % (prelibPath))
	#sqlName = libPath.split("\\")[-1]
	NISTlibdata = SqliteDict(sqlName, flag='r')

	#firstKey = list(NISTlibdata)[0] # don't do this! It's slow!
	firstKey = next(iter(NISTlibdata)) # we just want to grab the first key in the dictionary to see whether it's low res or not
	libType = "loRes" if firstKey[-1] == "_" else "hiRes"
	print("Library is %s!" % libType)

if __name__ == '__main__':
	root = Tk.Tk()  # Makes the window
	root.wm_title("NIST MSDC Proteomics Spectral Searcher")  # Makes the title that will appear in the top left
	root.geometry("1650x900")  # this is probably the secret to preventing the jumpy stuff for the navigation bar!!!!
	root.config(background="#FFFFFF")  # sets background color to white

	menubar = Tk.Menu(root)

	filemenu = Tk.Menu(menubar, tearoff=0)
	filemenu.add_command(label="Open MSP spectrum file", command=openMSPfile)
	filemenu.add_command(label="Process Thermo RAW file", command=rawExtract)
	filemenu.add_command(label="Change spectral library", command=changeLib)
	filemenu.add_separator()
	filemenu.add_command(label="Exit", command=root.quit)
	menubar.add_cascade(label="File", menu=filemenu)

	root.config(menu=menubar)

	rightFrame = Tk.Frame(root, width=1000, height=900)
	rightFrame.pack(side="right", fill="both", expand=True)
	#alsoRightFrame = Frame(root, width=900, height=900)
	#alsoRightFrame.pack(side="bottom",fill="x", expand=False)

	#rightFrame.grid(row=0, column=1, padx=10, pady=2)


	canvas = FigureCanvasTkAgg(f, master=rightFrame)

	canvas.get_tk_widget().pack(side="top", fill="both")

	theTool = NavigationToolbar2Tk(canvas, rightFrame)
	theTool.pack(side=Tk.TOP)
	theTool.update()

	canvas.draw()

	specData = Tk.Text(rightFrame, width=80, height=20, takefocus=0)
	#specData.grid(row=4, column=0, padx=10, pady=2)
	specData.pack(side="left", in_=rightFrame,fill="both", expand=True)

	specDataR = Tk.Text(rightFrame, width=80, height=20, takefocus=0)
	specDataR.pack(side="right", in_=rightFrame, fill="both", expand=True)


	# --------------------- top tree

	# put widgets here
	topleftFrame = Tk.Frame(root, width=100, height=1200)  # this is ABSOLUTELY the secret to preventing a jumpy window!!! expand MUST BE FALSE and the fill MUST BE "y"!!!
	topleftFrame.pack(side="top", fill="y", expand=True)

	right_tree = ttk.Treeview(topleftFrame, height="25")

	ysb1 = ttk.Scrollbar(orient=Tk.VERTICAL, command=right_tree.yview)
	right_tree['yscroll'] = ysb1.set


	ysb1.pack(side="right", in_=topleftFrame,fill="y")

	right_tree.heading("#0", text="Query")
	right_tree.column('#0', stretch="yes", minwidth=120, width=160)
	right_tree.pack(side="left", fill="both", expand=True)


	right_tree.bind('<<TreeviewSelect>>', itemClickedQueryTree)

	right_tree["columns"] = ("prec", "something")
	right_tree.column("prec", width=90)
	right_tree.column("something", width=50)
	right_tree.heading("prec", text="Precursor m/z")
	right_tree.heading("something", text="something")

	#-------------------- bottom treee

	leftbottomFrame = Tk.Frame(root, width=100, height = 600) # this is ABSOLUTELY the secret to preventing a jumpy window!!! expand MUST BE FALSE and the fill MUST BE "y"!!!
	leftbottomFrame.pack(side="bottom", fill="y", expand=True)

	left_tree = ttk.Treeview(leftbottomFrame, height="25")

	ysb2 = ttk.Scrollbar(orient=Tk.VERTICAL, command=left_tree.yview)
	left_tree['yscroll'] = ysb2.set


	ysb2.pack(side="right", in_=leftbottomFrame, fill="y")

	left_tree.heading("#0", text="Hits")
	left_tree.column('#0', stretch="yes", width=200)
	left_tree.pack(side="left", fill="both", expand=True)


	left_tree.bind('<<TreeviewSelect>>', itemClickedHitsTree)

	left_tree["columns"] = ("mf", "delta")
	left_tree.column("mf", width=40)
	left_tree.column("delta", width=60)
	left_tree.heading("mf", text="MF")
	left_tree.heading("delta", text="Delta m/z")

	#theSeparator = os.sep
	changeLib()

	#root.wm_title("NIST MSDC Proteomics Spectral Searcher [%s]" % (libPath))

	root.mainloop()  # start monitoring and updating the GUI. Nothing below here runs.

	print("Complete!")

