import subprocess, os, re, threading, time
from sqlitedict import SqliteDict
import tkinter as Tk
import tkinter.ttk as ttk
from tkinter import filedialog

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from rdkit import Chem
from rdkit.Chem import rdFMCS, Descriptors, rdinchi, rdDepictor
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D

print("Loading data...")
NIST17data = SqliteDict('./tandemNIST17_gui_rdkit.sqlite', flag='r')
#msgpack_file = open("tandemNIST17_gui.msgpack", 'rb')
#NIST17data = msgpack.unpack(msgpack_file, use_list=False, encoding='utf-8')
print("Complete!")

#print(list(NIST17data)[:10])

precPatt = re.compile(r'PrecursorMZ:\s(\d+\.\d+)')
spectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)')

processedResults = dict()
theOverPlot = []
theRange = 0

# this is the path for the mspepsearch library
libPath = os.path.join("C:/", "LIBS", "nist_msms")
zppm = 20

f = Figure(figsize=(10, 6), dpi=100, facecolor='white')
a = f.add_subplot(111)
f.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
a.stem([0], [0], markerfmt=' ')
a.spines['top'].set_visible(False)
a.spines['right'].set_visible(False)


def moltosvg(mol, SMARTSinfo, picSize=(400, 400)):
	'''
	mc = Chem.Mol(mol.ToBinary())
	try:
		Chem.Kekulize(mc)
	except:
		mc = Chem.Mol(mol.ToBinary())

	mc.UpdatePropertyCache(strict=False)
	Chem.SanitizeMol(mc, Chem.SanitizeFlags.SANITIZE_KEKULIZE,catchErrors=True)
	'''
	rdDepictor.Compute2DCoords(mol)

	drawer = rdMolDraw2D.MolDraw2DSVG(picSize[0], picSize[1])
	#drawer.SetFontSize(6)
	#drawer.DrawMolecule(mc)
	drawer.DrawMolecule(mol, highlightAtoms=mol.GetSubstructMatch(Chem.MolFromSmarts(SMARTSinfo)))
	drawer.FinishDrawing()
	svg = drawer.GetDrawingText()
	'''
	smartMol = Chem.MolFromSmarts(SMARTSinfo)
	Chem.Kekulize(smartMol)
	rdDepictor.Compute2DCoords(smartMol)
	drawer2 = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
	drawer2.DrawMolecule(smartMol)
	drawer2.FinishDrawing()
	svg2 = drawer2.GetDrawingText()
	'''
	#return svg.replace('svg:',''), svg2.replace('svg:','')
	return svg


def mcsCalc(allTheMs):
	res = rdFMCS.FindMCS(allTheMs, ringMatchesRingOnly=True, completeRingsOnly=True, bondCompare=rdFMCS.BondCompare.CompareOrderExact)

	SMARTSinfo = res.smartsString
	if (SMARTSinfo != "") & (res.canceled != True):
		for i,anM in enumerate(allTheMs):
			svgText = moltosvg(anM, SMARTSinfo)
			text_file = open(str(i)+'.svg', "w")
			text_file.write(svgText)
			text_file.close()


def doMCS():
	MCSwindow = Tk.Toplevel(root)
	MCSwindow.wm_title("MCS analysis")  # Makes the title that will appear in the top left
	MCSwindow.geometry("700x700")  # this is probably the secret to preventing the jumpy stuff for the navigation bar!!!!
	MCSwindow.config(background="#FFFFFF")  # sets background color to white

	l = Tk.Label(MCSwindow, text="MCS analysis started...")
	l.pack(side="top", fill="both", expand=True, padx=100, pady=100)

	#l.config(text = "Wait till I'm done...")

	l.update_idletasks()

	sortedHits = sorted([[processedResults[resultName][0][1], processedResults[resultName][0][4]] for resultName in processedResults], reverse=True)
	if len(sortedHits) > 5:
		truncatedSorted = sortedHits[:5]
	else:
		truncatedSorted = sortedHits

	allTheMs = []
	for aHit in truncatedSorted:
		InChIKey = aHit[1]

		rdkitM = NIST17data[InChIKey][1]

		allTheMs.append(rdkitM)

	thread = threading.Thread(target = mcsCalc, args=(allTheMs,))
	startTime = time.time()
	thread.start()


	while thread.is_alive():
		MCSwindow.update()
		elapsedTime = round(time.time() - startTime)
		l.config(text = "MCS analysis still going (%ss elapsed)..." % elapsedTime)
		time.sleep(0.01)

	l.config(text = "MCS analysis started...Complete!")



def openMSPfile():
	global processedResults, theRange, xData, yData, sortedSpectrum
	#print("Hello!")
	file_path = filedialog.askopenfilename()
	inputSpectrumReader = open(file_path, "r")
	inputSpectrum = []; inputAllText = ""

	_count = 0
	for line in inputSpectrumReader:
		spectrumData = re.match(spectrumPatt, line)
		if spectrumData is not None:
			inputSpectrum.append([float(spectrumData.group(1)), float(spectrumData.group(2))])
			#_textOut += theMZ + " " + str(round(theAb, 2))
		else:
			precData = re.match(precPatt, line)
			if precData is not None:
				precursorMZ = float(precData.group(1))
		inputAllText += line
		_count += 1

	specData.delete(1.0, Tk.END)
	specData.insert(0.0, inputAllText)

	inputSpectrumReader.close()

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
	subprocess.call(["MSPepSearch64.exe", "dyGailv", 
		"/ZI", "1.6", 
		"/ZPPM", str(zppm), 
		"/MPPM", "40", 
		"/MzLimits", "50", "-1", 
		"/DiffIK", "1", # excludes library match (1K means first segment of InChI)
		"/HITS", "100", 
		"/OutBestHitsOnly",
		"/LIB", libPath, "/LibInMem","/HiPri", "/INP",
		file_path, # this is the input file
		"/OUTTAB", "outfile.tsv", 
		"/OutNumComp", "/OutNumMP", "/OutMaxScore", 
		"/OutNumPk", "/OutSrchNumPk", "/OutInstrType", 
		"/OutCE", "/OutSrchCE", "/OutPrecursorMZ",
		"/OutPrecursorType", "/OutChemForm", "/OutSrchChemForm", "/OutIK", 
		"/OutSrchID", "/OutNISTrn", "/OutSrchNISTrn"])

	theResults = open("outfile.tsv", "r")

	processedResults = dict()
	lineCount = 0
	for line in theResults:
		parsedstuff = line.split("\t")

		if (lineCount >= 4) & (len(parsedstuff) > 5):
			theName = parsedstuff[24]
			theNISTNO = parsedstuff[33]
			InChIKey = parsedstuff[34]
			MFscore = parsedstuff[14]
			dotProduct = parsedstuff[15]
			deltaMass = parsedstuff[22]

			finalOutput = [theNISTNO, MFscore, dotProduct, deltaMass, InChIKey]
			if theName in processedResults:
				processedResults[theName].append(finalOutput)
			else:
				processedResults[theName] = [finalOutput]

		lineCount += 1

	#print(processedResults)
	left_tree.delete(*left_tree.get_children())
	for aName in processedResults:
		for aResult in processedResults[aName]:
			theNISTNO = aResult[0]
			MFscore = aResult[1]
			#print(NIST17data[theNISTNO])
			left_tree.insert("", "end", aName, text=aName, values=(MFscore, theNISTNO))


def itemClicked(event):
	selectedHit = str(left_tree.focus())

	global processedResults, theOverPlot, theRange, xData, yData, sortedSpectrum

	resultAllText = ""

	theNISTNO = processedResults[selectedHit][0][0]
	deltaMass = processedResults[selectedHit][0][3]
	deltaMZ = float(deltaMass)

	nistSpectrum = NIST17data[theNISTNO][16]
	exactMass = NIST17data[theNISTNO][1]
	formula = NIST17data[theNISTNO][2]
	InChIKey = NIST17data[theNISTNO][6]
	adductType = NIST17data[theNISTNO][7]
	nistPrecursorMZ = NIST17data[theNISTNO][8]
	SMILES = NIST17data[InChIKey][0]
	theM = NIST17data[InChIKey][1]

	resultAllText += "Name: %s\nInChIKey: %s\nExact mass: %s\nPrecursor m/z: %s\nFormula: %s\nAdduct type: %s\nDelta m/z: %s\n\n" % (selectedHit, InChIKey, exactMass, nistPrecursorMZ, formula, adductType, deltaMass)

	#print(nistSpectrum)

	matchSpectrum = []

	resultAllText += "m/z\tIntensity\n---------------\n"
	for peak in nistSpectrum:
		resultAllText += "%s\t%s\n" % (peak[0], peak[1])
		matchSpectrum.append([float(peak[0]), peak[1]])


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
		yData2.append(-1 * anAttr[1] / 10)

	theRange = (max(xData + xData2) + 10) - (min(xData + xData2) - 10)
	a.set_xlim(min(xData + xData2) - 10, max(xData + xData2) + 10)

	xDataShift = []; yDataShift = []
	for anAttr in shiftedMatches:
		xDataShift.append(anAttr[0])
		yDataShift.append(-1 * anAttr[1] / 10)

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
		_ycoor = mz_intensity[1] / 10.0


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



root = Tk.Tk()  # Makes the window
root.wm_title("NIST Spectral Searcher")  # Makes the title that will appear in the top left
root.geometry("1700x1000")  # this is probably the secret to preventing the jumpy stuff for the navigation bar!!!!
root.config(background="#FFFFFF")  # sets background color to white

menubar = Tk.Menu(root)

filemenu = Tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Open MSP spectrum file", command=openMSPfile)
filemenu.add_command(label="MCS", command=doMCS)
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
leftFrame = Tk.Frame(root, width=100, height=1200)  # this is ABSOLUTELY the secret to preventing a jumpy window!!! expand MUST BE FALSE and the fill MUST BE "y"!!!
leftFrame.pack(side="top", fill="y", expand=True)

left_tree = ttk.Treeview(leftFrame, height="25")

ysb = ttk.Scrollbar(orient=Tk.VERTICAL, command=left_tree.yview)
left_tree['yscroll'] = ysb.set


ysb.pack(side="right", in_=leftFrame, fill="y")

left_tree.heading("#0", text="Name")
left_tree.column('#0', stretch="yes", minwidth=120, width=160)
left_tree.pack(side="left", fill="both", expand=True)


left_tree.bind('<<TreeviewSelect>>', itemClicked)

left_tree["columns"] = ("mf", "NISTNO")
left_tree.column("mf", width=90)
left_tree.column("NISTNO", width=50)
left_tree.heading("mf", text="MatchFactor")
left_tree.heading("NISTNO", text="NISTNO")



root.mainloop()  # start monitoring and updating the GUI. Nothing below here runs.

print("Complete!")

