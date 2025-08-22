All the functions have descriptive names or comments. Here is a brief overview of all the files.

!! Important !!
Make sure all these files and folders are in the same folder, some are importing others


---HEDS - Folder
Needed to connect to SLM, multiple programs import this

---HG_modes.ipynb
* Important jupyter notebook, contains functions to do knife edge, create and display mode masks
* Very unorganised. Breifly describing steps to do knife edge or mode mask generation
	
	First Connect to SLM, navigate to the part of the notebook with "Main" as the heading
	Connect to SLM using the first code block, a preview window should open
	Now run all the cells below it so that functions to capture and analyse images and creation of phasemasks is done
	!! Important !! If you want to display any phasemask, generate from the function and use slm.showImageData(phase_mask_function())
	use a for loop to move the phasemask around, there are a few cells doing that
	The knife edge function takes orientation as an input to specify horizontal or vertical


---mode_decomp.py
* Contains all the functions created in HG_modes.ipynb in a more ordered way
* Has comments, you can directly look at it. Main purpose is to do knife edge and produce mode masks for HG modes that we can show using sl,.displayImageData()


---test.ipynb
* Imports mode_decomp.py we can access all the functions here
* make sure to run the import and then connect to slm and camera before using any other functions
* the function to find beam property stores the data collected for knife edge
* Preferred notebook to use

---rigol_ch3_reader.py
* contains functions to read data from rigol channel 3
* collets 3 data points of dc voltage and averages it and returns it
* use if trying to use photodiode. Import this in the notebook to use its functions, alternatively they can be copied to the notebook


---analysis.ipynb
* reads csv data collected by knife edge function
* produces image with both the gaussian plots perpendicular to each other
