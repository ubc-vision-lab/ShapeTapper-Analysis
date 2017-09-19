# Creates a results list which matches the spec of what's used for the test code for medial axis computation.

#config file indices
config_indices = { #... here we go
    "block":0,
    "trialno":1,
    "feedback":2,
    #ignore 3,4,5? should check but don't need it at the moment. trial.js should have it
    "response_time":6,
    "tooSlow_image":7,
    "dyn_mask1":8,
    "dyn_mask2":9,
    "dyn_mask3":10,
    "dyn_mask4":11,
    "dyn_mask5":12,
    "repeat_trial":13,
    "e1_xpos":14, #x position as percentage of screen, from bottom left, for image 1, and for all images if not oddball
    "e1_ypos":15, #y position as percentage of scree, from bottom left, for image 1, and for all images if not oddball
    "onset":16,
    #here follows the image1 stuff
    "e1_image":17,
    "e1_rotation":18,
    "e1_safety":19, #this is actually the percentage of the screen the image's diagonal will take (scaling factor, relative to screen height)
    "e1_target":20,
    "e1_dyn_mask_flag":21,
    "e1_display_time":22,
    "e1_mask_display_time":23,
    "e1_time_before_e2":25,
    #image2
    "e2_image":25,
    "e2_rotation":26,
    "e2_safety":27, #this is actually the percentage of the screen the image's diagonal will take (scaling factor, relative to screen height)
    "e2_target":28,
    "e2_dyn_mask_flag":29,
    "e2_display_time":30,
    "e2_mask_display_time":31,
    "e2_time_before_e3":32,
    #image3
    "e3_image":33,
    "e3_rotation":34,
    "e3_safety":35, #this is actually the percentage of the screen the image's diagonal will take (scaling factor, relative to screen height)
    "e3_target":36,
    "e3_dyn_mask_flag":37,
    "e3_display_time":38,
    "e3_mask_display_time":39,
    #oddball stuff
    "oddball_flag":40,
    "e2_xpos":41,
    "e2_ypos":42,
    "e3_xpos":43,
    "e3_ypos":44
}

demographic_indices = {
    "uID":0,
    "Gender":1,
    "Handedness":2,
    "Age":3,
    "ScreenWidth":4,
    "ScreenHeight":5,
    "DPI":6,
    "configFile":7
}

trial_result_indices = {
    "trial_num":0,
    "badflag":1,
    "reaction_time":2,
    "touch_x":3,
    "touch_y":4
}

def getSubjectInfo(filepath):
    demo_file = open(filepath, 'r')
    return demo_file.readlines()[1].split(',') #hard coded second line
    

# given a valid config file, returns the list of trials.
def getConfigLines(filepath):
    #open the file
    f = open(filepath,'r')
    configLines = [line.split(',') for line in f.readlines()]
    return configLines, getBlockSizes(configLines)
    #split by line
    #return split

#given the set of trials, gives the size of every block
def getBlockSizes(configLines):
    previous_block = 1 # start at 1 cause never block 0
    previous_trial_number = 1
    trials_in_block = []
    
    for line in configLines:
        trial_number = int(line[config_indices['trialno']])
        curr_block = int(line[config_indices['block']])
        if previous_block != curr_block:
            trials_in_block.append(previous_trial_number)
            previous_block = curr_block
        previous_trial_number = trial_number
        
    trials_in_block.append(trial_number)
    return trials_in_block

# given a block and trial number, grabs the appropriate trial from the experiment configuration
def getTrial(block, trial_number, block_sizes, exp_config):
    index = sum(block_sizes[:block]) + trial_number - 1 #minus one because list is zero indexed
    return exp_config[index]

# given a trial, returns the target image name
def getImageName(trial):
    return trial[config_indices['e1_image']]

# given a trial, returns the rotation of the target image, in degrees
def getRotation(trial):
    return float(trial[config_indices['e1_rotation']])

#given the data from the shapetapper experiment, returns the list of valid data (badflag == 0)
#essentially, it's getConfigLines where we skip the first row.
def getData(filepath):
    #open file directory
    #grab the block data
    #initiate array
    #filter out bad ones
    #add block number to front of result
    #add to array
    #return the result
    return null

#given a row of data, checks if the row was good
#maybe ignore this for now
def goodResult(data):
    return data[trial_result_indices['badflag']] == 0

if __name__ == '__main__':
    
    experiment, block_sizes = getConfigLines("solo_blake_8_64_1.txt")
    print(block_sizes)
    # subjectID = "" #argv?
    # #maybe make this eval pwd?
    # data_dir = "/" #some directory
    # config_dir = "/" #some directory
    
    # demographic_file = subjectID + "_demographic.txt"
    # subject_info = getSubjectInfo(demographic_file)
    # experiment, block_sizes = getConfigLines(config_dir + demographic_info[demographic_indices['configFile']])
    
    # results_data = getData(data_dir, subjectID)
    
    # #all data has been acquired. time to stitch images together
    

    # for data in subject_data:
        # fetchTrial(data[0], data[1],)
