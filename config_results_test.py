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
    "e4_ypos":44
}
# given a valid config file, returns the list of trials.
def getConfigLines(filepath):
    #open the file
    #split by line
    #return split

#given the set of trials, gives the size of the block
def getBlockSize(configLines):
    for line in configLines:


# given a block and trial number, grabs the appropriate trial from the experiment configuration
# assumes that all blocks have the same number of trials
def getTrial(block, trial, trials_per_config, exp_config):
    #zero indexed
    #index = (block-1)*trials_per_config + trial-1
    #return exp_config[index]

# given a trial, returns the target image name
def getImageName(trial):
    #return trial[index of image1]

# given a trial, returns the rotation of the target image, in degrees
def getRotation(trial):
    #return float(trial[index of image1 rotation])

#given the data from the shapetapper experiment, returns the list of valid data (badflag == 0)
#essentially, it's getConfigLines where we skip the first row.
def getData(filepath, subject_ID):
    #open file directory
    #grab the block data
    #initiate array
    #filter out bad ones
    #add block number to front of result
    #add to array
    #return the result
    return null

#given a row of data, checks if the row was good
def goodResult(data):
    #return trial[index of badflag]

if __name__ == '__main__':
    config_filename = "config.txt"
    subject_ID = "something" #subject

    #maybe make this eval pwd?
    data_dir = "" #some directory
    config_dir = "" #some directory

    subject_data = getData(data_dir, subject)

    for data in subject_data:
        fetchTrial(data[0], data[1],)
