import cv2
import time
import os
import numpy as np
import cantrips as can
import matplotlib.pyplot as plt

def video_to_frames(input_video_file, read_dir, save_dir = None, save_file_type = 'png', max_n_frames_to_convert = np.inf):
    """Function to extract frames from input video file
    and save them as separate frames in an output directory.
    Args:
        input_loc: Input video file.
        output_loc: Output directory to save the frames.
    Returns:
        None
    """
    if save_dir == None:
        save_dir = read_dir
    # Log the time
    time_start = time.time()
    # Start capturing the feed
    cap = cv2.VideoCapture(read_dir + input_video_file)
    # Find the number of frames
    video_length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)) - 1
    print ("Number of frames: ", video_length)
    if max_n_frames_to_convert < video_length:
        print ('WARNING: There are more frames (' + str(video_length) + ') than the specified maximum number of frames (' + str(max_n_frames_to_convert) + ').  We will only save the first ' + str(max_n_frames_to_convert) + 'frames')
    count = 0
    print ("Converting video..\n")
    # Start converting the video
    while cap.isOpened():
        print ('Hello 1? ')
        # Extract the frame
        ret, frame = cap.read()
        if not ret:
            print ('Hello 2? ')
            continue
        # Write the results back to output location.
        print ('Hello 3?')
        save_name =  input_video_file.split('.')[0] + '_' + "%#05d."  % (count+1) + save_file_type
        if save_file_type.lower() == 'fits':
            print ('frame[:, :, 0] = ' + str(frame[:, :, 0] ))
            bw_frame = np.array(frame[:, :, 0]) + np.array(frame[:, :, 1]) + np.array(frame[:, :, 1])
            bw_frame = np.flip(bw_frame, axis = 0)
            #bw_frame = bw_frame * 0.0 + 1
            #bw_frame = np.zeros((100, 100)) + 2
            #bw_frame = np.array([[bw_frame[i,j] * i + j for i in range(np.shape(bw_frame)[0])] for j in range(np.shape(bw_frame)[1])])
            can.saveDataToFitsFile(np.transpose(bw_frame), save_name, save_dir, header = 'default')
        else:
            print ('save_name = ' + str(save_name))
            cv2.imwrite(save_dir + save_name, frame)
        count = count + 1
        # If there are no more frames left
        if (count > (video_length-1) or count > max_n_frames_to_convert):

            # Log the time again
            time_end = time.time()
            # Release the feed
            cap.release()
            # Print stats
            print ("Done extracting frames.\n%d frames extracted" % count)
            print ("It took %d seconds forconversion." % (time_end-time_start))
            break

if __name__=="__main__":

    input_loc = '/path/to/video/00009.MTS'
    output_loc = '/path/to/output/frames/'
    video_to_frames(input_loc, output_loc)
