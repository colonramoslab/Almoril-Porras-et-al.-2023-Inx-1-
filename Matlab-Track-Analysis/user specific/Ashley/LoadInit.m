function [init] = LoadInit ()
%function [init] = LoadInit ()
%
%This function initializes the variables for extracting the maggot info,
%extracting the maggot tracks, etc.


%file io params... these get set when run program
init.MatlabPath=cd;
%init.DataPath='C:\Users\Ashley\Documents\Data\09\July\1545\';
%init.FileName='1545_1.jgp


%params to zero image except for where plate is
init.plate_startx=350;
init.plate_starty=75;
init.plate_endx=2200;
init.plate_endy=1900;
init.cam_y=1944;
init.cam_x=2592;

%params to get maggots information from a thresholded image
init.min_mag_area=10; %area has to be at least 10 pixels to be found
init.max_mag_area=200; %area over which we look for two maggots
init.maxideal_mag_area=250;
init.minideal_mag_area=125;



%params for calculating track and stitching parameters
init.frame_rate=4;   %frames per second (Hz) for the data
init.sigma=1*init.frame_rate;
init.num_frames_fortrack=2*init.frame_rate; %has to last 1 seconds
init.microns_per_pixel=66; %calibration specific to camera 
init.sep_frame=10*init.frame_rate;   %num of frames that tracks can be separated for stitching
                                    %can lose track for 10 secs in collision
init.sep_dist_perframe=400/init.microns_per_pixel/init.frame_rate;    %num of pixels that tracks can be 
                                                                      %separated, ~1000 um/s is max vel
init.sep_dist=init.sep_frame*init.sep_dist_perframe;
init.sep_dist_leeway=init.sep_dist_perframe+5;  %if collision occurs then track can be off by that much


%params for calculating run parameters
init.orientation=[1 0 0]; %vector that describes what direction = 0 degrees
init.angle_thresh=40/360*2*pi; %in radians, angle threshold for a headsweep
init.angle_back=20/360*2*pi;    %in radians, angle min to get out of headsweep
init.vel_thresh=120/init.frame_rate/init.microns_per_pixel;   %in pixels/frame, 120 microns/sec 
init.pts_to_calc=3;   %number of frames over which to calculate run parameters
init.num_frames_forrun=round(500/init.microns_per_pixel/init.vel_thresh); %worm = 1000 microns 
                                                                   %so has to move at least half a body length
                                                                   %has to be >1
        if(init.num_frames_forrun<2)
            disp('ERROR: Frame rate not high enough to get enough frames per run.');
        end


        
        
end