clear cond;

j = 1;
cond(j).name = 'low intensity directional proj high cs';
cond(j).basedir = 'e:\Phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Ones\CS';
cond(j).calibration = {'gershow2010-11-15a', 'Kane2010-12-20'};
cond(j).whichcalibration = [1 1 1 1 2 2 2 2 2];
cond(j).esetname = 'csDir45mid';

j = j+1;
cond(j).name = 'low intensity directional proj high rh5rh6';
cond(j).basedir = 'e:\Phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Ones\Rh5Rh6';
cond(j).calibration = {'gershow2010-11-15a', 'Kane2010-12-20'};
cond(j).whichcalibration = [1 2 2 2 2 2];
cond(j).esetname = 'rh5rh6Dir45mid';

j = j+1;
cond(j).name = 'low intensity directional proj low cs';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\All Zeros';
cond(j).calibration = {'Kane2010-12-20', 'Kane2011-02-09'};
cond(j).whichcalibration = [1 1 2 2 2 2];
cond(j).esetname = 'csDir45low';

j = j+1;
cond(j).name = 'high intensity directional cs';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\CantonS';
cond(j).calibration = {'Liz2011-03-31'};
cond(j).whichcalibration = [1 1 1 1];
cond(j).esetname = 'wcsDir45high';


j = j+1;
cond(j).name = 'high intensity directional wcs';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\wCs';
cond(j).calibration = {'Liz2011-03-31'};
cond(j).whichcalibration = [1 1 1];
cond(j).esetname = 'wcsDir45high';


j = j+1;
cond(j).name = 'high intensity directional gmr-hid';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\GMR-hid';
cond(j).calibration = {'Liz2011-03-31'};
cond(j).whichcalibration = [1 1 1 1];
cond(j).esetname = 'gmrhidDir45high';

j = j+1;
cond(j).name = 'high intensity directional rh5rh6';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\High Intensity\Directional\ywrh5rh6';
cond(j).calibration = {'Liz2011-03-31'};
cond(j).whichcalibration = [1 1 1 1];
cond(j).esetname = 'rh5rh6Dir45high';

j = j+1;
cond(j).name = 'projector off control';
cond(j).basedir = 'E:\phototaxis\Extracted Phototaxis Data rdiff\Low Intensity\Directional Gradients\Projector Off Control';
cond(j).calibration = {'Kane2011-02-02'};
cond(j).whichcalibration = [1 1 1 1];
cond(j).esetname = 'csDirProfOff';
