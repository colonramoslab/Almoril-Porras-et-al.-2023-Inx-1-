%generates documentation
%must run this script from matlab track analysis base directory
existsAndDefault('makeFrame', true); %whether to use the frame format for documentation
addpath([pwd '\m2html']);
[~,thisdir] = fileparts(pwd);
doc_directory = [thisdir filesep 'documentation'];
cd ..
rmdir(doc_directory,'s');
if (makeFrame)
    marc_m2html('mfiles',thisdir, 'htmldir', doc_directory, 'recursive','on', 'global','on',... 
          'template','frame', 'index','menu', 'graph','on');
else
    marc_m2html('mfiles',thisdir, 'htmldir', doc_directory, 'recursive','on', 'global','on',... 
          'graph','on');
end
cd (thisdir);
      
