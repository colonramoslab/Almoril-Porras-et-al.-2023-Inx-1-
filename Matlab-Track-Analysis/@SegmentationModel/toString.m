function str = toString(sm)
%converts segmentation model or group of segmentation models to 
%a string representation that can be saved to disk
%function str = toString(sc)

import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
str = yaml.dump(sm.toHashMap);
