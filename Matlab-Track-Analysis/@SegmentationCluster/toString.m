function str = toString(sc)
%converts segmentation cluster or group of segmentation clusters to 
%a string representation that can be saved to disk
%function str = toString(sc)

import('org.yaml.snakeyaml.Yaml');
yaml = Yaml();
str = yaml.dump(sc.toHashMap);
